from collections import defaultdict

from events import Events

from VMSfunctions.Chromatograms import *
from VMSfunctions.Common import *
import pandas as pd


# controller sends scan request
# mass spec generates scans (is an iterator over scans)
# scan contains: mz list, intensity list, rt, ms_level, precursor_mass, window
# simplest controller: just generates ms1 data
from VMSfunctions.DataGenerator import Peak


class Scan(object):
    def __init__(self, scan_id, mzs, intensities, ms_level, rt, density=None, isolation_windows = None):
        assert len(mzs) == len(intensities)
        self.scan_id = scan_id

        # ensure that mzs and intensites are sorted by their mz values
        p = mzs.argsort()
        self.mzs = mzs[p]
        self.intensities = intensities[p]

        self.ms_level = ms_level
        self.rt = rt
        self.num_peaks = len(mzs)

        self.density = density
        self.scan_duration = self._get_scan_duration()
        self.isolation_windows = isolation_windows

    def filter_intensity(self, min_intensity):
        keep = self.intensities > min_intensity
        self.mzs = self.mzs[keep]
        self.intensities = self.intensities[keep]
        self.num_peaks = len(self.mzs)

    def _get_scan_duration(self):
        if self.density is not None:
            scan_duration = self.density.scan_durations(self.ms_level, 1).flatten()
            return scan_duration[0]
        else: # TODO: improve this
            return 1

    def __repr__(self):
        return 'Scan %d num_peaks=%d rt=%.2f ms_level=%d' % (self.scan_id, self.num_peaks, self.rt, self.ms_level)


class ScanParameters(object):
    MS_LEVEL = 'ms_level'
    ISOLATION_WINDOWS = 'isolation_windows'
    PRECURSOR = 'precursor'

    def __init__(self):
        self.params = {}

    def set(self, key, value):
        self.params[key] = value

    def get(self, key):
        if key in self.params:
            return self.params[key]
        else:
            return None

    def __repr__(self):
        return 'ScanParameters %s' % (self.params)


class FragmentationEvent(object): # for benchmarking purpose
    def __init__(self, chem, query_rt, ms_level, peaks, scan_id):
        self.chem = chem
        self.query_rt = query_rt
        self.ms_level = ms_level
        self.peaks = peaks
        self.scan_id = scan_id

    def __repr__(self):
        return 'MS%d FragmentationEvent for %s at %f' % (self.ms_level, self.chem, self.query_rt)


class MassSpectrometer(LoggerMixin):
    MS_SCAN_ARRIVED = 'MsScanArrived'
    ACQUISITION_STREAM_OPENING = 'AcquisitionStreamOpening'
    ACQUISITION_STREAM_CLOSING = 'AcquisitionStreamClosing'

    def __init__(self, ionisation_mode):
        self.ionisation_mode = ionisation_mode

        # following IAPI events
        self.events =  Events((self.MS_SCAN_ARRIVED, self.ACQUISITION_STREAM_OPENING, self.ACQUISITION_STREAM_CLOSING,))
        self.event_dict = {
            self.MS_SCAN_ARRIVED: self.events.MsScanArrived,
            self.ACQUISITION_STREAM_OPENING: self.events.AcquisitionStreamOpening,
            self.ACQUISITION_STREAM_CLOSING: self.events.AcquisitionStreamClosing
        }

    def get_next_scan(self):
        raise NotImplementedError()

    def fire_event(self, event_name, arg=None):
        if event_name not in self.event_dict:
            raise ValueError('Unknown event name')

        # fire the event
        e = self.event_dict[event_name]
        if arg is not None:
            e(arg)
        else:
            e()

    def register(self, event_name, handler):
        if event_name not in self.event_dict:
            raise ValueError('Unknown event name')
        e = self.event_dict[event_name]
        e += handler # register a new event handler for e

    def clear(self, event_name):
        if event_name not in self.event_dict:
            raise ValueError('Unknown event name')
        e = self.event_dict[event_name]
        e.targets = []


# Independent here refers to how the intensity of each peak in a scan is independent of each other
# i.e. there's no ion supression effect
class IndependentMassSpectrometer(MassSpectrometer):
    def __init__(self, ionisation_mode, chemicals, density=None):
        super().__init__(ionisation_mode)
        self.chemicals = chemicals
        self.idx = 0
        self.time = 0
        self.queue = []
        self.repeating_scan_parameters = None
        self.precursor_information = defaultdict(list) # key: Precursor object, value: ms2 scans
        self.density = density # a PeakDensityEstimator object
        self.fragmentation_events = [] # which chemicals produce which peaks

        try: # for EmpiricalChromatogram, we have the min_rt and max_rt
            self.chrom_min_rts = np.array([chem.chrom.min_rt for chem in self.chemicals])
            self.chrom_max_rts = np.array([chem.chrom.max_rt for chem in self.chemicals])
        except AttributeError: # TODO: ask vinny how to set min_rt and max_rt for functional chromatograms
            self.chrom_min_rts = np.array([chem.rt-np.inf for chem in self.chemicals])
            self.chrom_max_rts = np.array([chem.rt+np.inf for chem in self.chemicals])

    def run(self, min_time, max_time, pbar=None):
        self.time = min_time
        self.fire_event(MassSpectrometer.ACQUISITION_STREAM_OPENING)
        try:
            while self.time < max_time:
                initial_time = self.time

                # if the processing queue is empty, then just do the repeating scan
                if len(self.queue) == 0:
                    param = self.repeating_scan_parameters
                else:
                    # otherwise pop the parameter for the next scan from the queue
                    param = self.queue.pop(0)
                scan = self.get_next_scan(param)
                # logger.debug('time %f scan %s' % (self.time, scan))

                # if MS2 and above, and the controller tells us which precursor ion the scan is coming from, store it
                precursor = param.get(ScanParameters.PRECURSOR)
                if scan.ms_level >= 2 and precursor is not None and scan.num_peaks > 0:
                    self.precursor_information[precursor].append(scan)

                # print a progress bar if provided
                if pbar is not None:
                    elapsed = self.time - initial_time
                    pbar.update(elapsed)

        finally:
            self.fire_event(MassSpectrometer.ACQUISITION_STREAM_CLOSING)
            if pbar is not None:
                pbar.close()

    def get_next_scan(self, param):
        if param is not None:
            scan = self._get_scan(self.time, param)
            self.fire_event(self.MS_SCAN_ARRIVED, scan)
            self.idx += 1
            self.time += scan.scan_duration
            return scan
        else:
            return None

    def add_to_queue(self, param):
        self.queue.append(param)

    def disable_repeating_scan(self):
        self.set_repeating_scan(None)

    def set_repeating_scan(self, params):
        self.repeating_scan_parameters = params

    def reset(self):
        for key in self.event_dict: # clear event handlers
            self.clear(key)
        self.time = 0 # reset internal time and index to 0
        self.idx = 0

    def _get_scan(self, scan_time, param):
        """
        Constructs a scan at a particular timepoint
        :param time: the timepoint
        :return: a mass spectrometry scan at that time
        """
        scan_mzs = []  # all the mzs values in this scan
        scan_intensities = []  # all the intensity values in this scan
        ms_level = param.get(ScanParameters.MS_LEVEL)
        isolation_windows = param.get(ScanParameters.ISOLATION_WINDOWS)
        scan_id = self.idx

        # for all chemicals that come out from the column coupled to the mass spec
        for i in self._get_chem_indices(scan_time):
            chemical = self.chemicals[i]

            # mzs is a list of (mz, intensity) for the different adduct/isotopes combinations of a chemical            
            mzs = self._get_all_mz_peaks(chemical, scan_time, ms_level, isolation_windows)
            peaks = []
            if mzs is not None:
                chem_mzs = []
                chem_intensities = []
                for peak_mz, peak_intensity in mzs:
                    if peak_intensity > 0:
                        chem_mzs.append(peak_mz)
                        chem_intensities.append(peak_intensity)
                        p = Peak(peak_mz, scan_time, peak_intensity, ms_level)
                        peaks.append(p)

                scan_mzs.extend(chem_mzs)
                scan_intensities.extend(chem_intensities)

            # for benchmarking purpose
            if len(peaks) > 0:
                frag = FragmentationEvent(chemical, scan_time, ms_level, peaks, scan_id)
                self.fragmentation_events.append(frag)

        scan_mzs = np.array(scan_mzs)
        scan_intensities = np.array(scan_intensities)
        return Scan(scan_id, scan_mzs, scan_intensities, ms_level, scan_time,
                    density=self.density, isolation_windows=isolation_windows)

    def _get_chem_indices(self, query_rt):
        rtmin_check = self.chrom_min_rts <= query_rt
        rtmax_check = query_rt <= self.chrom_max_rts
        idx = np.nonzero(rtmin_check & rtmax_check)[0]
        return idx

    def _get_all_mz_peaks(self, chemical, query_rt, ms_level, isolation_windows):
        if not self._rt_match(chemical, query_rt):
            return None
        mz_peaks = []
        for which_isotope in range(len(chemical.isotopes)):
            for which_adduct in range(len(self._get_adducts(chemical))):
                mz_peaks.extend(self._get_mz_peaks(chemical, query_rt, ms_level, isolation_windows, which_isotope, which_adduct))       
        if mz_peaks == []:
            return None
        else:
            return mz_peaks

    def _get_mz_peaks(self, chemical, query_rt, ms_level, isolation_windows, which_isotope, which_adduct):
        # Vinny says:
        # so the first "if" return the ms1 peaks if it is an ms1 scan
        # the first "elif" stops it returning any ms2+ fragments when its not monoisotopic
        # the second "elif" returns ms2 fragments if ms2 scan, ms3 fragments if ms3 scan
        # etc etc
        # the "else" checks isolation window matches
        # for example. if we wants to do an ms2 scan on a chemical. we would first have ms_level=2 and the chemicals
        # ms_level =1. So we would go to the "else". We then check the ms1 window matched. It then would loop through
        # the children who have ms_level = 2. So we then go to second elif and return the mz and intensity of each ms2 fragment
        # so in the last "else". if we check the isolation window match, then we query the children of that chemical
        # (which will have a different ms2 level). else we will return no fragments because the window isnt right
        mz_peaks = []
        if ms_level == 1 and chemical.ms_level == 1: # fragment ms1 peaks
            if not (which_isotope > 0 and which_adduct > 0):
                if self._isolation_match(chemical, query_rt, isolation_windows[0], which_isotope, which_adduct):
                    intensity = self._get_intensity(chemical, query_rt, which_isotope, which_adduct)
                    mz = self._get_mz(chemical, query_rt, which_isotope, which_adduct)
                    mz_peaks.extend([(mz, intensity)])
        elif ms_level > 1 and which_isotope > 0:
            pass # TODO: we need to deal with msN fragmentation of non-monoisotopic peaks?
        elif ms_level == chemical.ms_level: # fragment MSn peaks
            intensity = self._get_intensity(chemical, query_rt, which_isotope, which_adduct)
            mz = self._get_mz(chemical, query_rt, which_isotope, which_adduct)
            return [(mz, intensity)]
        else:
            if self._isolation_match(chemical, query_rt, isolation_windows[chemical.ms_level - 1], which_isotope,
                                     which_adduct) and chemical.children is not None:
                for i in range(len(chemical.children)):
                    mz_peaks.extend(self._get_mz_peaks(chemical.children[i], query_rt, ms_level, isolation_windows,
                                                       which_isotope, which_adduct))
            else:
                return []
        return mz_peaks

    def _get_adducts(self, chemical):
        if chemical.ms_level == 1:
            return chemical.adducts
        else:
            return self._get_adducts(chemical.parent)

    def _rt_match(self, chemical, query_rt):
        return chemical.ms_level != 1 or chemical.chromatogram._rt_match(query_rt - chemical.rt)

    def _get_intensity(self, chemical, query_rt, which_isotope, which_adduct):
        if chemical.ms_level == 1:
            intensity = chemical.isotopes[which_isotope][1] * self._get_adducts(chemical)[which_adduct][1] * \
                        chemical.max_intensity
            return intensity * chemical.chromatogram.get_relative_intensity(query_rt - chemical.rt)
        else:               
            return self._get_intensity(chemical.parent, query_rt, which_isotope, which_adduct) * \
                   chemical.parent_mass_prop * chemical.prop_ms2_mass

    def _get_mz(self, chemical, query_rt, which_isotope, which_adduct):
        if chemical.ms_level == 1:
            return (adduct_transformation(chemical.isotopes[which_isotope][0],
                                          self._get_adducts(chemical)[which_adduct][0]) +
                    chemical.chromatogram.get_relative_mz(query_rt - chemical.rt))
        else:
            return adduct_transformation(chemical.isotopes[which_isotope][0], self._get_adducts(chemical)[which_adduct][0])

    def _isolation_match(self, chemical, query_rt, isolation_windows, which_isotope, which_adduct):
        # assumes list is formated like:
        # [(min_1,max_1),(min_2,max_2),...]
        for window in isolation_windows:
            if window[0] < self._get_mz(chemical, query_rt, which_isotope, which_adduct) <= window[1]:
                return True
        return False

class DsDAMassSpec(IndependentMassSpectrometer):
    
    def run(self, schedule, pbar=None):
        self.schedule = schedule
        self.time = schedule["targetTime"][0]
        initial_time = schedule["targetTime"][0]
        self.fire_event(MassSpectrometer.ACQUISITION_STREAM_OPENING)
        try:
            for time_index in range(len(self.schedule["targetTime"])):
                initial_time = self.time
                if len(self.queue) == 0:
                    param = self.repeating_scan_parameters
                else:
                    # otherwise pop the parameter for the next scan from the queue
                    param = self.queue.pop(0)
                scan = self.get_next_scan(param, time_index)
                # logger.debug('time %f scan %s' % (self.time, scan))

                # if MS2 and above, and the controller tells us which precursor ion the scan is coming from, store it
                precursor = param.get(ScanParameters.PRECURSOR)
                if scan.ms_level >= 2 and precursor is not None and scan.num_peaks > 0:
                    self.precursor_information[precursor].append(scan)

                # print a progress bar if provided
                if pbar is not None and self.time is not None:
                    elapsed = self.time - initial_time
                    pbar.update(elapsed)
                    
        finally:
            self.fire_event(MassSpectrometer.ACQUISITION_STREAM_CLOSING)
            if pbar is not None:
                pbar.close()

    def get_next_scan(self, param, time_index):
        if param is not None:
            scan = self._get_scan(self.time, param)
            self.fire_event(self.MS_SCAN_ARRIVED, scan)
            self.idx += 1
            if time_index + 1 < len(self.schedule["targetTime"]):
                self.time = self.schedule["targetTime"][time_index + 1]
            else:
                self.time = None
            return scan
        else:
            return None

# class ThermoFusionMassSpectrometer:

#     def __next__(self):
#         raise NotImplementedError()