from events import Events

from VMSfunctions.Chemicals import *
from VMSfunctions.Chromatograms import *


# controller sends scan request
# mass spec generates scans (is an iterator over scans)
# scan contains: mz list, intensity list, rt, ms_level, precursor_mass, window
# simplest controller: just generates ms1 data

class Scan(object):
    def __init__(self, scan_id, mzs, intensities, ms_level, rt):
        assert len(mzs) == len(intensities)
        self.scan_id = scan_id

        # ensure that mzs and intensites are sorted by their mz values
        p = mzs.argsort()
        self.mzs = mzs[p]
        self.intensities = intensities[p]

        self.ms_level = ms_level
        self.rt = rt
        self.num_peaks = len(mzs)
        self.scan_duration = self._get_scan_duration()

    def _get_scan_duration(self):
        # TODO: improve this
        return 1

    def __repr__(self):
        return 'Scan %d -- num_peaks=%d rt=%.2f ms_level=%d' % (self.scan_id, self.num_peaks, self.rt, self.ms_level)


class MassSpectrometer(object):
    MS_SCAN_ARRIVED = 'MsScanArrived'
    ACQUISITION_STREAM_OPENING = 'AcquisitionStreamOpening'
    ACQUISITION_STREAM_CLOSING = 'AcquisitionStreamClosing'

    def __init__(self):
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


# Independent here refers to how the intensity of each peak in a scan is independent of each other
# i.e. there's no ion supression effect
class IndependentMassSpectrometer(MassSpectrometer):
    def __init__(self, chemicals, scan_parameters):
        super().__init__()
        self.chemicals = chemicals
        self.idx = 0
        self.time = 0
        self.scan_parameters = scan_parameters

    def run(self, max_time):
        self.fire_event(MassSpectrometer.ACQUISITION_STREAM_OPENING)
        try:
            while self.time < max_time:
                scan = self.get_next_scan(self.scan_parameters)
                self.fire_event(self.MS_SCAN_ARRIVED, scan)
        finally:
            self.fire_event(MassSpectrometer.ACQUISITION_STREAM_CLOSING)

    def get_next_scan(self, scan_parameters):
        ms_level = scan_parameters['ms_level']
        isolation_windows = scan_parameters['isolation_windows']
        scan = self._get_scan(self.time, ms_level, isolation_windows)
        self.idx += 1
        self.time += scan.scan_duration
        return scan

    def _get_scan(self, scan_time, scan_level, isolation_windows):
        """
        Constructs a scan at a particular timepoint
        :param time: the timepoint
        :return: a mass spectrometry scan at that time
        """
        if scan_level > 1:
            raise NotImplementedError()  # TODO: add ms2 support
        scan_mzs = []  # all the mzs values in this scan
        scan_intensities = []  # all the intensity values in this scan

        # for all chemicals that come out from the column coupled to the mass spec
        for i in range(len(self.chemicals)):
            chemical = self.chemicals[i]

            # mzs is a list of (mz, intensity) for the different adduct/isotopes combinations of a chemical
            mzs = self._get_all_mz_peaks(chemical, scan_time, scan_level, isolation_windows)

            if mzs is not None:
                chem_mzs = [x[0] for x in mzs]
                chem_intensities = [x[1] for x in mzs]
                scan_mzs.extend(chem_mzs)
                scan_intensities.extend(chem_intensities)

        scan_mzs = np.array(scan_mzs)
        scan_intensities = np.array(scan_intensities)
        return Scan(self.idx, scan_mzs, scan_intensities, scan_level, scan_time)

    def _get_all_mz_peaks(self, chemical, query_rt, ms_level, isolation_windows):
        if ms_level == 1:
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
        mz_peaks = []
        if ms_level == 1 and chemical.ms_level == 1:
            if not (which_isotope > 0 and which_adduct > 0):
                if self._isolation_match(chemical, query_rt, isolation_windows[0], which_isotope, which_adduct):
                    intensity = self._get_intensity(chemical, query_rt, which_isotope, which_adduct)
                    mz = self._get_mz(chemical, query_rt, which_isotope, which_adduct)
                    mz_peaks.extend([(mz, intensity)])
        elif ms_level > 1 and which_isotope > 0:
            pass
        elif ms_level == chemical.ms_level:
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
        if chemical.ms_level == 1:
            if chemical.chromatogram._rt_match(query_rt - chemical.rt):
                return True
            else:
                return False
        else:
            True

    def _get_intensity(self, chemical, query_rt, which_isotope, which_adduct):
        if chemical.ms_level == 1:
            intensity = chemical.isotopes[which_isotope][1] * self._get_adducts(chemical)[which_adduct][1] * \
                        chemical.max_intensity
            return intensity * chemical.chromatogram.get_relative_intensity(query_rt - chemical.rt)
        else:
            return self._get_intensity(chemical.parent, query_rt, which_isotope, which_adduct) * \
                   chemical.parent_mass_prop

    def _get_mz(self, chemical, query_rt, which_isotope, which_adduct):
        if chemical.ms_level == 1:
            return (adductTransformation(chemical.isotopes[which_isotope][0],
                                         self._get_adducts(chemical)[which_adduct][0]) +
                    chemical.chromatogram.get_relative_mz(query_rt - chemical.rt))
        else:
            return adductTransformation(chemical.isotopes[which_isotope][0], self._get_adducts(chemical)[which_adduct][0])

    def _isolation_match(self, chemical, query_rt, isolation_windows, which_isotope, which_adduct):
        # assumes list is formated like:
        # [(min_1,max_1),(min_2,max_2),...]
        for window in isolation_windows:
            if window[0] < self._get_mz(chemical, query_rt, which_isotope, which_adduct) <= window[1]:
                return True
        return False

# class ThermoFusionMassSpectrometer:

#     def __next__(self):
#         raise NotImplementedError()