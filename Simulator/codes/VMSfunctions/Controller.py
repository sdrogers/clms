import sys
from collections import namedtuple

import pylab as plt
from tqdm import tqdm

from VMSfunctions.MassSpec import *
from VMSfunctions.MzmlWriter import *


class Controller(LoggerMixin):
    def __init__(self, mass_spec):
        self.scans = defaultdict(list)  # key: ms level, value: list of scans for that level
        self.mass_spec = mass_spec
        self.make_plot = False

    def handle_scan(self, scan):
        self.scans[scan.ms_level].append(scan)
        self._process_scan(scan)
        self._update_parameters(scan)

    def handle_acquisition_open(self):
        raise NotImplementedError()

    def handle_acquisition_closing(self):
        raise NotImplementedError()

    def write_mzML(self, analysis_name, outfile):
        writer = MzmlWriter(analysis_name, self.scans, precursor_information=self.mass_spec.precursor_information)
        writer.write_mzML(outfile)

    def _process_scan(self, scan):
        raise NotImplementedError()

    def _update_parameters(self, scan):
        raise NotImplementedError()

    def run(self, min_time, max_time, progress_bar=True):
        raise NotImplementedError()

    def _plot_scan(self, scan):
        if self.make_plot:
            plt.figure()
            for i in range(scan.num_peaks):
                x1 = scan.mzs[i]
                x2 = scan.mzs[i]
                y1 = 0
                y2 = scan.intensities[i]
                a = [[x1, y1], [x2, y2]]
                plt.plot(*zip(*a), marker='', color='r', ls='-', lw=1)
            plt.title('Scan {0} {1}s -- {2} peaks'.format(scan.scan_id, scan.rt, scan.num_peaks))
            plt.show()


class SimpleMs1Controller(Controller):
    def __init__(self, mass_spec):
        super().__init__(mass_spec)
        default_scan = ScanParameters()
        default_scan.set(ScanParameters.MS_LEVEL, 1)
        default_scan.set(ScanParameters.ISOLATION_WINDOWS, [[DEFAULT_MS1_SCAN_WINDOW]])

        mass_spec.reset()
        mass_spec.set_repeating_scan(default_scan)
        mass_spec.register(MassSpectrometer.MS_SCAN_ARRIVED, self.handle_scan)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_OPENING, self.handle_acquisition_open)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_CLOSING, self.handle_acquisition_closing)

    def run(self, min_time, max_time, progress_bar=True):
        if progress_bar:
            with tqdm(total=max_time - min_time, initial=0) as pbar:
                self.mass_spec.run(min_time, max_time, pbar=pbar)
        else:
            self.mass_spec.run(min_time, max_time)

    def handle_acquisition_open(self):
        self.logger.info('Acquisition open')

    def handle_acquisition_closing(self):
        self.logger.info('Acquisition closing')

    def _process_scan(self, scan):
        if scan.num_peaks > 0:
            self.logger.info('Time %f Received %s' % (self.mass_spec.time, scan))
            self._plot_scan(scan)

    def _update_parameters(self, scan):
        pass  # do nothing


class Precursor(object):
    def __init__(self, precursor_mz, precursor_intensity, precursor_charge, precursor_scan_id):
        self.precursor_mz = precursor_mz
        self.precursor_intensity = precursor_intensity
        self.precursor_charge = precursor_charge
        self.precursor_scan_id = precursor_scan_id

    def __str__(self):
        return 'Precursor mz %f intensity %f charge %d scan_id %d' % (self.precursor_mz, self.precursor_intensity, self.precursor_charge, self.precursor_scan_id)

class TopNController(Controller):
    def __init__(self, mass_spec, N, isolation_window, mz_tol, rt_tol, min_ms1_intensity):
        super().__init__(mass_spec)
        self.last_ms1_scan = None
        self.N = N
        self.isolation_window = isolation_window  # the isolation window (in Dalton) to select a precursor ion
        self.mz_tol = mz_tol  # the m/z window (ppm) to prevent the same precursor ion to be fragmented again
        self.rt_tol = rt_tol  # the rt window to prevent the same precursor ion to be fragmented again
        self.min_ms1_intensity = min_ms1_intensity  # minimum ms1 intensity to fragment

        mass_spec.reset()
        default_scan = ScanParameters()
        default_scan.set(ScanParameters.MS_LEVEL, 1)
        default_scan.set(ScanParameters.ISOLATION_WINDOWS, [[DEFAULT_MS1_SCAN_WINDOW]])
        mass_spec.set_repeating_scan(default_scan)

        # register new event handlers under this controller
        mass_spec.register(MassSpectrometer.MS_SCAN_ARRIVED, self.handle_scan)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_OPENING, self.handle_acquisition_open)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_CLOSING, self.handle_acquisition_closing)

    def run(self, min_time, max_time, progress_bar=True):
        if progress_bar:
            with tqdm(total=max_time - min_time, initial=0) as pbar:
                self.mass_spec.run(min_time, max_time, pbar=pbar)
        else:
            self.mass_spec.run(min_time, max_time)

    def handle_acquisition_open(self):
        self.logger.info('Time %f Acquisition open' % self.mass_spec.time)

    def handle_acquisition_closing(self):
        self.logger.info('Time %f Acquisition closing' % self.mass_spec.time)

    def _process_scan(self, scan):
        self.logger.info('Time %f Received from mass spec %s' % (self.mass_spec.time, scan))
        if scan.ms_level == 1:  # we get an ms1 scan, if it has a peak, then store it for fragmentation next time
            if scan.num_peaks > 0:
                self.last_ms1_scan = scan
            else:
                self.last_ms1_scan = None

        elif scan.ms_level == 2:  # if we get ms2 scan, then do something with it
            # scan.filter_intensity(self.min_ms2_intensity)
            if scan.num_peaks > 0:
                self._plot_scan(scan)

    def _update_parameters(self, scan):

        # if there's a previous ms1 scan to process
        if self.last_ms1_scan is not None:

            mzs = self.last_ms1_scan.mzs
            intensities = self.last_ms1_scan.intensities
            rt = self.last_ms1_scan.rt

            # loop over points in decreasing intensity
            fragmented_count = 0
            idx = np.argsort(intensities)[::-1]
            for i in idx:
                mz = mzs[i]
                intensity = intensities[i]

                # stopping criteria is after we've fragmented N ions or we found ion < min_intensity
                if fragmented_count >= self.N:
                    self.logger.debug('Time %f Top-%d ions have been selected' % (self.mass_spec.time, self.N))
                    break

                if intensity < self.min_ms1_intensity:
                    self.logger.debug(
                        'Time %f Minimum intensity threshold %f reached at %f, %d' % (self.mass_spec.time, self.min_ms1_intensity, intensity, fragmented_count))
                    break

                # skip ion in the dynamic exclusion list of the mass spec
                if self.mass_spec.exclude(mz, rt):
                    continue

                # send a new ms2 scan parameter to the mass spec
                dda_scan_params = ScanParameters()
                dda_scan_params.set(ScanParameters.MS_LEVEL, 2)

                # create precursor object, assume it's all singly charged
                precursor_charge = +1 if (self.mass_spec.ionisation_mode == POSITIVE) else -1
                precursor = Precursor(precursor_mz=mz, precursor_intensity=intensity,
                                      precursor_charge=precursor_charge, precursor_scan_id=self.last_ms1_scan.scan_id)
                mz_lower = mz - self.isolation_window # Da
                mz_upper = mz + self.isolation_window # Da
                isolation_windows = [[(mz_lower, mz_upper)]]
                dda_scan_params.set(ScanParameters.ISOLATION_WINDOWS, isolation_windows)
                dda_scan_params.set(ScanParameters.PRECURSOR, precursor)

                # save dynamic exclusion parameters too
                dda_scan_params.set(ScanParameters.DYNAMIC_EXCLUSION_MZ_TOL, self.mz_tol)
                dda_scan_params.set(ScanParameters.DYNAMIC_EXCLUSION_RT_TOL, self.rt_tol)

                # push this dda scan parameter to the mass spec queue
                self.mass_spec.add_to_queue(dda_scan_params)
                fragmented_count += 1

            for param in self.mass_spec.queue:
                precursor = param.get(ScanParameters.PRECURSOR)
                if precursor is not None:
                    self.logger.debug('- %s' % str(precursor))

                # set this ms1 scan as has been processed
            self.last_ms1_scan = None


class TreeController(Controller):
    def __init__(self, mass_spec, dia_design, window_type, kaufmann_design, extra_bins, num_windows=None):
        super().__init__(mass_spec)
        self.last_ms1_scan = None
        self.dia_design = dia_design
        self.window_type = window_type
        self.kaufmann_design = kaufmann_design
        self.extra_bins = extra_bins
        self.num_windows = num_windows

        mass_spec.reset()
        default_scan = ScanParameters()
        default_scan.set(ScanParameters.MS_LEVEL, 1)
        default_scan.set(ScanParameters.ISOLATION_WINDOWS, [[DEFAULT_MS1_SCAN_WINDOW]])
        mass_spec.set_repeating_scan(default_scan)

        mass_spec.register(MassSpectrometer.MS_SCAN_ARRIVED, self.handle_scan)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_OPENING, self.handle_acquisition_open)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_CLOSING, self.handle_acquisition_closing)

    def run(self, min_time, max_time, progress_bar=True):
        if progress_bar:
            with tqdm(total=max_time - min_time, initial=0) as pbar:
                self.mass_spec.run(min_time, max_time, pbar=pbar)
        else:
            self.mass_spec.run(min_time, max_time)

    def handle_acquisition_open(self):
        self.logger.info('Acquisition open')

    def handle_acquisition_closing(self):
        self.logger.info('Acquisition closing')

    def _process_scan(self, scan):
        self.logger.info('Received scan {}'.format(scan))
        if scan.ms_level == 1:  # if we get a non-empty ms1 scan
            if scan.num_peaks > 0:
                self.last_ms1_scan = scan
            else:
                self.last_ms1_scan = None

        elif scan.ms_level == 2:  # if we get ms2 scan, then do something with it
            if scan.num_peaks > 0:
                self._plot_scan(scan)

    def _update_parameters(self, scan):

        # if there's a previous ms1 scan to process
        if self.last_ms1_scan is not None:

            rt = self.last_ms1_scan.rt

            # then get the last ms1 scan, select bin walls and create scan locations
            mzs = self.last_ms1_scan.mzs
            default_range = [DEFAULT_MS1_SCAN_WINDOW]  # TODO: this should maybe come from somewhere else?
            locations = DiaWindows(mzs, default_range, self.dia_design, self.window_type, self.kaufmann_design,
                                   self.extra_bins, self.num_windows).locations
            self.logger.debug('Window locations {}'.format(locations))
            for i in range(len(locations)):  # define isolation window around the selected precursor ions
                isolation_windows = locations[i]
                dda_scan_params = ScanParameters()
                dda_scan_params.set(ScanParameters.MS_LEVEL, 2)
                dda_scan_params.set(ScanParameters.ISOLATION_WINDOWS, isolation_windows)
                self.mass_spec.add_to_queue(dda_scan_params)  # push this dda scan to the mass spec queue

            # set this ms1 scan as has been processed
            self.last_ms1_scan = None


class KaufmannWindows(object):
    """
    Method for creating window designs based on Kaufmann paper - https://www.ncbi.nlm.nih.gov/pubmed/27188447
    """

    def __init__(self, bin_walls, bin_walls_extra, kaufmann_design, extra_bins=0):
        self.locations = []
        if kaufmann_design == "nested":
            n_locations_internal = 4
            for i in range(0, 8):
                self.locations.append([[(bin_walls[(0 + i * 8)], bin_walls[(8 + i * 8)])]])
        elif kaufmann_design == "tree":
            n_locations_internal = 3
            self.locations.append([[(bin_walls[0], bin_walls[32])]])
            self.locations.append([[(bin_walls[32], bin_walls[64])]])
            self.locations.append([[(bin_walls[16], bin_walls[48])]])
            self.locations.append([[(bin_walls[8], bin_walls[24]), (bin_walls[40], bin_walls[56])]])
        else:
            raise ValueError("not a valid design")
        locations_internal = [[[]] for i in range(n_locations_internal + extra_bins)]
        for i in range(0, 4):
            locations_internal[0][0].append((bin_walls[(4 + i * 16)], bin_walls[(12 + i * 16)]))
            locations_internal[1][0].append((bin_walls[(2 + i * 16)], bin_walls[(6 + i * 16)]))
            locations_internal[1][0].append((bin_walls[(10 + i * 16)], bin_walls[(14 + i * 16)]))
            locations_internal[2][0].append((bin_walls[(1 + i * 16)], bin_walls[(3 + i * 16)]))
            locations_internal[2][0].append((bin_walls[(9 + i * 16)], bin_walls[(11 + i * 16)]))
            if kaufmann_design == "nested":
                locations_internal[3][0].append((bin_walls[(5 + i * 16)], bin_walls[(7 + i * 16)]))
                locations_internal[3][0].append((bin_walls[(13 + i * 16)], bin_walls[(15 + i * 16)]))
            else:
                locations_internal[2][0].append((bin_walls[(5 + i * 16)], bin_walls[(7 + i * 16)]))
                locations_internal[2][0].append((bin_walls[(13 + i * 16)], bin_walls[(15 + i * 16)]))
        if extra_bins > 0:  # TODO: fix this
            for j in range(extra_bins):
                for i in range(64 * (2 ** j)):
                    locations_internal[n_locations_internal + j][0].append((bin_walls_extra[int(
                        0 + i * ((2 ** extra_bins) / (2 ** j)))], bin_walls_extra[int(
                        ((2 ** extra_bins) / (2 ** j)) / 2 + i * ((2 ** extra_bins) / (2 ** j)))]))
        self.locations.extend(locations_internal)


class DiaWindows(object):
    """
    Create DIA window design
    """

    def __init__(self, ms1_mzs, ms1_range, dia_design, window_type, kaufmann_design, extra_bins, num_windows=None,
                 range_slack=0.01):
        ms1_range_difference = ms1_range[0][1] - ms1_range[0][0]
        # set the number of windows for kaufmann method
        if dia_design == "kaufmann":
            num_windows = 64
        # dont allow extra bins for basic method
        if dia_design == "basic" and extra_bins > 0:
            sys.exit("Cannot have extra bins with 'basic' dia design.")
        # find bin walls and extra bin walls
        if window_type == "even":
            internal_bin_walls = [ms1_range[0][0]]
            for window_index in range(0, num_windows):
                internal_bin_walls.append(ms1_range[0][0] + ((window_index + 1) / num_windows) * ms1_range_difference)
            internal_bin_walls[0] = internal_bin_walls[0] - range_slack * ms1_range_difference
            internal_bin_walls[-1] = internal_bin_walls[-1] + range_slack * ms1_range_difference
            internal_bin_walls_extra = None
            if extra_bins > 0:
                internal_bin_walls_extra = [ms1_range[0][0]]
                for window_index in range(0, num_windows * (2 ** extra_bins)):
                    internal_bin_walls_extra.append(ms1_range[0][0] + (
                            (window_index + 1) / (num_windows * (2 ** extra_bins))) * ms1_range_difference)
                internal_bin_walls_extra[0] = internal_bin_walls_extra[0] - range_slack * ms1_range_difference
                internal_bin_walls_extra[-1] = internal_bin_walls_extra[-1] + range_slack * ms1_range_difference
        elif window_type == "percentile":
            internal_bin_walls = np.percentile(ms1_mzs,
                                               np.arange(0, 100 + 100 / num_windows, 100 / num_windows)).tolist()
            internal_bin_walls[0] = internal_bin_walls[0] - range_slack * ms1_range_difference
            internal_bin_walls[-1] = internal_bin_walls[-1] + range_slack * ms1_range_difference
            internal_bin_walls_extra = None
            if extra_bins > 0:
                internal_bin_walls_extra = np.percentile(ms1_mzs,
                                                         np.arange(0, 100 + 100 / (num_windows * (2 ** extra_bins)),
                                                                   100 / (num_windows * (2 ** extra_bins)))).tolist()
                internal_bin_walls_extra[0] = internal_bin_walls_extra[0] - range_slack * ms1_range_difference
                internal_bin_walls_extra[-1] = internal_bin_walls_extra[-1] + range_slack * ms1_range_difference
        else:
            raise ValueError("Incorrect window_type selected. Must be 'even' or 'percentile'.")
            # convert bin walls and extra bin walls into locations to scan
        if dia_design == "basic":
            self.locations = []
            for window_index in range(0, num_windows):
                self.locations.append([[(internal_bin_walls[window_index], internal_bin_walls[window_index + 1])]])
        elif dia_design == "kaufmann":
            self.locations = KaufmannWindows(internal_bin_walls, internal_bin_walls_extra, kaufmann_design,
                                             extra_bins).locations
        else:
            raise ValueError("Incorrect dia_design selected. Must be 'basic' or 'kaufmann'.")


class DsDAController(Controller):
    def __init__(self, mass_spec, N, isolation_window, rt_tol, min_ms1_intensity, exclusion_list=None):
        super().__init__(mass_spec)
        self.last_ms1_scan = None
        self.N = N
        self.isolation_window = isolation_window  # the isolation window (in Dalton) around a precursor ion to be fragmented
        self.rt_tol = rt_tol  # the rt window to prevent the same precursor ion to be fragmented again
        self.min_ms1_intensity = min_ms1_intensity  # minimum ms1 intensity to fragment

        mass_spec.reset()

        # register new event handlers under this controller
        mass_spec.register(MassSpectrometer.MS_SCAN_ARRIVED, self.handle_scan)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_OPENING, self.handle_acquisition_open)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_CLOSING, self.handle_acquisition_closing)

    def run(self, schedule_file, progress_bar=True):
        self.schedule = pd.read_csv(schedule_file)
        for idx, row in self.schedule.iterrows():
            target_mass = row.targetMass
            target_time = row.targetTime

            if np.isnan(target_mass):
                ms_level = 1
                isolation_windows = [[(0, 1000)]]
                precursor = None
            else:
                ms_level = 2
                mz_lower = target_mass - self.isolation_window
                mz_upper = target_mass + self.isolation_window
                isolation_windows = [[(mz_lower, mz_upper)]]
                precursor_charge = +1 if (self.mass_spec.ionisation_mode == POSITIVE) else -1
                scan_id = 0
                precursor = Precursor(precursor_mz=target_mass, precursor_intensity=0,
                                      precursor_charge=precursor_charge, precursor_scan_id=scan_id)

            dda_scan_params = ScanParameters()
            dda_scan_params.set(ScanParameters.MS_LEVEL, ms_level)
            dda_scan_params.set(ScanParameters.ISOLATION_WINDOWS, isolation_windows)
            dda_scan_params.set(ScanParameters.TIME, target_time)
            if precursor:
                dda_scan_params.set(ScanParameters.PRECURSOR, precursor)
            self.mass_spec.add_to_queue(dda_scan_params)  # push this scan to the mass spec queue

        if progress_bar:
            with tqdm(total=target_time,
                      initial=0) as pbar:
                self.mass_spec.run(self.schedule, pbar=pbar)
        else:
            self.mass_spec.run(self.schedule)

    def handle_acquisition_open(self):
        self.logger.info('Acquisition open')

    def handle_acquisition_closing(self):
        self.logger.info('Acquisition closing')

    def _process_scan(self, scan):
        self.logger.info('Received {}'.format(scan))
        if scan.ms_level == 1:  # we get an ms1 scan, store it for fragmentation next time
            self.last_ms1_scan = scan
        elif scan.ms_level == 2:  # if we get ms2 scan, then do something with it
            # scan.filter_intensity(self.min_ms2_intensity)
            if scan.num_peaks > 0:
                self._plot_scan(scan)

    def _update_parameters(self, scan):
        pass




