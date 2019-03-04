from collections import defaultdict

import pylab as plt

from VMSfunctions.MassSpec import *


class Controller(object):
    def __init__(self, chemicals, mass_spec):
        self.scans = defaultdict(list) # key: ms level, value: list of scans for that level
        self.chemicals = chemicals
        self.mass_spec = mass_spec

    def handle_scan(self, scan):
        self.scans[scan.ms_level].append(scan)
        self.update_parameters(scan)

    def handle_acquisition_open(self):
        raise NotImplementedError()

    def handle_acquisition_closing(self):
        raise NotImplementedError()

    def update_parameters(self, scan):
        raise NotImplementedError()

    def _plot_scan(self, scan):
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
    def __init__(self, chemicals, mass_spec):
        super().__init__(chemicals, mass_spec)
        default_scan = ScanParameters()
        default_scan.set(ScanParameters.MS_LEVEL, 1)
        default_scan.set(ScanParameters.ISOLATION_WINDOWS, [[(0, 1e3)]])
        mass_spec.set_repeating_scan(default_scan)
        mass_spec.register(MassSpectrometer.MS_SCAN_ARRIVED, self.handle_scan)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_OPENING, self.handle_acquisition_open)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_CLOSING, self.handle_acquisition_closing)

    def run(self, max_time):
        self.mass_spec.run(max_time)

    def handle_scan(self, scan):
        super().handle_scan(scan)
        if scan.num_peaks > 0:
            self._plot_scan(scan)
            for mz, intensity in zip(scan.mzs, scan.intensities):
                print(mz, intensity)

    def handle_acquisition_open(self):
        print('Acquisition open')

    def handle_acquisition_closing(self):
        print('Acquisition closing')

    def update_parameters(self, scan):
        pass # do nothing


class TopNController(Controller):
    def __init__(self, chemicals, mass_spec, N, mz_tol):
        super().__init__(chemicals, mass_spec)
        self.last_ms1_scan = None
        self.N = N
        self.mz_tol = mz_tol

        default_scan = ScanParameters()
        default_scan.set(ScanParameters.MS_LEVEL, 1)
        default_scan.set(ScanParameters.ISOLATION_WINDOWS, [[(0, 1e3)]])
        mass_spec.set_repeating_scan(default_scan)

        mass_spec.register(MassSpectrometer.MS_SCAN_ARRIVED, self.handle_scan)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_OPENING, self.handle_acquisition_open)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_CLOSING, self.handle_acquisition_closing)

    def run(self, max_time):
        self.mass_spec.run(max_time)

    def handle_scan(self, scan):
        super().handle_scan(scan)
        if scan.num_peaks > 0:
            print(scan)

            if scan.ms_level == 1: # if we get a non-empty ms1 scan
                if scan.num_peaks > 0:
                    self.last_ms1_scan = scan
                else:
                    self.last_ms1_scan = None

            elif scan.ms_level == 2: # if we get ms2 scan, then just show it
                self._plot_scan(scan)
                for mz, intensity in zip(scan.mzs, scan.intensities):
                    print(mz, intensity)
                print()

        self.update_parameters(scan)

    def handle_acquisition_open(self):
        print('Acquisition open')

    def handle_acquisition_closing(self):
        print('Acquisition closing')

    def update_parameters(self, scan):

        # if there's a previous ms1 scan to process
        if self.last_ms1_scan is not None:

            # then get the last ms1 scan and select its top-N precursor ions
            intensities = self.last_ms1_scan.intensities
            largest_indices = intensities.argsort()[-self.N:][::-1]
            largest_mzs = self.last_ms1_scan.mzs[largest_indices]

            for i in range(len(largest_mzs)): # define isolation window around the selected precursor ions
                mz = largest_mzs[i]
                mz_lower = mz * (1 - self.mz_tol / 1e6)
                mz_upper = mz * (1 + self.mz_tol / 1e6)
                isolation_windows = [[(mz_lower, mz_upper)]]
                print('Isolated precursor ion', mz, 'window', isolation_windows)
                dda_scan_params = ScanParameters()
                dda_scan_params.set(ScanParameters.MS_LEVEL, 2)
                dda_scan_params.set(ScanParameters.ISOLATION_WINDOWS, isolation_windows)
                self.mass_spec.add_to_queue(dda_scan_params) # push this dda scan to the mass spec queue
            print()

            # set this ms1 scan as has been processed
            self.last_ms1_scan = None

# class Controller:

# class DIAController(Controller):

# class TopNController(Controller):