import pylab as plt

from VMSfunctions.MassSpec import *


class Controller(object):
    def __init__(self, chemicals, mass_spec):
        self.scans = []
        self.chemicals = chemicals
        self.mass_spec = mass_spec

    def handle_scan(self, scan):
        raise NotImplementedError()

    def handle_acquisition_open(self):
        raise NotImplementedError()

    def handle_acquisition_closing(self):
        raise NotImplementedError()

    def update_parameters(self):
        raise NotImplementedError()


class SimpleMs1Controller(Controller):
    def __init__(self, chemicals, mass_spec, plot_peaks=25):
        super().__init__(chemicals, mass_spec)
        self.scan_parameters = {
            'isolation_windows': [[(0, 1e3)]],  # TODO: change to dictionary?
            'ms_level': 1
        }
        self.plot_peaks = plot_peaks
        mass_spec.register(MassSpectrometer.MS_SCAN_ARRIVED, self.handle_scan)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_OPENING, self.handle_acquisition_open)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_CLOSING, self.handle_acquisition_closing)

    def run(self, max_time):
        self.mass_spec.run(max_time)

    def handle_scan(self, scan):
        self.scans.append(scan)
        if scan.num_peaks > self.plot_peaks:
            self._plot_scan(scan)
            for mz, intensity in zip(scan.mzs, scan.intensities):
                print(mz, intensity)
        self.update_parameters()

    def handle_acquisition_open(self):
        print('Acquisition open')

    def handle_acquisition_closing(self):
        print('Acquisition closing')

    def update_parameters(self):
        pass # do nothing

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


class TopNController(Controller):
    def __init__(self, chemicals, mass_spec, plot_peaks=25):
        super().__init__(chemicals, mass_spec)
        self.scan_parameters = {
            'isolation_windows': [[(0, 1e3)]],  # TODO: change to dictionary?
            'ms_level': 1
        }
        self.plot_peaks = plot_peaks
        mass_spec.register(MassSpectrometer.MS_SCAN_ARRIVED, self.handle_scan)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_OPENING, self.handle_acquisition_open)
        mass_spec.register(MassSpectrometer.ACQUISITION_STREAM_CLOSING, self.handle_acquisition_closing)

    def run(self, max_time):
        self.mass_spec.run(max_time)

    def handle_scan(self, scan):
        self.scans.append(scan)
        if scan.num_peaks > 0:
            print(scan)
            if scan.ms_level == 2:
                self._plot_scan(scan)
                for mz, intensity in zip(scan.mzs, scan.intensities):
                    print(mz, intensity)
        self.update_parameters(scan)

    def handle_acquisition_open(self):
        print('Acquisition open')

    def handle_acquisition_closing(self):
        print('Acquisition closing')

    def update_parameters(self, scan):
        if len(self.scans) % 5 == 0:
            N = 5
            last_scan = self.scans[-1]
            if scan.num_peaks > 0:
                # print('topN')
                intensities = last_scan.intensities
                largest_indices = intensities.argsort()[-N:][::-1]
                largest_mzs = last_scan.mzs[largest_indices]
                largst_intensities = intensities[largest_indices]
                # print(largest_mzs)
                # print(largst_intensities)
                for i in range(len(largest_mzs)):
                    mz = largest_mzs[i]
                    isolation_windows = [[(mz-1, mz+1)]]
                    params = {
                        'isolation_windows': isolation_windows,  # TODO: change to dictionary?
                        'ms_level': 2
                    }
                    # print(params)
                    self.mass_spec.add_to_queue(params)
        else:
            # print('MS1 scan')
            # TODO: make the mass spec runs a default scan parameter without being told by the controller
            self.scan_parameters = {
                'isolation_windows': [[(0, 1e3)]],  # TODO: change to dictionary?
                'ms_level': 1
            }

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


# class Controller:

# class DIAController(Controller):

# class TopNController(Controller):