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


# class Controller:

# class DIAController(Controller):

# class TopNController(Controller):