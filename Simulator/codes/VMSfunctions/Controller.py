from collections import defaultdict, namedtuple

import pylab as plt

from VMSfunctions.MassSpec import *
from VMSfunctions.Common import *
from psims.mzml.writer import MzMLWriter


class Controller(object):
    def __init__(self, mass_spec):
        self.scans = defaultdict(list) # key: ms level, value: list of scans for that level
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
    def __init__(self, mass_spec):
        super().__init__(mass_spec)
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


ExclusionItem = namedtuple('ExclusionItem', 'from_mz to_mz from_rt to_rt')

Precursor = namedtuple('Precursor', 'precursor_mz precursor_intensity precursor_charge precursor_scan_id')

class TopNController(Controller):
    def __init__(self, mass_spec, N, mz_tol, rt_tol, exclusion_list=[], min_ms2_intensity=0):
        super().__init__(mass_spec)
        self.last_ms1_scan = None
        self.N = N
        self.mz_tol = mz_tol # the m/z window around a precursor ion to be fragmented
        self.rt_tol = rt_tol # the rt window to prevent the same precursor ion to be fragmented again
        self.exclusion_list = exclusion_list # a list of ExclusionItem
        self.min_ms2_intensity = min_ms2_intensity

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

        if scan.ms_level == 1: # if we get a non-empty ms1 scan
            if scan.num_peaks > 0:
                print(scan)
                self.last_ms1_scan = scan
            else:
                self.last_ms1_scan = None

        elif scan.ms_level == 2: # if we get ms2 scan, then do something with it
            scan.filter_intensity(self.min_ms2_intensity)
            if scan.num_peaks > 0:
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

            rt = self.last_ms1_scan.rt

            # then get the last ms1 scan and select its top-N precursor ions
            intensities = self.last_ms1_scan.intensities
            largest_indices = intensities.argsort()[-self.N:][::-1]
            largest_mzs = self.last_ms1_scan.mzs[largest_indices]
            largest_intensities = self.last_ms1_scan.intensities[largest_indices]

            for i in range(len(largest_mzs)): # define isolation window around the selected precursor ions

                mz = largest_mzs[i]
                intensity = largest_intensities[i]
                if self._exclude(mz, rt, self.exclusion_list):
                    continue

                if self.mass_spec.ionisation_mode == POSITIVE:
                    precursor_charge = +1 # assume it's all +1 if positive
                elif self.mass_spec.ionisation_mode == NEGATIVE:
                    precursor_charge = -1

                precursor = Precursor(precursor_mz=mz, precursor_intensity=intensity,
                                      precursor_charge=precursor_charge, precursor_scan_id=self.last_ms1_scan.scan_id)

                mz_lower = mz * (1 - self.mz_tol / 1e6)
                mz_upper = mz * (1 + self.mz_tol / 1e6)
                isolation_windows = [[(mz_lower, mz_upper)]]
                print('Isolated precursor ion', mz, 'window', isolation_windows)
                dda_scan_params = ScanParameters()
                dda_scan_params.set(ScanParameters.MS_LEVEL, 2)
                dda_scan_params.set(ScanParameters.ISOLATION_WINDOWS, isolation_windows)
                dda_scan_params.set(ScanParameters.PRECURSOR, precursor)
                self.mass_spec.add_to_queue(dda_scan_params) # push this dda scan to the mass spec queue

                # dynamic exclusion: prevent the same precursor ion being fragmented multiple times in the same windows
                rt_lower = rt - self.rt_tol
                rt_upper = rt + self.rt_tol
                x = ExclusionItem(from_mz=mz_lower, to_mz=mz_upper, from_rt=rt_lower, to_rt=rt_upper)
                self.exclusion_list.append(x)

            print()

            # set this ms1 scan as has been processed
            self.last_ms1_scan = None

            # remove expired items from exclusion list
            self.exclusion_list = list(filter(lambda x: x.to_rt > rt, self.exclusion_list))


    def write_mzML(self, analysis_name, out_file):
        ms1_id_to_scan = {x.scan_id: x for x in self.scans[1]}
        spectrum_count = self._get_spectrum_count(ms1_id_to_scan)

        ms1_id_to_precursors = defaultdict(list)
        for p in self.mass_spec.precursor_information:
            ms1_id_to_precursors[p.precursor_scan_id].append(p)

        with MzMLWriter(open(out_file, 'wb')) as out:

            # add default controlled vocabularies
            out.controlled_vocabularies()

            # TODO: correctly set these fields
            out.file_description(
                file_contents={
                    'MSn spectrum': '',
                    'centroid spectrum': ''
                },
                source_files=[]
            )
            out.sample_list(samples=[])
            out.software_list(software_list={
                'id': 'VMS',
                'version': '1.0.0'
            })
            out.scan_settings_list(scan_settings=[])
            out.instrument_configuration_list(instrument_configurations={
                'id': 'VMS',
                'component_list': []
            })
            out.data_processing_list({'id': 'VMS'})

            # open the run
            with out.run(id=analysis_name):

                # open spectrum list sections
                with out.spectrum_list(count=spectrum_count):
                    for ms1_id in sorted(ms1_id_to_precursors.keys()):
                        # write ms1 scan
                        ms1_scan = ms1_id_to_scan[ms1_id]
                        out.write_spectrum(
                            ms1_scan.mzs, ms1_scan.intensities,
                            id=ms1_scan.scan_id, params=[
                                'MS1 Spectrum',
                                {'ms level': 1},
                                {'total ion current': np.sum(ms1_scan.intensities)}
                            ])

                        # get all precursor ions in this ms1 scan
                        precursors = ms1_id_to_precursors[ms1_id]
                        for precursor in precursors:
                            # get all ms2 scans produced from this precursor ion
                            ms2_scans = self.mass_spec.precursor_information[precursor]
                            for prod in ms2_scans: # write ms2 scan information
                                out.write_spectrum(
                                    prod.mzs, prod.intensities,
                                    id=prod.scan_id, params=[
                                        'MSn Spectrum',
                                        {'ms level': 2},
                                        {'total ion current': np.sum(prod.intensities)}
                                    ],
                                    # also include precursor information
                                    precursor_information={
                                        "mz": precursor.precursor_mz,
                                        "intensity": precursor.precursor_intensity,
                                        "charge": precursor.precursor_charge,
                                        "scan_id": precursor.precursor_scan_id
                                    })

                # open chromatogram list sections
                with out.chromatogram_list(count=1):
                        tic_rts, tic_intensities = self._get_tic_chromatogram()
                        out.write_chromatogram(tic_rts, tic_intensities, id='tic',
                                               chromatogram_type='total ion current chromatogram',
                                               time_unit='second')


        out.close()

    def _exclude(self, mz, rt, exclusion_list): # TODO: make this faster?
        for x in exclusion_list:
            exclude_mz = x.from_mz < mz < x.to_mz
            exclude_rt = x.from_rt < rt < x.to_rt
            if exclude_mz and exclude_rt:
                print('Excluded precursor ion', mz, rt, x)
                return True
        return False

    def _get_spectrum_count(self, ms1_id_to_scan):
        spectrum_count = 0
        for precursor, scans in self.mass_spec.precursor_information.items():
            ms1_scan = ms1_id_to_scan[precursor.precursor_scan_id]
            ms2_scans = scans
            spectrum_count += ms1_scan.num_peaks
            spectrum_count += np.sum([x.num_peaks for x in scans])
        return spectrum_count

    def _get_tic_chromatogram(self):
        time_array = []
        intensity_array = []
        for ms1_scan in self.scans[1]:
            time_array.append(ms1_scan.rt)
            intensity_array.append(np.sum(ms1_scan.intensities))
        time_array = np.array(time_array)
        intensity_array = np.array(intensity_array)
        return time_array, intensity_array

# class Controller:

# class DIAController(Controller):

# class TopNController(Controller):