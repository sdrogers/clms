from collections import defaultdict

import numpy as np
from psims.mzml.writer import MzMLWriter as PsimsMzMLWriter

from VMSfunctions.Common import DEFAULT_MS1_SCAN_WINDOW


class MzmlWriter(object):
    """A class to write peak data to mzML file"""

    def __init__(self, analysis_name, scans, precursor_information=None):
        """
        Initialises the mzML writer class.
        :param analysis_name: Name of the analysis.
        :param scans: A dictionary where key is scan level, value is a list of Scans object for that level.
        :param precursor_information: A dictionary where key is Precursor object, value is a list of ms2 scans only
        """
        self.analysis_name = analysis_name
        self.scans = scans
        self.precursor_information = precursor_information


    def write_mzML(self, out_file):
        with PsimsMzMLWriter(open(out_file, 'wb')) as writer:

            # add default controlled vocabularies
            writer.controlled_vocabularies()

            # write other fields like sample list, software list, etc.
            self._write_info(writer)

            # open the run
            with writer.run(id=self.analysis_name):
                self._write_spectra(writer, self.scans, self.precursor_information)

                # open chromatogram list sections
                with writer.chromatogram_list(count=1):
                    tic_rts, tic_intensities = self._get_tic_chromatogram(self.scans)
                    writer.write_chromatogram(tic_rts, tic_intensities, id='tic',
                                           chromatogram_type='total ion current chromatogram',
                                           time_unit='second')

        writer.close()


    def _write_info(self, out):
        # check file contains what kind of spectra
        has_ms1_spectrum = 1 in self.scans
        has_msn_spectrum = 1 in self.scans and len(self.scans) > 1
        file_contents = [
            'centroid spectrum'
        ]
        if has_ms1_spectrum:
            file_contents.append('MS1 spectrum')
        if has_msn_spectrum:
            file_contents.append('MSn spectrum')
        out.file_description(
            file_contents=file_contents,
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

    def sort_filter(self, all_scans):
        all_scans = sorted(all_scans, key=lambda x: x.scan_id)
        all_scans = [x for x in all_scans if x.num_peaks > 0]
        return all_scans

    def _write_spectra(self, writer, scans, precursor_information):
        assert len(scans) <= 3 # NOTE: we only support writing up to ms2 scans for now

        # get all scans across different ms_levels and sort them by scan_id
        all_scans = []
        for ms_level in scans:
            all_scans.extend(scans[ms_level])
        all_scans = self.sort_filter(all_scans)
        spectrum_count = len(all_scans)

        # get precursor information for each scan, if available
        scan_precursor = {}
        for precursor, ms2_scans in precursor_information.items():
            assert len(ms2_scans) == 1
            ms2_scan = ms2_scans[0]
            scan_precursor[ms2_scan.scan_id] = precursor

        # write scans
        with writer.spectrum_list(count=spectrum_count):
            for scan in all_scans:
                precursor = None
                if scan.scan_id in scan_precursor:
                    precursor = scan_precursor[scan.scan_id]
                self._write_scan(writer, scan, precursor)

    def _get_scan_id(self, scan_id):
        return scan_id

    def _write_scan(self, out, scan, precursor):
        assert scan.num_peaks > 0
        label = 'MS1 Spectrum' if scan.ms_level == 1 else 'MSn Spectrum'
        precursor_information = None
        if precursor is not None:
            precursor_information = {
                "mz": precursor.precursor_mz,
                "intensity": precursor.precursor_intensity,
                "charge": precursor.precursor_charge,
                "spectrum_reference": self._get_scan_id(precursor.precursor_scan_id),
                "activation": ["HCD", {"collision energy": 25.0}]
            }
        lowest_observed_mz = min(scan.mzs)
        highest_observed_mz = max(scan.mzs)
        bp_pos = np.argmax(scan.intensities)
        bp_intensity = scan.intensities[bp_pos]
        bp_mz = scan.mzs[bp_pos]

        out.write_spectrum(
            scan.mzs, scan.intensities,
            id= self._get_scan_id(scan.scan_id),
            centroided=True,
            scan_start_time=scan.rt / 60.0,
            scan_window_list=[DEFAULT_MS1_SCAN_WINDOW],
            params=[
                {label: ''},
                {'ms level': scan.ms_level},
                {'total ion current': np.sum(scan.intensities)},
                {'lowest observed m/z': lowest_observed_mz},
                {'highest observed m/z': highest_observed_mz},
                # {'base peak m/z', bp_mz},
                # {'base peak intensity', bp_intensity}
            ],
            precursor_information=precursor_information
        )

    def _get_tic_chromatogram(self, scans):
        time_array = []
        intensity_array = []
        for ms1_scan in scans[1]:
            time_array.append(ms1_scan.rt)
            intensity_array.append(np.sum(ms1_scan.intensities))
        time_array = np.array(time_array)
        intensity_array = np.array(intensity_array)
        return time_array, intensity_array
