from collections import defaultdict

import numpy as np
from psims.mzml.writer import MzMLWriter as PsimsMzMLWriter

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
        with PsimsMzMLWriter(open(out_file, 'wb')) as out:

            # add default controlled vocabularies
            out.controlled_vocabularies()

            # write other fields like sample list, software list, etc.
            self._write_info(out)

            # open the run
            with out.run(id=self.analysis_name):
                if len(self.precursor_information) > 0: # write spectra based on precursor information
                    self._write_spectra_with_precursors(out, self.precursor_information, self.scans)
                else: # write spectra without precursor information
                    self._write_spectra(out, self.scans)

                # open chromatogram list sections
                with out.chromatogram_list(count=1):
                    tic_rts, tic_intensities = self._get_tic_chromatogram(self.scans)
                    out.write_chromatogram(tic_rts, tic_intensities, id='tic',
                                           chromatogram_type='total ion current chromatogram',
                                           time_unit='second')

        out.close()


    def _write_info(self, out):
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


    def _write_spectra_with_precursors(self, out, precursor_information, scans):
        ms1_id_to_scan = {x.scan_id: x for x in scans[1]}
        ms1_id_to_precursors = defaultdict(list)
        for p in precursor_information:
            ms1_id_to_precursors[p.precursor_scan_id].append(p)
        spectrum_count = len(ms1_id_to_precursors) + sum([len(products) for _, products in precursor_information.items()])

        # open spectrum list sections
        with out.spectrum_list(count=spectrum_count):
            for ms1_id in sorted(ms1_id_to_precursors.keys()):
                # write ms1 scan
                ms1_scan = ms1_id_to_scan[ms1_id]
                self._write_scan(out, ms1_scan)

                # get all precursor ions in this ms1 scan
                precursors = ms1_id_to_precursors[ms1_id]
                for precursor in precursors:
                    # get all ms2 scans produced from this precursor ion
                    ms2_scans = precursor_information[precursor]
                    for prod in ms2_scans:  # write ms2 scan information
                        self._write_scan(out, prod, precursor=precursor)

    def _write_spectra(self, out, scans):
        # get all scans across different ms_levels and sort them by scan_id
        all_scans = []
        for ms_level in scans:
            all_scans.extend(scans[ms_level])
        all_scans = sorted(all_scans, key=lambda x: x.scan_id)
        spectrum_count = len(all_scans)

        # write scans
        with out.spectrum_list(count=spectrum_count):
            for scan in all_scans:
                self._write_scan(out, scan)

    def _write_scan(self, out, scan, precursor=None):
        if scan.num_peaks > 0:
            label = 'MS1 Spectrum' if scan.ms_level == 1 else 'MSn Spectrum'
            precursor_information = None if precursor is None else {
                "mz": precursor.precursor_mz,
                "intensity": precursor.precursor_intensity,
                "charge": precursor.precursor_charge,
                "scan_id": precursor.precursor_scan_id
            }
            out.write_spectrum(
                scan.mzs, scan.intensities,
                id=scan.scan_id,
                scan_start_time=scan.rt / 60.0,
                params=[
                    label,
                    {'ms level': scan.ms_level},
                    {'total ion current': np.sum(scan.intensities)}
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
