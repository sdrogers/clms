import math

import numpy as np

from VMSfunctions.Chemicals import ChemicalCreator, UnknownChemical
from VMSfunctions.Common import *

logger = get_logger('Noise')


class NoisyRoiCreator(ChemicalCreator):
    def __init__(self, peak_sampler, data_source, filename=None):
        super().__init__(peak_sampler)

        # if filename is specified, then we use the ROIs only from that file
        # otherwise we use the ROIs from all files found in the data source
        if filename is not None:
            rois_data = data_source.all_rois[filename]
        else: # combine the extracted ROIs for all files
            rois_data = []
            for filename in data_source.all_rois:
                rois_data.extend(data_source.all_rois[filename])

        # collect the regions of interest that contain no peaks
        self.false_rois = [roi for roi in rois_data['rois'] if not roi.pickedPeak]

    def sample(self, N, ms_levels=2, min_num_scans=None):
        self.ms_levels = ms_levels
        self.crp_samples = [[] for i in range(self.ms_levels)]
        self.crp_index = [[] for i in range(self.ms_levels)]
        self.alpha = math.inf
        self.counts = [[] for i in range(self.ms_levels)]
        if self.ms_levels > 2:
            logger.warning("Warning ms_level > 3 not implemented properly yet. Uses scaled ms_level = 2 information for now")

        chemicals = []
        while len(chemicals) < N:
            # pick one roi and check that it's valid
            selected = np.random.choice(len(self.false_rois), 1)[0]
            roi = self.false_rois[selected]
            if self._valid_roi(roi, min_num_scans):
                # if yes, then try to turn this into a chromatogram and unknown chemical
                chrom = roi.to_chromatogram()
                if chrom is not None:
                    chem = self.to_unknown_chemical(chrom)
                    chem.children = self._get_children(1, chem)
                    chemicals.append(chem)
        return chemicals

    def _valid_roi(self, roi, min_scans):
        if min_scans is not None:
            if roi.num_scans() >= min_scans:
                return True
            else:
                return False
        return True

    def to_unknown_chemical(self, chrom):
        idx = np.argmax(chrom.raw_intensities)
        mz = chrom.raw_mzs[idx]
        rt = chrom.raw_rts[idx]
        max_intensity = chrom.raw_intensities[idx]
        chem = UnknownChemical(mz, rt, max_intensity, chrom, None)
        return chem
