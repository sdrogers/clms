import math

import numpy as np
import pylab as plt

from VMSfunctions.Chemicals import ChemicalCreator, UnknownChemical
from VMSfunctions.Common import *


class RoiToChemicalCreator(ChemicalCreator):
    def __init__(self, peak_sampler, data_source, filename=None, min_ms1_intensity=None):
        super().__init__(peak_sampler)

        # if filename is specified, then we use the ROIs only from that file
        # otherwise we use the ROIs from all files found in the data source
        if filename is not None:
            rois_data = data_source.all_rois[filename]['rois']
        else: # combine the extracted ROIs for all files
            rois_data = []
            for filename in data_source.all_rois:
                rois_data.extend(data_source.all_rois[filename]['rois'])

        self.ms_levels = 2
        self.crp_samples = [[] for i in range(self.ms_levels)]
        self.crp_index = [[] for i in range(self.ms_levels)]
        self.alpha = math.inf
        self.counts = [[] for i in range(self.ms_levels)]
        if self.ms_levels > 2:
            self.logger.warning("Warning ms_level > 3 not implemented properly yet. Uses scaled ms_level = 2 information for now")

        # collect the regions of interest that contain no peaks
        false_rois = [roi for roi in rois_data if not roi.pickedPeak]
        self.chromatograms = []
        self.chemicals = []
        for i in range(len(false_rois)):
            if i % 50000 == 0:
                self.logger.debug('%6d/%6d' % (i, len(false_rois)))
            # if yes, then try to turn this into a chromatogram and unknown chemical
            roi = false_rois[i]

            # raise numpy warning as exception, see https://stackoverflow.com/questions/15933741/how-do-i-catch-a-numpy-warning-like-its-an-exception-not-just-for-testing
            chrom = None
            with np.errstate(divide='raise'):
                try:
                    chrom = roi.to_chromatogram()
                except FloatingPointError:
                    self.logger.debug('Invalid chromatogram {}'.format(i))
                except ZeroDivisionError:
                    self.logger.debug('Invalid chromatogram {}'.format(i))

            if chrom is not None and self._valid_roi(roi, min_ms1_intensity=min_ms1_intensity):
                chem = self.to_unknown_chemical(chrom)
                try:
                    chem.children = self._get_children(1, chem) # TODO: this should happen inside the mass spec class
                except KeyError:
                    pass
                self.chromatograms.append(chrom)
                self.chemicals.append(chem)
        assert len(self.chromatograms) == len(self.chemicals)
        self.logger.info('Found %d ROIs above thresholds' % len(self.chromatograms))

    def sample(self, N, min_num_scans=None):
        chemicals = []
        while len(chemicals) < N:
            # pick one chemical with roi and check that it's valid
            chem = np.random.choice(len(self.chemicals), 1)[0]
            roi = chem.chromatogram.roi
            if roi is not None and self._valid_roi(roi, min_num_scans=min_num_scans):
                chemicals.append(chem)
        return chemicals

    def _valid_roi(self, roi, min_num_scans=None, min_ms1_intensity=None):
        if min_num_scans is not None and roi.num_scans() < min_num_scans:
            return False
        if min_ms1_intensity is not None and max(roi.intensities()) < min_ms1_intensity:
            return False
        return True

    def to_unknown_chemical(self, chrom):
        idx = np.argmax(chrom.raw_intensities) # find intensity apex
        mz = chrom.raw_mzs[idx]

        # In the MassSpec, we assume that chemical starts eluting from chem.rt + chem.chromatogram.rts (normalised to start from 0)
        # So here, we have to set set chemical rt to start from the minimum of chromatogram raw rts, so it elutes correct.
        # rt = chrom.raw_rts[idx]
        rt = min(chrom.raw_rts)

        max_intensity = chrom.raw_intensities[idx]
        chem = UnknownChemical(mz, rt, max_intensity, chrom, None)
        chem.type = CHEM_NOISE
        return chem

    def plot_chems(self, n_plots, reverse=False):
        sorted_chems = sorted(self.chemicals, key=lambda chem: chem.chromatogram.roi.num_scans())
        if reverse:
            sorted_chems.reverse()
        for c in sorted_chems[0:n_plots]:
            chrom = c.chromatogram
            plt.plot(chrom.raw_rts, chrom.raw_intensities)
            plt.show()