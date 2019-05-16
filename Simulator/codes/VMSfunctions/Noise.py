import math
import numpy as np
import pylab as plt

from VMSfunctions.Chemicals import ChemicalCreator, UnknownChemical
from VMSfunctions.Common import *


class RoiToChemicalCreator(ChemicalCreator):
    def __init__(self, peak_sampler, all_roi):
        super().__init__(peak_sampler)
        self.rois_data = all_roi
        self.ms_levels = 2
        self.crp_samples = [[] for i in range(self.ms_levels)]
        self.crp_index = [[] for i in range(self.ms_levels)]
        self.alpha = math.inf
        self.counts = [[] for i in range(self.ms_levels)]
        if self.ms_levels > 2:
            self.logger.warning(
                "Warning ms_level > 3 not implemented properly yet. Uses scaled ms_level = 2 information for now")

        self.chromatograms = []
        self.chemicals = []
        for i in range(len(self.rois_data)):
            if i % 50000 == 0:
                self.logger.debug('%6d/%6d' % (i, len(self.rois_data)))
            roi = self.rois_data[i]

            # raise numpy warning as exception, see https://stackoverflow.com/questions/15933741/how-do-i-catch-a-numpy-warning-like-its-an-exception-not-just-for-testing
            chrom = None
            with np.errstate(divide='raise'):
                try:
                    chrom = roi.to_chromatogram()
                except FloatingPointError:
                    self.logger.debug('Invalid chromatogram {}'.format(i))
                except ZeroDivisionError:
                    self.logger.debug('Invalid chromatogram {}'.format(i))

            if chrom is not None:
                chem = self.to_unknown_chemical(chrom)
                try:
                    chem.children = self._get_children(1, chem, n_peaks=1)  # TODO: this should happen inside the mass spec class
                except KeyError:
                    pass
                self.chromatograms.append(chrom)
                self.chemicals.append(chem)
        assert len(self.chromatograms) == len(self.chemicals)
        self.logger.info('Found %d ROIs above thresholds' % len(self.chromatograms))


    def sample(self, chromatogram_creator, mz_range, rt_range, min_ms1_intensity, n_ms1_peaks, ms_levels=2,
               chemical_type=None,
               formula_list=None, compound_list=None, alpha=math.inf, fixed_mz=False, adduct_proportion_cutoff=0.05):
        return NotImplementedError()


    def sample_from_chromatograms(self, chromatogram_creator, min_rt, max_rt, min_ms1_intensity, ms_levels=2):
        return NotImplementedError()


    def to_unknown_chemical(self, chrom):
        idx = np.argmax(chrom.raw_intensities)  # find intensity apex
        mz = chrom.raw_mzs[idx]

        # In the MassSpec, we assume that chemical starts eluting from chem.rt + chem.chromatogram.rts (normalised to start from 0)
        # So here, we have to set set chemical rt to start from the minimum of chromatogram raw rts, so it elutes correct.
        # rt = chrom.raw_rts[idx]
        rt = min(chrom.raw_rts)

        max_intensity = chrom.raw_intensities[idx]
        mz = mz - PROTON_MASS
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
