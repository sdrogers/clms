import numpy as np
import pandas as pd
import scipy
import scipy.stats

from VMSfunctions.Common import chromatogramDensityNormalisation


class Chromatogram(object):

    def get_relative_intensity(self, query_rt):
        raise NotImplementedError()

    def get_relative_mz(self, query_rt):
        raise NotImplementedError()

    def _rt_match(self, rt):
        raise NotImplementedError()


class EmpiricalChromatogram(Chromatogram):
    """
    Empirical Chromatograms to be used within Chemicals
    """

    def __init__(self, rts, mzs, intensities):
        self.rts = [x - min(rts) for x in rts]
        self.mzs = [x - sum(mzs) / len(mzs) for x in mzs]  # may want to just set this to 0 and remove from input
        self.intensities = chromatogramDensityNormalisation(rts, intensities)
        self.raw_rts = rts
        self.raw_mzs = mzs
        self.raw_intensities = intensities

    def get_relative_intensity(self, query_rt):
        if self._rt_match(query_rt) == False:
            return None
        else:
            return ((self.intensities[self._get_rt_neighbours_which(query_rt)[0]] +
                     (self.intensities[self._get_rt_neighbours_which(query_rt)[1]]
                      - self.intensities[self._get_rt_neighbours_which(query_rt)[0]]) * self._get_distance(query_rt)))

    def get_relative_mz(self, query_rt):
        if self._rt_match(query_rt) == False:
            return None
        else:
            return ((self.mzs[self._get_rt_neighbours_which(query_rt)[0]] +
                     (self.mzs[self._get_rt_neighbours_which(query_rt)[1]]
                      - self.mzs[self._get_rt_neighbours_which(query_rt)[0]]) * self._get_distance(query_rt)))

    def _get_rt_neighbours(self, query_rt):
        rt_below = max(x for x in self.rts if x <= query_rt)
        rt_above = min(x for x in self.rts if x >= query_rt)
        return ([rt_below, rt_above])

    def _get_rt_neighbours_which(self, query_rt):
        which_rt_below = self.rts.index(self._get_rt_neighbours(query_rt)[0])
        which_rt_above = self.rts.index(self._get_rt_neighbours(query_rt)[1])
        return ([which_rt_below, which_rt_above])

    def _get_distance(self, query_rt):
        return ((query_rt - self._get_rt_neighbours(query_rt)[0]) /
                (self._get_rt_neighbours(query_rt)[1] - self._get_rt_neighbours(query_rt)[0]))

    def _rt_match(self, query_rt):
        if query_rt < min(self.rts) or query_rt > max(self.rts):
            return False
        else:
            return True


# Make this more generalisable. Make scipy.stats... as input, However this makes it difficult to do the cutoff
class FunctionalChromatogram(Chromatogram):
    """
    Functional Chromatograms to be used within Chemicals
    """

    def __init__(self, distribution, parameters, cutoff=0.01):
        self.cutoff = cutoff
        self.mz = 0
        if distribution == "normal":
            self.distrib = scipy.stats.norm(parameters[0], parameters[1])
        elif distribution == "gamma":
            self.distrib = scipy.stats.gamma(parameters[0], parameters[1], parameters[2])
        elif distribution == "uniform":
            self.distrib = scipy.stats.uniform(parameters[0], parameters[1])
        else:
            raise NotImplementedError("distribution not implemented")

    def get_relative_intensity(self, query_rt):
        if self._rt_match(query_rt) == False:
            return None
        else:
            return (self.distrib.pdf(query_rt + self.distrib.ppf(self.cutoff / 2)) * (1 / (1 - self.cutoff)))

    def get_relative_mz(self, query_rt):
        if self._rt_match(query_rt) == False:
            return None
        else:
            return self.mz

    def _rt_match(self, query_rt):
        if query_rt < 0 or query_rt > self.distrib.ppf(1 - (self.cutoff / 2)) - self.distrib.ppf(self.cutoff / 2):
            return False
        else:
            return True


class ChromatogramCreator(object):
    def __init__(self, xcms_output=None):
        self.xcms_output = xcms_output
        if self.xcms_output != None:
            self.chromatograms = self._load_chromatograms(self.xcms_output)
        else:
            self.chromatograms = None

    def sample(self):
        if self.chromatograms != None:
            selected = np.random.choice(len(self.chromatograms), 1)[0]
            return self.chromatograms[selected]
        else:
            NotImplementedError("Functional Chromatograms not implemented here yet")

    def _load_chromatograms(self, xcms_output):
        return self._load_xcms_df(xcms_output)

    def _load_xcms_df(self, df_file):
        """
        Load CSV file of chromatogram information exported by the XCMS script 'process_data.R'
        :param df_file: the input csv file exported by the script (in gzip format)
        :return: a list of Chromatogram objects
        """
        df = pd.read_csv(df_file, compression='gzip')
        peak_ids = df.id.unique()
        groups = df.groupby('id')
        chroms = []
        for i in range(len(peak_ids)):
            if i % 5000 == 0:
                print(i)
            pid = peak_ids[i]
            chrom = self._get_xcms_chromatograms(groups, pid)
            if chrom is not None:
                chroms.append(chrom)
        return chroms

    def _get_xcms_chromatograms(self, groups, pid):
        selected = groups.get_group(pid)
        rts = self._get_values(selected, 'rt_values')
        mzs = self._get_values(selected, 'mz_values')
        intensities = self._get_values(selected, 'intensity_values')
        assert len(rts) == len(mzs)
        assert len(rts) == len(intensities)
        if len(rts) > 1:
            chrom = EmpiricalChromatogram(rts, mzs, intensities)
        else:
            chrom = None
        return chrom

    def _get_values(self, df, column_name):
        return df[column_name].values


class ChromatogramLoader(object):
    """
    Loads unknown chemicals and chromatogram data exported from the R script.
    Assumes 1 unknown chemical == 1 chromatogram.
    """
    def __init__(self, xcms_output, min_ms1_intensity, min_rt, max_rt):
        self.min_ms1_intensity = min_ms1_intensity
        self.min_rt = min_rt
        self.max_rt = max_rt
        self.observed_chemicals = self._load_xcms_df(xcms_output)
        self.observed_chromatograms = [x.chromatogram for x in self.observed_chemicals]

    def _load_xcms_df(self, df_file):
        """
        Load CSV file of chromatogram information exported by the XCMS script 'process_data.R'
        :param df_file: the input csv file exported by the script (in gzip format)
        :return: a list of Chromatogram objects
        """
        df = pd.read_csv(df_file, compression='gzip')
        peak_ids = df.id.unique()
        groups = df.groupby('id')
        chemicals = []
        print('Processing exported chromatograms')
        for i in range(len(peak_ids)):
            if i % 5000 == 0:
                print(i)
            pid = peak_ids[i]
            chem = self._get_chemical(groups, pid)
            if chem is not None and self._valid_chem(chem):
                chemicals.append(chem)
        print('Loaded %d UnknownChemicals/Chromatograms' % len(chemicals))
        return chemicals

    def sample(self):
        if self.chromatograms != None:
            selected = np.random.choice(len(self.chromatograms), 1)[0]
            return self.chromatograms[selected]
        else:
            NotImplementedError("Functional Chromatograms not implemented here yet")

    def _get_chemical(self, groups, pid):
        """
        Constructs an EmpiricalChromarogram object from groups in the dataframe.
        :param groups: pandas group object, produced from df.groupby('id'), i.e. each group is a set of rows
        grouped by the 'id' column in the dataframe.
        :param pid: the peak id
        :return: an MS1 chromatographic peak object
        """
        selected = groups.get_group(pid)
        mz = self._get_value(selected, 'mz')
        rt = self._get_value(selected, 'rt')
        max_intensity = self._get_value(selected, 'maxo')
        rts = self._get_values(selected, 'rt_values')
        mzs = self._get_values(selected, 'mz_values')
        intensities = self._get_values(selected, 'intensity_values')
        assert len(rts) == len(mzs)
        assert len(rts) == len(intensities)
        if len(rts) > 1:
            chrom = EmpiricalChromatogram(rts, mzs, intensities)
            chem = UnknownChemical(mz, rt, max_intensity, chrom, None)
        else:
            chem = None
        return chem

    def _valid_chem(self, chem):
        if chem.max_intensity < self.min_ms1_intensity:
            return False
        elif chem.rt < self.min_rt:
            return False
        elif chem.rt > self.max_rt:
            return False
        return True

    def _get_value(self, df, column_name):
        return self._get_values(df, column_name)[0]

    def _get_values(self, df, column_name):
        return df[column_name].values
