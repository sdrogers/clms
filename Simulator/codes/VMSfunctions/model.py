import pandas as pd
import scipy
import numpy as np
import sys
import scipy.stats
from VMSfunctions.Common import chromatogramDensityNormalisation, aductTransformation

class Formula(object):
    def __init__(self, formula_string, mz):
        self.formula_string = formula_string
        self.mz = mz # calculate this
        
    def _get_mz(self):
        return self.mz
    
    def _get_n_element(self, element):
        return None
    
class Isotopes(object):
    def __init__(self, formula, isotope_proportions, isotope_names):
        self.formula = formula
        self.isotope_proportions = isotope_proportions # remove eventually
        self.isotope_names = isotope_names # remove eventually
        self.mz_diff = 1.0033548378
        
    def get_isotopes(self):
        # update this to work properly
        peaks = []
        for i in range(len(self.isotope_names)):
            peaks.extend([(self._get_isotope_mz(self._get_isotope_names()[i]), self._get_isotope_proportions()[i], self._get_isotope_names()[i])])
        return peaks
    # outputs [(mz_1, intensity_proportion_1, isotope_name_1),...,(mz_n, intensity_proportion_n, isotope_name_n)]
    
    def _get_isotope_proportions(self):
        return self.isotope_proportions
    
    def _get_isotope_names(self):
        return self.isotope_names
    
    def _get_isotope_mz(self, isotope):
        if isotope == "Mono":
            return self.formula._get_mz()
        elif isotope == "1C13":
            return self.formula._get_mz() + self.mz_diff
        elif isotope == "2C13":
            self.formula.get_mz() + 2 * self.mz_diff
        else:
            return None
            # turn this into a proper function
            
class Aducts(object):
    def __init__(self, formula, aduct_names, aduct_proportions):
        self.aduct_names = aduct_names # remove eventually
        self.aduct_proportions = aduct_proportions # remove eventually
        self.formula = formula
        
    def get_aducts(self):
        aducts = []
        for j in range(len(self.aduct_names)):
            aducts.extend([(self._get_aduct_names()[j], self._get_aduct_proportions()[j])])
        return aducts
    
    def _get_aduct_proportions(self):
        return self.aduct_proportions
    
    def _get_aduct_names(self):
        return self.aduct_names

class Chemical(object):
    
    def __repr__(self):
        raise NotImplementedError()
        
    def get_all_mz_peaks(self, query_rt, ms_level, isolation_windows):
        if ms_level == 1:
            if not self._rt_match(query_rt):
                return None
        mz_peaks = []
        for which_isotope in range(len(self.isotopes)):
            for which_aduct in range(len(self._get_aducts()[which_isotope])):
                mz_peaks.extend(self._get_mz_peaks(query_rt, ms_level, isolation_windows, which_isotope, which_aduct))
        if mz_peaks == []:
            return None
        else:
            return mz_peaks

    def _get_mz_peaks(self, query_rt, ms_level, isolation_windows, which_isotope, which_aduct):
        mz_peaks = []
        if ms_level ==1 and self.ms_level == 1:
            if self._isolation_match(query_rt, isolation_windows[0], which_isotope, which_aduct): # check just first set of windows
                intensity = self._get_intensity(query_rt, which_isotope, which_aduct)
                mz = self._get_mz(query_rt, which_isotope, which_aduct)
                mz_peaks.extend([(mz, intensity)])
        elif ms_level > 1 and which_isotope > 0:
            pass
        elif ms_level == self.ms_level:
            intensity = self._get_intensity(query_rt, which_isotope, which_aduct)
            mz = self._get_mz(query_rt, which_isotope, which_aduct)
            return [(mz, intensity)]
        else:
            if self._isolation_match(query_rt, isolation_windows[self.ms_level-1], which_isotope, which_aduct) and self.children != None:
                for i in range(len(self.children)):
                    mz_peaks.extend(self.children[i]._get_mz_peaks(query_rt, ms_level, isolation_windows, which_isotope, which_aduct))
            else:
                return []
        return mz_peaks
        
    def _get_aducts(self):
        if self.ms_level == 1:
            return self.aducts
        else:
            return self.parent._get_aducts()
        
    def _rt_match(self, query_rt):
        if self.ms_level == 1:
            if self.chromatogram._rt_match(query_rt - self.rt) == True:
                return True
            else:
                return False
        else:
            True
        
    def _get_intensity(self, query_rt, which_isotope, which_aduct):
        if self.ms_level == 1:
            intensity = self.isotopes[which_isotope][1] * self._get_aducts()[which_isotope][which_aduct][1] * self.max_intensity
            return (intensity * self.chromatogram.get_relative_intensity(query_rt - self.rt))
        else:
            return (self.parent._get_intensity(query_rt, which_isotope, which_aduct) * self.parent_mass_prop)

    def _get_mz(self, query_rt, which_isotope, which_aduct):
        if self.ms_level == 1:
            return (aductTransformation(self.isotopes[which_isotope][0], self._get_aducts()[which_isotope][which_aduct][0]) + self.chromatogram.get_relative_mz(query_rt - self.rt))
        else:
            return (aductTransformation(self.isotopes[which_isotope][0], self._get_aducts()[which_isotope][which_aduct][0]))
            
    def _isolation_match(self, query_rt, isolation_windows, which_isotope, which_aduct):                        
        # assumes list is formated like:
        # [(min_1,max_1),(min_2,max_2),...]
        for window in isolation_windows:
            if (self._get_mz(query_rt, which_isotope, which_aduct) > window[0] and self._get_mz(query_rt, which_isotope, which_aduct) <= window[1]):
                return True
        return False
    
class UnknownChemical(Chemical):
    """
    Chemical from an unknown chemical formula
    """
    def __init__(self, mz, rt, max_intensity, chromatogram, children):
        self.max_intensity = max_intensity
        self.isotopes = [(mz, 1, "Mono")] # [(mz, intensity_proportion, isotope,name)]
        self.aducts = [[("M+H",1)]]
        self.rt = rt
        self.chromatogram = chromatogram
        self.children = children
        self.ms_level = 1
        
    def __repr__(self):
         return 'UnknownChemical mz=%.4f rt=%.2f max_intensity=%.2f' % (self.isotopes[0][0], self.rt, self.isotopes[0][1])

class KnownChemical(Chemical):
    """
    Chemical from an known chemical formula
    """
    def __init__(self, formula, isotopes, aducts, rt, max_intensity, chromatogram, children):
        self.formula = formula.formula_string
        self.isotopes = isotopes.get_isotopes()
        self.aducts = aducts.get_aducts()
        self.rt = rt
        self.max_intensity = max_intensity
        self.chromatogram = chromatogram
        self.children = children
        self.ms_level = 1
    
    def __repr__(self):
         return 'KnownChemical - %r' % (self.formula)
        
class MSN(Chemical):
    """
    ms2+ fragments
    """
    def __init__(self, mz, ms_level, parent_mass_prop, children=None, parent= None):
        self.isotopes = [(mz,None,"MSN")]
        self.ms_level = ms_level
        self.parent_mass_prop = parent_mass_prop
        self.children = children
        self.parent = parent
        
    def __repr__(self):
         return 'MSN Fragment mz=%.4f ms_level=%d' % (self.isotopes[0][0], self.ms_level)

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
        self.mzs = [x - sum(mzs)/len(mzs) for x in rts] # may want to just set this to 0 and remove from input
        self.intensities = chromatogramDensityNormalisation(rts, intensities)

    def get_relative_intensity(self, query_rt):
        if self._rt_match(query_rt) == False:
            return None
        else:
            return((self.intensities[self._get_rt_neighbours_which(query_rt)[0]] + 
                    (self.intensities[self._get_rt_neighbours_which(query_rt)[1]]
                     -self.intensities[self._get_rt_neighbours_which(query_rt)[0]]) * self._get_distance(query_rt)))
        
    def get_relative_mz(self, query_rt):
        if self._rt_match(query_rt) == False:
            return None
        else:
            return((self.mzs[self._get_rt_neighbours_which(query_rt)[0]] + 
                    (self.mzs[self._get_rt_neighbours_which(query_rt)[1]]
                     -self.mzs[self._get_rt_neighbours_which(query_rt)[0]]) * self._get_distance(query_rt)))
        
    def _get_rt_neighbours(self, query_rt):
        rt_below = max(x for x in self.rts if x <= query_rt)
        rt_above = min(x for x in self.rts if x >= query_rt)
        return([rt_below, rt_above])
    
    def _get_rt_neighbours_which(self, query_rt):
        which_rt_below = self.rts.index(self._get_rt_neighbours(query_rt)[0])
        which_rt_above = self.rts.index(self._get_rt_neighbours(query_rt)[1])
        return([which_rt_below, which_rt_above])
        
    def _get_distance(self, query_rt):
        return((query_rt - self._get_rt_neighbours(query_rt)[0]) / 
               (self._get_rt_neighbours(query_rt)[0] - self._get_rt_neighbours(query_rt)[1]))
    
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
    def __init__(self, distribution, parameters, cutoff = 0.01):
        self.cutoff = cutoff
        self.mz = 0
        if distribution == "normal":
            self.distrib = scipy.stats.norm(parameters[0],parameters[1])
        elif distribution == "gamma":
            self.distrib = scipy.stats.gamma(parameters[0],parameters[1],parameters[2])
        elif distribution == "uniform":
            self.distrib = scipy.stats.uniform(parameters[0],parameters[1])
        else:
            raise NotImplementedError("distribution not implemented")
            
    def get_relative_intensity(self, query_rt):
        if self._rt_match(query_rt) == False:
            return None
        else:
            return(self.distrib.pdf(query_rt + self.distrib.ppf(self.cutoff/2)) * ( 1 / (1 - self.cutoff)))
        
    def get_relative_mz(self, query_rt):
        if self._rt_match(query_rt) == False:
            return None
        else:
            return self.mz

    def _rt_match(self, query_rt):
        if query_rt < 0 or query_rt > self.distrib.ppf(1-(self.cutoff/2)) - self.distrib.ppf(self.cutoff/2):
            return False
        else:
            return True

class Column(object):
    def __init__(self, type, data_file=None):
        self.type = type
        self.chromatograms = []
        if data_file is not None:
            self.chromatogram = self._load_xcms_df(data_file)

    def getChromatograms(self):
        return self.chromatograms

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
        noises = []
        for i in range(len(peak_ids)):
            if i % 5000 == 0:
                print(i)
            pid = peak_ids[i]
            p = self._get_chrom_peak(groups, pid)
            chroms.append(p)
        return chroms

    def _get_chrom_peak(self, groups, pid):
        """
        Constructs a Peak object from groups in the dataframe.
        :param groups: pandas group object, produced from df.groupby('id'), i.e. each group is a set of rows
        grouped by the 'id' column in the dataframe.
        :param pid: the peak id
        :return: an MS1 chromatographic peak object
        """
        selected = groups.get_group(pid)
        rts = self._get_values(selected, 'rt_values')
        mzs = self._get_values(selected, 'mz_values')
        intensities = self._get_values(selected, 'intensity_values')
        ec = EmpiricalChromatogram(rts, mzs, intensities)
        return ec

    def _get_values(self, df, column_name):
        return df[column_name].values


class Compound(object):
    def __init__(self, name, chemical_formula, monisotopic_molecular_weight, smiles, inchi, inchikey):
        self.name = name
        self.chemical_formula = chemical_formula
        self.monisotopic_molecular_weight = monisotopic_molecular_weight
        self.smiles = smiles
        self.inchi = inchi
        self.inchikey = inchikey

# do we need a class which incorporates both KnownChemical and UnknownChemical
# for now would add no noise, but would make future addition of noise easier
# we would need to tell Ronan if we were to do this, so he could adjust how the MassSpectrometer stores Chemicals

# method for generating list chemicals

# class MassSpectrometer:
#     # Make a generator

#     def __next__(self):
#         raise NotImplementedError()

# class IndependentMassSpectrometer(MassSpectrometer):

# class ThermoFusionMassSpectrometer:

#     def __next__(self):
#         raise NotImplementedError()

# class Formula(object):

#     def __init__(self, formula_string):
#         self.formula_string = formula_string

#     def get_isotope_distribution(self, proportion):
#         raise NotImplementedError()

# class Controller:

# class DIAController(Controller):

# class TopNController(Controller):

# # PeakGenerator and RestrictedCRP will either be used within MassSpectrometer or be obsolete
# class PeakGenerator:

# class RestrictedCRP(PeakGenerator):

# class Peak(object):
#     def __init__(self, mz, rt,
#                  intensity,                         # the maximum intensity value, 'maxo' in xcms
#                  ms_level,                          # the ms-level: 1, 2, ...
#                  parent=None,                       # the parent peak, it's another peak object
#                  filename=None,                     # the name of the originating file
#                  scan_number=None):                 # the scan number where this peak comes from
#         self.mz = mz
#         self.rt = rt
#         self.intensity = intensity
#         self.ms_level = ms_level
#         self.filename = filename
#         self.scan_number = scan_number
#         self.parent = parent

#     def __repr__(self):
#         return 'Peak mz=%.4f rt=%.2f intensity=%.2f ms_level=%d' % (self.mz, self.rt, self.intensity, self.ms_level)

#     def get(self, ms_level, rt, isolation_windows):
#         if not ms_level == self.ms_level:
#             return None
#         if not self._rt_match(rt):
#             return None
#         if ms_level == 1:
#             if self._isolation_match(isolation_windows):
#                 return (self._get_mz(rt), self._get_intensity())
#             else:
#                 return None
#         # if we get here, it's ms_level >1 and rt is ok
#         if (self.parent._isolation_match(isolation_windows)):
#             return (self._get_mz(rt), self._get_intensity(rt))

#     def _get_mz(self, rt):
#         return self.mz

#     def _get_intensity(self, rt):
#         return self.intensity

#     def _rt_match(self, rt):
#         if not rt:
#             return True
#         if self.ms_level == 1:
#             return True  # eventually, some checking here and possibly returning false
#         else:
#             return self.parent._rt_match(rt)

#     def _isolation_match(self, isolation_windows):
#         # assumes list is formated like:
#         # [(min_1,max_1),(min_2,max_2),...],
#         if not isolation_windows:
#             return True
#         for window in isolation_windows:
#             if (self.mz > window[0] and self.mz <= window[1]):
#                 return True
#         return False
#         # for i in range(0,len(isolation_window[0])):
#         #         if self.parent.mz > isolation_window[0][i] and self.parent.mz <= isolation_window[1][i]:


# do we also need to write some form of function which extracts peaks from raw data files
# potentially just a wrapper around Joes R function
# the cleaning of the data would need to be fully automated in this case though
# would need this to publish unless we can make extract beer peaks available

# class ChromatographicPeak(Peak):
#     def __init__(self, mz, rt,
#                  intensity,                         # the maximum intensity value, 'maxo' in xcms
#                  ms_level,                          # the ms-level: 1, 2, ...
#                  parent=None,                       # the parent peak, it's another peak object
#                  filename=None,                     # the filename where this peak comes from
#                  scan_number=None,                  # the scan number where this peak comes from
#                  pid=None,                          # the peak identifier, some string or number
#                  sample_idx=None,                   # the index of the originating file in the dataset
#                  sn=None,                           # the signal-to-noise ratio
#                  mz_range=None,                     # the mz range of the signal
#                  rt_range=None,                     # the rt range of the signal
#                  integrated_intensity=None,         # the integrated intensity value, 'into' in xcms
#                  mz_values=np.array([]),            # the mz values of the signal
#                  rt_values=np.array([]),            # the rt values of the signal
#                  intensity_values=np.array([])):    # the intensity values of the signal
#         super(ChromatographicPeak, self).__init__(mz, rt, intensity, ms_level, parent=parent, filename=filename,
#                                                   scan_number=scan_number)
#         self.pid = pid
#         self.sample_idx = sample_idx
#         self.sn = sn
#         self.mz_range = mz_range
#         self.rt_range = rt_range
#         self.integrated_intensity = integrated_intensity
#         self.mz_values = mz_values
#         self.rt_values = rt_values
#         self.intensity_values = intensity_values

#     def __repr__(self):
#         return 'ChromatographicPeak mz=%.4f rt=%.2f intensity=%.2f ms_level=%d' % (self.mz, self.rt, self.intensity, self.ms_level)

#     def get(self, ms_level, rt, isolation_windows):
#         raise NotImplementedError

#     def _get_mz(self, rt):
#         raise NotImplementedError

#     def _get_intensity(self, rt):
#         raise NotImplementedError

#     def _rt_match(self, rt):
#         raise NotImplementedError

#     def _isolation_match(self, isolation_windows):
#         raise NotImplementedError


# class NoisyPeak(Peak):
#     def __init__(self, ms2_mz_noise_sd, ms2_intensity_noise_sd, mz, rt, intensity, ms_level, parent=None, filename=None, scan_number=None):
#         super().__init__(mz, rt, intensity, ms_level, parent=None, filename=None, scan_number=None)
#         self.ms2_mz_noise_sd = ms2_mz_noise_sd
#         self.ms2_intensity_noise_sd = ms2_intensity_noise_sd
#     def get(self, ms_level, rt, isolation_windows):
#         if not ms_level == self.ms_level:
#             return None
#         if not self._rt_match(rt):
#             return None
#         if ms_level == 1:
#             if self._isolation_match(isolation_windows):
#                 return (self._get_mz(rt), self._get_intensity())
#             else:
#                 return None
#         # if we get here, it's ms_level >1 and rt is ok
#         if (self.parent._isolation_match(isolation_windows)):
#             if self.ms_level==2:
#                 noisy_mz = self._get_mz(rt) + np.random.normal(0,self.ms2_mz_noise_sd,1)
#                 noisy_intensity = self._get_intensity(rt) + np.random.normal(0,self.ms2_intensity_noise_sd,1)
#                 return (noisy_mz, noisy_intensity)
#             else:
#                 return(self._get_mz(rt), self._get_intensity(rt))
