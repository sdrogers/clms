import pickle
from sklearn.externals import joblib
import numpy as np
import math

# some useful constants
MZ = 'mz'
INTENSITY = 'intensity'
MZ_INTENSITY = MZ + '_' + INTENSITY
RT = 'rt'
N_PEAKS = 'n_peaks'


def save_obj(obj, filename, use_joblib=True):
    """
    Save object to file
    :param obj: the object to save
    :param filename: the output file
    :param sklearn: if true, use joblib to dump the model
    :return: None
    """
    with open(filename, 'wb') as f:
        if use_joblib:
            # Supposedly using joblib works better for numpy array, see:
            # - https://scikit-learn.org/stable/modules/model_persistence.html
            # - http://gael-varoquaux.info/programming/new_low-overhead_persistence_in_joblib_for_big_data.html
            joblib.dump(obj, f, compress=3)
        else:
            pickle.dump(obj, f)


def load_obj(filename, use_joblib=True):
    """
    Load saved object from file
    :param filename: The file to load
    :param joblib: If true, use joblib to load the model
    :return: the loaded object
    """
    with open(filename, 'rb') as f:
        if use_joblib:
            return joblib.load(f)
        else:
            return pickle.load(f)



# Formulae C12H6N2

# this needs completely redesigning - is basically a dummy function in a workable format
class Formula(object):

    def __init__(self, formula_string, mz, max_intensity, isotope_proportions):
        self.formula_string = formula_string
        self.mz = mz # remove eventually
        self.max_intensity = max_intensity # remove eventually
        self.isotope_proportions = isotope_proportions # remove eventually     
        
    def get_mz_peaks(self):
        # update this to work properly
        peaks = []
        for i in range(len(self.mz)):
            peaks.extend([(self.mz[i], self.max_intensity * self.isotope_proportions[i])])
        return peaks
    # outputs [[(mz_1,max_intensity_1),...,(mz_n,max_intensity_n)], isoptopes, isoptope_proportions]
    
    def get_names(self):
        return ["M+H","Other"]
    
    def get_proportions(self):
        return self.isotope_proportions

class Chemical(object):
    
    def __repr__(self):
        raise NotImplementedError()

    def get_mz_peaks(self, rt, ms_level, isolation_windows):
        raise NotImplementedError()
        
    def get_mz_peaks(self, query_rt, ms_level, isolation_windows):
        if ms_level == 1:
            if not self._rt_match(query_rt):
                return None
        mz_peaks = []
        for which_isotope in range(len(self.isotopes)):
            if ms_level ==1 and self.ms_level == 1:
                if self._isolation_match(query_rt, isolation_windows[0], which_isotope): # check just first set of windows
                    intensity = self._get_intensity(query_rt, which_isotope)
                    mz = self._get_mz(query_rt, which_isotope)
                    return [(mz, intensity)]
            elif ms_level == self.ms_level:
                intensity = self._get_intensity(query_rt, which_isotope)
                mz = self._get_mz(query_rt, which_isotope)
                return [(mz, intensity)]
            else:
                if self._isolation_match(query_rt, isolation_windows[self.ms_level-1], which_isotope) and self.children != None:
                    for i in range(len(self.children)):
                        mz_peaks.extend(self.children[i].get_mz_peaks(query_rt,ms_level,isolation_windows))
                else:
                    return []
        if mz_peaks == []:
            return None
        else:
            return mz_peaks
        
    def _rt_match(self, query_rt):
        if self.ms_level == 1:
            if self.chromatogram._rt_match(query_rt - self.rt) == True:
                return True
            else:
                return False
        else:
            True
        
    def _get_intensity(self, query_rt, which_isotope):
        if self.ms_level == 1:
            return self.isotopes[which_isotope][1] * self.chromatogram.get_relative_intensity(query_rt - self.rt)
        else:
            return (self.parent._get_intensity(query_rt, which_isotope) * self.parent_mass_prop)

    def _get_mz(self, query_rt, which_isotope):
        if self.ms_level == 1:
            return (self.isotopes[which_isotope][0] + self.chromatogram.get_relative_mz(query_rt - self.rt))
        else:
            return self.isotopes[0][0]
            
    def _isolation_match(self, query_rt, isolation_windows, which_isotope):                        
        # assumes list is formated like:
        # [(min_1,max_1),(min_2,max_2),...]
        for window in isolation_windows:
            if (self._get_mz(query_rt, which_isotope) > window[0] and self._get_mz(query_rt, which_isotope) <= window[1]):
                return True
        return False

class UnknownChemical(Chemical):
    """
    Chemical from an unknown chemical formula
    """
    def __init__(self, mz, rt, max_intensity, chromatogram, children):
        self.isotopes = [(mz, max_intensity)]
        self.rt = rt
        self.chromatogram = chromatogram
        self.children = children
        self.isotope_names = ["M+H"]
        self.isotope_proportions = [1]
        self.ms_level = 1
        
    def __repr__(self):
         return 'UnknownChemical mz=%.4f rt=%.2f max_intensity=%.2f' % (self.isotopes[0][0], self.rt, self.isotopes[0][1])

class KnownChemical(Chemical):
    """
    Chemical from an known chemical formula
    """
    def __init__(self, formula, rt, max_itensity, chromatogram, children):
        self.formula = formula.formula_string
        self.rt = rt
        self.chromatogram = chromatogram
        self.children = children
        self.isotopes = formula.get_mz_peaks()
        self.isotope_names = formula.get_names()
        self.isotope_proportions = formula.get_proportions()
        self.ms_level = 1
    
    def __repr__(self):
         return 'KnownChemical - %r' % (self.formula)


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

def chromatogramDensityNormalisation(rts, intensities):
    """
    Definition to standardise the area under a chromatogram to 1. Returns updated intensities
    """
    area = 0.0
    for rt_index in range(len(rts)-1):
        area += ((intensities[rt_index] + intensities[rt_index + 1])/2) / (rts[rt_index+1] - rts[rt_index])
    new_intensities = [x * (1 / area) for x in intensities]
    return new_intensities

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
       

class MSN(Chemical):
    """
    ms2+ fragments
    """
    def __init__(self, mz, ms_level, parent_mass_prop, children=None, parent= None):
        self.isotopes = [(mz,None)]
        self.ms_level = ms_level
        self.parent_mass_prop = parent_mass_prop
        self.children = children
        self.parent = parent
        
    def __repr__(self):
         return 'MSN Fragment mz=%.4f ms_level=%d' % (self.isotopes[0][0], self.ms_level)
        
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
            
class Compound(object):
    def __init__(self, name, chemical_formula, monisotopic_molecular_weight, smiles, inchi, inchikey):
        self.name = name
        self.chemical_formula = chemical_formula
        self.monisotopic_molecular_weight = monisotopic_molecular_weight
        self.smiles = smiles
        self.inchi = inchi
        self.inchikey = inchikey