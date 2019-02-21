import pandas as pd
import scipy
import numpy as np
import sys
import scipy.stats
from VMSfunctions.Common import chromatogramDensityNormalisation, adductTransformation
import re


# Compound was just something I wrote to fill with the stuff I extracted from the HMBD database
# wouldnt go into any of the stuff you are writing
# could potentially be used to fill 'Formula' which could then be used to create Chemicals
class Compound(object):
    def __init__(self, name, chemical_formula, monisotopic_molecular_weight, smiles, inchi, inchikey):
        self.name = name
        self.chemical_formula = chemical_formula
        self.monisotopic_molecular_weight = monisotopic_molecular_weight
        self.smiles = smiles
        self.inchi = inchi
        self.inchikey = inchikey


class Formula(object):
    def __init__(self, formula_string, mz):
        self.formula_string = formula_string
        self.mz = mz # TODO: calculate this later from the formula_string
        
    def _get_mz(self):
        return self.mz
    
    def _get_n_element(self, element):
        if self.formula_string == element:
            return 1
        split_formula = self.formula_string.split(element)
        if(len(split_formula)==1):
            return 0
        for i in range(len(split_formula)):
            if split_formula[i+1][0].isdigit():
                return float(re.split('[A-Z]+',split_formula[i+1])[0])
            else:
                if split_formula[i+1][0].islower():
                    pass
                else:
                    return 1
        return 0
    
class Isotopes(object):
    def __init__(self, formula):
        self.formula = formula
        self.C12_proportion = 0.989
        self.mz_diff = 1.0033548378
        
    def get_isotopes(self, total_proportion):
        # update this to work properly
        peaks = [() for i in range(len(self._get_isotope_proportions(total_proportion)))]
        for i in range(len(peaks)):
            peaks[i] += (self._get_isotope_mz(self._get_isotope_names(i)),)
            peaks[i] += (self._get_isotope_proportions(total_proportion)[i],)
            peaks[i] += (self._get_isotope_names(i),)
        return peaks
    # outputs [(mz_1, intensity_proportion_1, isotope_name_1),...,(mz_n, intensity_proportion_n, isotope_name_n)]
    
    def _get_isotope_proportions(self, total_proportion):
        proportions = [] 
        while sum(proportions) < total_proportion:
            proportions.extend([scipy.stats.binom.pmf(len(proportions),self.formula._get_n_element("C"),1-self.C12_proportion)])        
        normalised_proportions = [proportions[i]/sum(proportions) for i in range(len(proportions))]
        return normalised_proportions
    
    def _get_isotope_names(self, isotope_number):
        if isotope_number == 0:
            return "Mono"
        else:
            return str(isotope_number) + "C13"
    
    def _get_isotope_mz(self, isotope):
        if isotope == "Mono":
            return self.formula._get_mz()
        elif isotope[-3:] == "C13":
            return self.formula._get_mz() - float(isotope.split("C13")[0]) * self.mz_diff
        else:
            return None
            # turn this into a proper function

class Adducts(object):
    def __init__(self, formula):
        self.adduct_names = ["M+H", "[M+ACN]+H", "[M+CH3OH]+H", "[M+NH3]+H"] # remove eventually
        self.formula = formula
        
    def get_adducts(self):
        adducts = []
        proportions = self._get_adduct_proportions()
        for j in range(len(self.adduct_names)):
            if proportions[j] != 0:
                adducts.extend([(self._get_adduct_names()[j], proportions[j])])
        return adducts
   
    def _get_adduct_proportions(self):
        # replace this with something proper
        proportions = np.random.binomial(1,0.1,3) * np.random.uniform(0.1,0.2,3)
        proportions = [1-sum(proportions)] + proportions.tolist()
        return proportions
    
    def _get_adduct_names(self):
        return self.adduct_names

class Chemical(object):
    
    def __repr__(self):
        raise NotImplementedError()
        
    def get_all_mz_peaks(self, query_rt, ms_level, isolation_windows):
        if ms_level == 1:
            if not self._rt_match(query_rt):
                return None
        mz_peaks = []
        for which_isotope in range(len(self.isotopes)):
            for which_adduct in range(len(self._get_adducts())):
                mz_peaks.extend(self._get_mz_peaks(query_rt, ms_level, isolation_windows, which_isotope, which_adduct))
        if mz_peaks == []:
            return None
        else:
            return mz_peaks

    def _get_mz_peaks(self, query_rt, ms_level, isolation_windows, which_isotope, which_adduct):
        mz_peaks = []
        if ms_level ==1 and self.ms_level == 1:
            if not (which_isotope > 0 and which_adduct > 0): # dont give non-mono isotopes adducts
                if self._isolation_match(query_rt, isolation_windows[0], which_isotope, which_adduct): # check just first set of windows
                    intensity = self._get_intensity(query_rt, which_isotope, which_adduct)
                    mz = self._get_mz(query_rt, which_isotope, which_adduct)
                    mz_peaks.extend([(mz, intensity)])
        elif ms_level > 1 and which_isotope > 0:
            pass
        elif ms_level == self.ms_level:
            intensity = self._get_intensity(query_rt, which_isotope, which_adduct)
            mz = self._get_mz(query_rt, which_isotope, which_adduct)
            return [(mz, intensity)]
        else:
            if self._isolation_match(query_rt, isolation_windows[self.ms_level-1], which_isotope, which_adduct) and self.children != None:
                for i in range(len(self.children)):
                    mz_peaks.extend(self.children[i]._get_mz_peaks(query_rt, ms_level, isolation_windows, which_isotope, which_adduct))
            else:
                return []
        return mz_peaks
        
    def _get_adducts(self):
        if self.ms_level == 1:
            return self.adducts
        else:
            return self.parent._get_adducts()
        
    def _rt_match(self, query_rt):
        if self.ms_level == 1:
            if self.chromatogram._rt_match(query_rt - self.rt) == True:
                return True
            else:
                return False
        else:
            True
        
    def _get_intensity(self, query_rt, which_isotope, which_adduct):
        if self.ms_level == 1:
            intensity = self.isotopes[which_isotope][1] * self._get_adducts()[which_adduct][1] * self.max_intensity 
            return (intensity * self.chromatogram.get_relative_intensity(query_rt - self.rt))
        else:
            return (self.parent._get_intensity(query_rt, which_isotope, which_adduct) * self.parent_mass_prop)

    def _get_mz(self, query_rt, which_isotope, which_adduct):
        if self.ms_level == 1:
            return (adductTransformation(self.isotopes[which_isotope][0], self._get_adducts()[which_adduct][0]) + self.chromatogram.get_relative_mz(query_rt - self.rt))
        else:
            return (adductTransformation(self.isotopes[which_isotope][0], self._get_adducts()[which_adduct][0]))
            
    def _isolation_match(self, query_rt, isolation_windows, which_isotope, which_adduct):                        
        # assumes list is formated like:
        # [(min_1,max_1),(min_2,max_2),...]
        for window in isolation_windows:
            if (self._get_mz(query_rt, which_isotope, which_adduct) > window[0] and self._get_mz(query_rt, which_isotope, which_adduct) <= window[1]):
                return True
        return False
    
class UnknownChemical(Chemical):
    """
    Chemical from an unknown chemical formula
    """
    def __init__(self, mz, rt, max_intensity, chromatogram, children = None):
        self.max_intensity = max_intensity
        self.isotopes = [(mz, 1, "Mono")] # [(mz, intensity_proportion, isotope,name)]
        self.adducts = [("M+H",1)]
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
    def __init__(self, formula, isotopes, adducts, rt, max_intensity, chromatogram, children = None, total_proportion = 0.99):
        self.formula = formula
        self.isotopes = isotopes.get_isotopes(total_proportion)
        self.adducts = adducts.get_adducts()
        self.rt = rt
        self.max_intensity = max_intensity
        self.chromatogram = chromatogram
        self.children = children
        self.ms_level = 1
    
    def __repr__(self):
         return 'KnownChemical - %r' % (self.formula.formula_string)
        
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
        self.mzs = [x - sum(mzs)/len(mzs) for x in mzs] # may want to just set this to 0 and remove from input
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
    def __init__(self, type, peak_sampler):
        self.type = type
        self.peak_sampler = peak_sampler
        self.chromatograms = []
        self.chemicals = []

    def get_chromatograms(self):
        return self.chromatograms

    def get_chemicals(self):
        return self.chemicals

    def _sample_peaks(self):
        """
        Samples peaks
        :return: the sampled peaks, where the number of peaks is also drawn from a density
        """
        return self.peak_sampler.sample(ms_level=1)

    def _sample_peak(self, n_peaks):
        """
        Samples just n_peaks from the densities
        :param n_peaks: the number of peaks to generate
        :return: the sampled peaks
        """
        return self.peak_sampler.sample(ms_level=1, n_peaks=n_peaks)


class KnownColumn(Column):
    """
    Represents a chromatographic column of known compounds
    """
    def __init__(self, type, formula_strings, peak_sampler):
        super().__init__(type, peak_sampler)
        for fs in formula_strings:
            p = self._sample_peak(1)[0]
            mz = p.mz # TODO: mz should be calculated from the formula string instead
            rt = p.rt
            max_intensity = p.intensity
            formula = Formula(fs, mz)
            isotopes = Isotopes(formula)
            adducts = Adducts(formula)
            chrom = FunctionalChromatogram("normal", [0, 1])
            chem = KnownChemical(formula, isotopes, adducts, rt, max_intensity, chrom, None)
            self.chemicals.append(chem)
            self.chromatograms.append(chrom)

class UnknownColumn(Column):
    """
    Represents a chromatographic column of unknown chemical compounds that go into the LC-MS instruments.
    """
    def __init__(self, type, num_chemicals, peak_sampler, chromatogram_loader):
        super().__init__(type, peak_sampler)

        # now generates UnknownChemical n-times
        observed_chromatograms = chromatogram_loader.observed_chromatograms
        for i in range(num_chemicals):
            p = self._sample_peak(1)[0]
            chrom = self._sample_chromatogram(observed_chromatograms)
            chem = UnknownChemical(p.mz, p.rt, p.intensity, chrom, None)
            self.chemicals.append(chem)
            self.chromatograms.append(chrom)

    def _sample_chromatogram(self, observed_chromatograms):
        # TODO: sample a chromatogram for this Chemical. For now, we just choose one randomly
        selected = np.random.choice(len(observed_chromatograms), 1)[0]
        return observed_chromatograms[selected]


# controller sends scan request
# mass spec generates scans (is an iterator over scans)
# scan contains: mz list, intensity list, rt, ms_level, precursor_mass, window
# simplest controller: just generates ms1 data

class MassSpectrometer(object):
    # Make a generator
    def __next__(self):
        raise NotImplementedError()

class Scan(object):
    def __init__(self, scan_id, mzs, intensities, ms_level, rt):
        self.scan_id = scan_id

        # ensure that mzs and intensites are sorted by their mz values
        p = mzs.argsort()
        self.mzs = mzs[p]
        self.intensities = intensities[p]

        self.ms_level = ms_level
        self.rt = rt
        assert len(mzs) == len(intensities)
        self.num_peaks = len(mzs)

    def __repr__(self):
        return 'Scan %d -- num_peaks=%d rt=%.2f ms_level=%d' % (self.scan_id, self.num_peaks, self.rt, self.ms_level)

# Independent here refers to how the intensity of each peak in a scan is independent of each other
# i.e. there's no ion supression effect
class IndependentMassSpectrometer(MassSpectrometer):
    def __init__(self, column, scan_times, scan_levels, isolation_windows):
        self.column = column
        self.scan_times = scan_times
        self.scan_levels = scan_levels
        self.isolation_windows = isolation_windows # TODO: isolation windows should be a parameter passed to __next__
        self.chromatograms = column.get_chromatograms()
        self.idx = 0

    def __iter__ (self):
        return self

    def __next__ (self):
        try:
            time = self.scan_times[self.idx]
            ms_level = self.scan_levels[self.idx]
            if ms_level > 1:
                raise NotImplementedError() # TODO: add ms2 support
            scan = self._get_scan(time, ms_level, self.isolation_windows)
        except IndexError:
            raise StopIteration()
        self.idx += 1
        return scan

    def _get_scan(self, scan_time, scan_level, isolation_windows):
        """
        Constructs a scan at a particular timepoint
        :param time: the timepoint
        :return: a mass spectrometry scan at that time
        """
        scan_mzs = []  # all the mzs values in this scan
        scan_intensities = []  # all the intensity values in this scan

        # for all chemicals that come out from the column coupled to the mass spec
        for i in range(len(self.column.chemicals)):
            chemical = self.column.chemicals[i]

            # mzs is a list of (mz, intensity) for the different adduct/isotopes combinations of a chemical
            mzs = chemical.get_all_mz_peaks(scan_time, scan_level, isolation_windows)
            if mzs is not None:
                chem_mzs = [x[0] for x in mzs]
                chem_intensities = [x[1] for x in mzs]
                scan_mzs.extend(chem_mzs)
                scan_intensities.extend(chem_intensities)

        scan_mzs = np.array(scan_mzs)
        scan_intensities = np.array(scan_intensities)
        return Scan(self.idx, scan_mzs, scan_intensities, scan_level, scan_time)


# class ThermoFusionMassSpectrometer:

#     def __next__(self):
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
