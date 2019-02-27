import numpy as np
import re
import scipy
import scipy.stats

from VMSfunctions.Common import adductTransformation


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
        self.mz = mz  # TODO: calculate this later from the formula_string

    def _get_mz(self):
        return self.mz

    def _get_n_element(self, element):
        if self.formula_string == element:
            return 1
        split_formula = self.formula_string.split(element)
        if (len(split_formula) == 1):
            return 0
        for i in range(len(split_formula)):
            if split_formula[i + 1][0].isdigit():
                return float(re.split('[A-Z]+', split_formula[i + 1])[0])
            else:
                if split_formula[i + 1][0].islower():
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
            proportions.extend(
                [scipy.stats.binom.pmf(len(proportions), self.formula._get_n_element("C"), 1 - self.C12_proportion)])
        normalised_proportions = [proportions[i] / sum(proportions) for i in range(len(proportions))]
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
        self.adduct_names = ["M+H", "[M+ACN]+H", "[M+CH3OH]+H", "[M+NH3]+H"]  # remove eventually
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
        proportions = np.random.binomial(1, 0.1, 3) * np.random.uniform(0.1, 0.2, 3)
        proportions = [1 - sum(proportions)] + proportions.tolist()
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
        if ms_level == 1 and self.ms_level == 1:
            if not (which_isotope > 0 and which_adduct > 0):  # dont give non-mono isotopes adducts
                if self._isolation_match(query_rt, isolation_windows[0], which_isotope,
                                         which_adduct):  # check just first set of windows
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
            if self._isolation_match(query_rt, isolation_windows[self.ms_level - 1], which_isotope, which_adduct)\
                    and self.children != None:
                for i in range(len(self.children)):
                    mz_peaks.extend(self.children[i]._get_mz_peaks(query_rt, ms_level, isolation_windows, which_isotope,
                                                                   which_adduct))
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
            return (adductTransformation(self.isotopes[which_isotope][0],
                                         self._get_adducts()[which_adduct][0]) + self.chromatogram.get_relative_mz(
                query_rt - self.rt))
        else:
            return (adductTransformation(self.isotopes[which_isotope][0], self._get_adducts()[which_adduct][0]))

    def _isolation_match(self, query_rt, isolation_windows, which_isotope, which_adduct):
        # assumes list is formated like:
        # [(min_1,max_1),(min_2,max_2),...]
        for window in isolation_windows:
            if (self._get_mz(query_rt, which_isotope, which_adduct) > window[0] and self._get_mz(query_rt,
                                                                                                 which_isotope,
                                                                                                 which_adduct) <=
                    window[1]):
                return True
        return False


class UnknownChemical(Chemical):
    """
    Chemical from an unknown chemical formula
    """

    def __init__(self, mz, rt, max_intensity, chromatogram, children=None):
        self.max_intensity = max_intensity
        self.isotopes = [(mz, 1, "Mono")]  # [(mz, intensity_proportion, isotope,name)]
        self.adducts = [("M+H", 1)]
        self.rt = rt
        self.chromatogram = chromatogram
        self.children = children
        self.ms_level = 1

    def __repr__(self):
        return 'UnknownChemical mz=%.4f rt=%.2f max_intensity=%.2f' % (
            self.isotopes[0][0], self.rt, self.isotopes[0][1])


class KnownChemical(Chemical):
    """
    Chemical from an known chemical formula
    """

    def __init__(self, formula, isotopes, adducts, rt, max_intensity, chromatogram, children=None,
                 total_proportion=0.99):
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

    def __init__(self, mz, ms_level, parent_mass_prop, children=None, parent=None):
        self.isotopes = [(mz, None, "MSN")]
        self.ms_level = ms_level
        self.parent_mass_prop = parent_mass_prop
        self.children = children
        self.parent = parent

    def __repr__(self):
        return 'MSN Fragment mz=%.4f ms_level=%d' % (self.isotopes[0][0], self.ms_level)


class ChemicalCreator(object):
    def __init__(self, peak_sampler, chromatograms):
        self.peak_sampler = peak_sampler
        self.chromatograms = chromatograms

    def sample(self, min_rt, max_rt, min_ms1_intensity, n_ms1_peaks, ms_levels=2, chemical_type=None, chromatogram_type="Empirical",
               formula_list=None, use_chrom_tuple=False):
        self.n_ms1_peaks = n_ms1_peaks
        self.ms_levels = ms_levels
        self.chemical_type = chemical_type
        self.chromatogram_type = chromatogram_type
        self.formula_list = formula_list
        self.use_chrom_tuple = use_chrom_tuple
        self.min_rt = min_rt
        self.max_rt = max_rt
        self.min_ms1_intensity = min_ms1_intensity
        self.n_ms1_peaks = n_ms1_peaks
        if self.ms_levels > 2:
            print("Warning ms_level > 3 not implemented properly yet. Uses ms_level = 2 information for now")
        n_ms1 = self._get_n(1)
        chemicals = []
        formula = None
        i = 0
        while len(chemicals) < n_ms1:
            sampled_peak = self.peak_sampler.sample(ms_level=1, n_peaks=1)
            chrom = self.chromatograms.sample()
            if self.chemical_type == "Known":
                formula = self.formula_list[i]
            chem = self._get_chemical(1, formula, chrom, sampled_peak[0])
            if chem is not None and self._valid_ms1_chem(chem):
                chem.children = self._get_children(1, chem)
                chemicals.append(chem)
                i += 1
        return chemicals

    # needs to standardise children intensities, such that they add up to parent intensity times scalign factor

    # need to add CRP

    def _get_children(self, parent_ms_level, parent):
        children_ms_level = parent_ms_level + 1
        n_peaks = self._get_n(children_ms_level)
        if n_peaks == None:
            return None
        elif children_ms_level == self.ms_levels:
            kids = []
            for index_children in range(n_peaks):
                kid = self._get_unknown_msn(children_ms_level, None, None, parent)
                kids.append(kid)
            return kids
        elif children_ms_level < self.ms_levels:
            kids = []
            for index_children in range(n_peaks):
                kid = self._get_unknown_msn(children_ms_level, None, None, parent)
                kid._get_children(children_ms_level, kid)
                kids.append()
            return kids
        else:
            return None

    def _get_n(self, ms_level):
        if ms_level == 1:
            return int(self.n_ms1_peaks)
            # TODO: give the option to sample this from a density
        elif ms_level == 2:
            return int(self.peak_sampler.density_estimator.n_peaks(2, 1))  # not sure this will work
        else:
            return int(self.peak_sampler.density_estimator.n_peaks(2, 1))

    def _get_chemical(self, ms_level, formula, chromatogram, sampled_peak):
        if formula != None:
            return self._get_known_ms1(formula, chromatogram, sampled_peak)
        else:
            return self._get_unknown_msn(ms_level, chromatogram, sampled_peak)

    def _get_known_ms1(self, formula, chromatogram, sampled_peak):
        # eventually get rid of mz here
        mz = self._get_mz(1, chromatogram, sampled_peak)
        rt = self._get_rt(chromatogram, sampled_peak)
        intensity = self.get_intensity(chromatogram, sampled_peak)
        formula = Formula(formula, mz)
        isotopes = Isotopes(formula)
        adducts = Adducts(formula)
        chrom = None  # TODO: fix this
        return KnownChemical(formula, isotopes, adducts, rt, intensity, chrom, None)

    def _get_unknown_msn(self, ms_level, chromatogram, sampled_peak, parent=None):
        if ms_level == 1:
            mz = self._get_mz(1, chromatogram, sampled_peak)
            rt = self._get_rt(chromatogram, sampled_peak)
            intensity = self._get_intensity(chromatogram, sampled_peak)
            return UnknownChemical(mz, rt, intensity, chromatogram, None)
        else:
            mz = self._get_mz(ms_level, chromatogram, sampled_peak)
            parent_mass_prop = self._get_parent_prop(ms_level)
            return MSN(mz, ms_level, parent_mass_prop, None, parent)

    def _get_parent_prop(self, ms_level):
        return np.random.uniform(0.2, 0.8, 1).tolist()[0]
        # this needs to come from a density

    def _get_mz(self, ms_level, chromatogram, sampled_peak):
        # not sure what I meant this to do
        if chromatogram == None and sampled_peak == None:
            if ms_level == 2:
                return self.peak_sampler.sample(ms_level, 1)[0].mz
            else:
                return self.peak_sampler.sample(2, 1)[0].mz
        elif self.use_chrom_tuple == False:
            return sampled_peak.mz
        else:
            NotImplementedError()
            # extract same stuff from chromatogram

    def _get_rt(self, chromatogram, sampled_peak):
        if self.use_chrom_tuple == False:
            return sampled_peak.rt
        else:
            NotImplementedError()
            # extract same stuff from chromatogram

    def _get_intensity(self, chromatogram, sampled_peak):
        if self.use_chrom_tuple == False:
            return sampled_peak.intensity
        else:
            NotImplementedError()
            # extract same stuff from chromatogram

    def _valid_ms1_chem(self, chem):
        if chem.max_intensity < self.min_ms1_intensity:
            return False
        elif chem.rt < self.min_rt:
            return False
        elif chem.rt > self.max_rt:
            return False
        return True
