import numpy as np
import re
import scipy
import scipy.stats
import math
import copy
from random import sample
import sys

from VMSfunctions.ChineseRestaurantProcess import *

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
    def __init__(self, formula_string):
        self.formula_string = formula_string
        self.atom_names = ['C', 'H', 'N', 'O', 'P', 'S', 'Cl', 'I', 'Br', 'Si', 'F', 'D']
        self.atoms = {}
        for atom in self.atom_names:
            self.atoms[atom] = self._get_n_element(atom)

    def _get_mz(self):
        return self.compute_exact_mass()

    def _get_n_element(self, atom_name):
        # Do some regex matching to find the numbers of the important atoms
        ex = atom_name + '(?![a-z])' + '\d*'
        m = re.search(ex, self.formula_string)
        if m == None:
            return 0
        else:
            ex = atom_name + '(?![a-z])' + '(\d*)'
            m2 = re.findall(ex, self.formula_string)
            total = 0
            for a in m2:
                if len(a) == 0:
                    total += 1
                else:
                    total += int(a)
            return total

    def compute_exact_mass(self):
        masses = {'C': 12.00000000000, 'H': 1.00782503214, 'O': 15.99491462210, 'N': 14.00307400524,
                  'P': 30.97376151200, 'S': 31.97207069000, 'Cl': 34.96885271000, 'I': 126.904468, 'Br': 78.9183376,
                  'Si': 27.9769265327, 'F': 18.99840320500, 'D': 2.01410177800}
        exact_mass = 0.0
        for a in self.atoms:
            exact_mass += masses[a] * self.atoms[a]
        return exact_mass

    def __repr__(self):
        return self.formula_string

    def __str__(self):
        return self.formula_string


class Isotopes(object):
    def __init__(self, formula):
        self.formula = formula
        self.C12_proportion = 0.989
        self.mz_diff = 1.0033548378
        # TODO: Add fucntionality for elements other than Carbon

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
        self.adduct_names = ["M+H", "[M+ACN]+H", "[M+CH3OH]+H", "[M+NH3]+H"]  # TODO: remove eventually
        self.formula = formula

    def get_adducts(self):
        adducts = []
        proportions = self._get_adduct_proportions()
        for j in range(len(self.adduct_names)):
            if proportions[j] != 0:
                adducts.extend([(self._get_adduct_names()[j], proportions[j])])
        return adducts

    def _get_adduct_proportions(self):
        # TODO: replace this with something proper
        proportions = np.random.binomial(1, 0.1, 3) * np.random.uniform(0.1, 0.2, 3)
        proportions = [1 - sum(proportions)] + proportions.tolist()
        return proportions

    def _get_adduct_names(self):
        return self.adduct_names


class Chemical(object):

    def __repr__(self):
        raise NotImplementedError()


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
            self.isotopes[0][0], self.rt, self.max_intensity)


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
        return 'KnownChemical - %r rt=%.2f max_intensity=%.2f' % (self.formula.formula_string, self.rt, self.max_intensity)


class MSN(Chemical):
    """
    ms2+ fragments
    """

    def __init__(self, mz, ms_level, prop_ms2_mass, parent_mass_prop, children=None, parent=None):
        self.isotopes = [(mz, None, "MSN")]
        self.ms_level = ms_level
        self.prop_ms2_mass = prop_ms2_mass
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
               formula_list=None, use_chrom_tuple=False, alpha=math.inf, compound_list=None):
        # TODO: Speed this up!
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
        self.crp_samples = [[] for i in range(self.ms_levels)]
        self.crp_index = [[] for i in range(self.ms_levels)]
        self.alpha = alpha
        if self.ms_levels > 2:
            print("Warning ms_level > 3 not implemented properly yet. Uses scaled ms_level = 2 information for now")
        n_ms1 = self._get_n(1)
        print(n_ms1, " ms1 peaks to be created.")
        if self.chemical_type == "Known" and compound_list != None:
            if len(compound_list)<n_ms1:
                sys.exit('compound_list not long enough')
            self.formula_list = []
            compound_sample = sample(range(len(compound_list)),n_ms1)
            for formula_index in range(n_ms1):
                self.formula_list.append(compound_list[compound_sample[formula_index]].chemical_formula)
        formula = None
        chemicals = []
        i = 0
        total = 0
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
            total += 1
            if (i/10 == math.floor(i/10)):
                print("i = ", i, "Total = ", total)
        return chemicals

    def _get_children(self, parent_ms_level, parent):
        children_ms_level = parent_ms_level + 1
        n_peaks = self._get_n(children_ms_level)
        if n_peaks == None:
            return None
        elif children_ms_level == self.ms_levels:
            kids = []
            kids_intensity_proportions = self._get_msn_proportions(children_ms_level, n_peaks)
            for index_children in range(n_peaks):
                next_crp = Restricted_Crp(self.alpha, self.crp_index[children_ms_level-1], index_children)
                self.crp_index[children_ms_level-1].append(next_crp)
                if next_crp == max(self.crp_index[children_ms_level-1]):
                    kid = self._get_unknown_msn(children_ms_level, None, None, parent)
                    kid.prop_ms2_mass = kids_intensity_proportions[index_children]
                    self.crp_samples[children_ms_level-1].append(kid)
                kid = self.crp_samples[children_ms_level-1][next_crp]
                kids.append(kid)
            return kids
        elif children_ms_level < self.ms_levels:
            kids = []
            kids_intensity_proportions = self._get_msn_proportions(children_ms_level, n_peaks)
            for index_children in range(n_peaks):
                next_crp = Restricted_Crp(self.alpha, self.crp_index[children_ms_level-1], index_children)
                self.crp_index[children_ms_level - 1].append(next_crp)
                if next_crp == max(self.crp_index[children_ms_level - 1]):
                    kid = self._get_unknown_msn(children_ms_level, None, None, parent)
                    kid.prop_ms2_mass = kids_intensity_proportions[index_children]
                    kid.children = self._get_children(children_ms_level, kid)
                    self.crp_samples[children_ms_level - 1].append(kid)
                kid = copy.deepcopy(self.crp_samples[children_ms_level - 1][next_crp])
                kid.parent_mass_prop = self._get_parent_mass_prop(children_ms_level)
                kids.append(kid)
            self.crp_samples[children_ms_level-1].extend(kids)
            return kids
        else:
            return None
        
    def _get_msn_proportions(self, children_ms_level,n_peaks):
        if children_ms_level == 2:
            kids_intensities = self.peak_sampler.sample(children_ms_level, n_peaks)
        else:
            kids_intensities = self.peak_sampler.sample(2, n_peaks)
        kids_intensities_total = sum([x.intensity for x in kids_intensities])
        kids_intensities_proportion = [x.intensity / kids_intensities_total for x in kids_intensities]
        return kids_intensities_proportion

    def _get_n(self, ms_level):
        if ms_level == 1:
            return int(self.n_ms1_peaks)
            # TODO: give the option to sample this from a density
        elif ms_level == 2:
            return int(self.peak_sampler.density_estimator.n_peaks(2, 1))  # not sure this will work
        else:
            return int(math.floor(self.peak_sampler.density_estimator.n_peaks(2, 1) / (5**(ms_level-2))))

    def _get_chemical(self, ms_level, formula, chromatogram, sampled_peak):
        if formula != None:
            return self._get_known_ms1(formula, chromatogram, sampled_peak)
        else:
            return self._get_unknown_msn(ms_level, chromatogram, sampled_peak)

    def _get_known_ms1(self, formula, chromatogram, sampled_peak):
        mz = self._get_mz(1, chromatogram, sampled_peak)
        rt = self._get_rt(chromatogram, sampled_peak)
        intensity = self._get_intensity(chromatogram, sampled_peak)
        formula = Formula(formula)
        isotopes = Isotopes(formula)
        adducts = Adducts(formula)
        return KnownChemical(formula, isotopes, adducts, rt, intensity, chromatogram, None)

    def _get_unknown_msn(self, ms_level, chromatogram, sampled_peak, parent=None):
        if ms_level == 1:
            mz = self._get_mz(1, chromatogram, sampled_peak)
            rt = self._get_rt(chromatogram, sampled_peak)
            intensity = self._get_intensity(chromatogram, sampled_peak)
            return UnknownChemical(mz, rt, intensity, chromatogram, None)
        else:
            mz = self._get_mz(ms_level, chromatogram, sampled_peak)
            parent_mass_prop = self._get_parent_mass_prop(ms_level)
            prop_ms2_mass = None
            return MSN(mz, ms_level, prop_ms2_mass, parent_mass_prop, None, parent)

    def _get_parent_mass_prop(self, ms_level):
        return np.random.uniform(0.5, 0.9, 1).tolist()[0]
        # TODO: this needs to come from a density

    def _get_mz(self, ms_level, chromatogram, sampled_peak):
        # TODO: work out what the hell this is
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
            # TODO: extract same stuff from chromatogram

    def _get_rt(self, chromatogram, sampled_peak):
        if self.use_chrom_tuple == False:
            return sampled_peak.rt
        else:
            NotImplementedError()
            # TODO: extract same stuff from chromatogram

    def _get_intensity(self, chromatogram, sampled_peak):
        if self.use_chrom_tuple == False:
            return sampled_peak.intensity
        else:
            NotImplementedError()
            # TODO: extract same stuff from chromatogram

    def _valid_ms1_chem(self, chem):
        if chem.max_intensity < self.min_ms1_intensity:
            return False
        elif chem.rt < self.min_rt:
            return False
        elif chem.rt > self.max_rt:
            return False
        return True
