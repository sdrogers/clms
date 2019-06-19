import copy
import re
import math

import scipy
import scipy.stats
import glob

from VMSfunctions.ChineseRestaurantProcess import *
from VMSfunctions.Common import *


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
        self.mass = self._get_mz()

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
            return self.formula._get_mz() + float(isotope.split("C13")[0]) * self.mz_diff
        else:
            return None


class Adducts(object):
    def __init__(self, formula, adduct_proportion_cutoff=0.05):
        self.adduct_names = list(POS_TRANSFORMATIONS.keys())
        self.formula = formula
        self.adduct_proportion_cutoff = adduct_proportion_cutoff

    def get_adducts(self):
        adducts = []
        proportions = self._get_adduct_proportions()
        for j in range(len(self.adduct_names)):
            if proportions[j] != 0:
                adducts.extend([(self._get_adduct_names()[j], proportions[j])])
        return adducts

    def _get_adduct_proportions(self):
        # TODO: replace this with something proper
        prior = np.ones(len(self.adduct_names)) * 0.1
        prior[0] = 1.0  # give more weight to the first one, i.e. M+H
        proportions = np.random.dirichlet(prior)
        while max(proportions) < 0.2:
            proportions = np.random.dirichlet(prior)
        proportions[np.where(proportions < self.adduct_proportion_cutoff)] = 0
        proportions = proportions / max(proportions)
        proportions.tolist()
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
        return 'KnownChemical - %r rt=%.2f max_intensity=%.2f' % (
        self.formula.formula_string, self.rt, self.max_intensity)


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


class ChemicalCreator(LoggerMixin):
    def __init__(self, peak_sampler, ROI_sources=None, database=None):
        self.peak_sampler = peak_sampler
        self.ROI_sources = ROI_sources
        self.database = database

    def sample(self, mz_range, rt_range, min_ms1_intensity, n_ms1_peaks, ms_levels, use_database=True, alpha=math.inf,
               fixed_mz=False, adduct_proportion_cutoff=0.05):
        self.mz_range = mz_range
        self.rt_range = rt_range
        self.min_ms1_intensity = min_ms1_intensity
        self.n_ms1_peaks = n_ms1_peaks
        self.ms_levels = ms_levels
        self.use_database = use_database
        self.alpha = alpha
        self.fixed_mz = fixed_mz
        self.adduct_proportion_cutoff = adduct_proportion_cutoff

        # set up some counters
        self.crp_samples = [[] for i in range(self.ms_levels)]
        self.crp_index = [[] for i in range(self.ms_levels)]
        self.counts = [[] for i in range(self.ms_levels)]

        # sample from kernel densities
        if self.ms_levels > 2:
            print("Warning ms_level > 3 not implemented properly yet. Uses scaled ms_level = 2 information for now")
        n_ms1 = self._get_n(1)
        self.logger.debug("{} ms1 peaks to be created.".format(n_ms1))
        sampled_peaks = self.peak_sampler.sample(1, n_ms1, self.mz_range[0][0], self.mz_range[0][1],
                                                 self.rt_range[0][0],
                                                 self.rt_range[0][1], self.min_ms1_intensity)
        # Get formulae from database and check there are enough of them
        if self.use_database == True:
            if len(self.database) < n_ms1 or self.database is None:
                self.logger.warning(
                    'Warning: Not enough compounds in the database. Compounds are not being used. No adducts and isoptopes will be created')
                self.use_database = False
            self.formula_list = self._sample_formulae(sampled_peaks)

        # Get file split information
        split = self._get_n_ROI_files()

        # create chemicals
        chemicals = []
        # load first ROI file
        current_ROI = 0
        ROIs = self._load_ROI_file(current_ROI)
        ROI_intensities = np.array([r.max_intensity for r in ROIs])
        for i in range(n_ms1):
            if i == sum(split[0:(current_ROI+1)]):
                current_ROI += 1
                ROIs = self._load_ROI_file(current_ROI)
                ROI_intensities = np.array([r.max_intensity for r in ROIs])
            if self.use_database == True:
                formula = self.formula_list[i]
            else:
                formula = None
            ROI = ROIs[self._get_ROI_idx(ROI_intensities, sampled_peaks[i].intensity)]
            chem = self._get_chemical(1, formula, ROI, sampled_peaks[i])
            if self.fixed_mz == True:
                chem.mzs = [0 for i in range(len(chem.chromatogram.raw_mzs))]  # not sure how this will work. Not sure why this is set to 0?
            if ms_levels > 1:
                chem.children = self._get_children(1, chem)
            chem.type = CHEM_DATA
            chemicals.append(chem)
            if i % 100 == 0:
                self.logger.debug("i = {}".format(i))
        return chemicals

    def _get_n_ROI_files(self):
        count = 0
        for i in range(len(self.ROI_sources)):
            count += len(glob.glob(self.ROI_sources[i] + '\\*.p'))
        split = np.array([int(np.floor(self.n_ms1_peaks / count)) for i in range(count)])
        split[0:int(self.n_ms1_peaks - sum(split))] += 1
        return split

    def _load_ROI_file(self, file_index):
        num_ROI = 0
        for i in range(len(self.ROI_sources)):
            len_ROI = len(glob.glob(self.ROI_sources[i] + '\\*.p'))
            if len_ROI > file_index:
                ROI_file = glob.glob(self.ROI_sources[i] + '\\*.p')[file_index - num_ROI]
                ROI = load_obj(ROI_file)
                self.logger.debug("Loaded {}".format(ROI_file))
                return ROI
            num_ROI += len(glob.glob(self.ROI_sources[i] + '\\*.p'))

    def _get_ROI_idx(self, ROI_intensities, intensity):
        return (np.abs(ROI_intensities - intensity)).argmin()

    def _sample_formulae(self, sampled_peaks):
        formula_list = []
        compound_mass_list = []
        for index_compound in range(len(self.database)):
            compound_mass_list.append(Formula(self.database[index_compound].chemical_formula).mass)
        sort_index = np.argsort(compound_mass_list)
        compound_mass_list = np.array(compound_mass_list)[sort_index].tolist()
        compound_list = np.array(self.database)[sort_index].tolist()
        for formula_index in range(len(sampled_peaks)):
            mz_peak_sample = sampled_peaks[formula_index].mz
            formula_list.append(compound_list[takeClosest(compound_mass_list, mz_peak_sample)].chemical_formula)
        return formula_list

    def _get_children(self, parent_ms_level, parent, n_peaks=None):  # TODO: this should be moved to the mass spec class
        children_ms_level = parent_ms_level + 1
        if n_peaks is None:
            n_peaks = self._get_n(children_ms_level)
        kids = []
        kids_intensity_proportions = self._get_msn_proportions(children_ms_level, n_peaks)
        if self.alpha < math.inf:
            # draws from here if using Chinese Restaurant Process (SLOW!!!)
            for index_children in range(n_peaks):
                next_crp, self.counts[children_ms_level - 1] = Restricted_Crp(self.alpha,
                                                                              self.counts[children_ms_level - 1],
                                                                              self.crp_index[children_ms_level - 1],
                                                                              index_children)
                self.crp_index[children_ms_level - 1].append(next_crp)
                if next_crp == max(self.crp_index[children_ms_level - 1]):
                    kid = self._get_unknown_msn(children_ms_level, None, None, parent)
                    kid.prop_ms2_mass = kids_intensity_proportions[index_children]
                    if children_ms_level < self.ms_levels:
                        kid.children = self._get_children(children_ms_level, kid)
                    self.crp_samples[children_ms_level - 1].append(kid)
                else:
                    kid = copy.deepcopy(self.crp_samples[children_ms_level - 1][next_crp])
                    kid.parent_mass_prop = self._get_parent_mass_prop(children_ms_level)
                    kid.parent = parent
                kids.append(kid)
            self.crp_samples[children_ms_level - 1].extend(kids)
        else:
            # Draws from here if children all independent
            for index_children in range(n_peaks):
                kid = self._get_unknown_msn(children_ms_level, None, None, parent)
                kid.prop_ms2_mass = kids_intensity_proportions[index_children]
                if children_ms_level < self.ms_levels:
                    kid.children = self._get_children(children_ms_level, kid)
                kids.append(kid)
        return kids

    def _get_msn_proportions(self, children_ms_level, n_peaks):
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
            return int(self.peak_sampler.density_estimator.n_peaks(2, 1))
        else:
            return int(math.floor(self.peak_sampler.density_estimator.n_peaks(2, 1) / (5 ** (ms_level - 2))))

    def _get_chemical(self, ms_level, formula, ROI, sampled_peak):
        # self._get_chemical(1, formula, ROI, sampled_peaks[i])
        if formula != None:
            return self._get_known_ms1(formula, ROI, sampled_peak)
        else:
            return self._get_unknown_msn(ms_level, ROI, sampled_peak)

    def _get_known_ms1(self, formula, ROI, sampled_peak):  # fix this
        rt = sampled_peak.rt
        #len_chrom = ROI.chromatogram.
        #adjusted_rt = rt
        intensity = sampled_peak.intensity
        formula = Formula(formula)
        isotopes = Isotopes(formula)
        adducts = Adducts(formula, self.adduct_proportion_cutoff)
        return KnownChemical(formula, isotopes, adducts, rt, intensity, ROI.chromatogram, None)

    def _get_unknown_msn(self, ms_level, ROI, sampled_peak, parent=None):  # fix this
        if ms_level == 1:
            mz = sampled_peak.mz
            rt = sampled_peak.rt
            intensity = sampled_peak.intensity
            return UnknownChemical(mz, rt, intensity, ROI.chromatogram, None)
        else:
            if ms_level == 2:
                mz = self.peak_sampler.sample(ms_level, 1)[0].mz
            else:
                mz = self.peak_sampler.sample(2, 1)[0].mz
            parent_mass_prop = self._get_parent_mass_prop(ms_level)
            prop_ms2_mass = None
            return MSN(mz, ms_level, prop_ms2_mass, parent_mass_prop, None, parent)

    def _get_parent_mass_prop(self, ms_level):
        return np.random.uniform(0.5, 0.9, 1).tolist()[0]
        # TODO: this needs to come from a density

    def _valid_ms1_chem(self, chem):
        if chem.max_intensity < self.min_ms1_intensity:
            return False
        elif chem.rt < self.min_rt:
            return False
        elif chem.rt > self.max_rt:
            return False
        return True


class MultiSampleCreator(LoggerMixin):

    def __init__(self, original_dataset, n_samples, classes, intensity_noise_sd,
                 change_probabilities, change_differences_means, change_differences_sds, dropout_probabilities=None,
                 experimental_classes=None, experimental_probabilitities=None, experimental_sds=None, save_location=None):
        self.original_dataset = original_dataset
        self.n_samples = n_samples
        self.classes = classes
        self.intensity_noise_sd = intensity_noise_sd
        self.change_probabilities = change_probabilities
        self.change_differences_means = change_differences_means
        self.change_differences_sds = change_differences_sds
        self.dropout_probabilities = dropout_probabilities
        self.experimental_classes = experimental_classes
        self.experimental_probabilitities = experimental_probabilitities
        self.experimental_sds = experimental_sds
        self.save_location = save_location

        self.sample_classes = []
        for index_classes in range(len(self.classes)):
            self.sample_classes.extend([self.classes[index_classes] for i in range(n_samples[index_classes])])
        self.chemical_statuses = self._get_chemical_statuses()
        self.chemical_differences_from_class1 = self._get_chemical_differences_from_class1()
        if self.experimental_classes is not None:
            self.sample_experimental_statuses = self._get_experimental_statuses()
            self.experimental_effects = self._get_experimental_effects()
        self.logger.debug("Classes, Statuses and Differences defined.")

        self.samples = []
        for index_sample in range(sum(self.n_samples)):
            self.logger.debug("Dataset {} of {} created.".format(index_sample + 1, sum(self.n_samples)))
            new_sample = copy.deepcopy(self.original_dataset)
            which_class = np.where(np.array(self.classes) == self.sample_classes[index_sample])
            for index_chemical in range(len(new_sample)):
                if not np.array(self.chemical_statuses)[which_class][0][index_chemical] == "missing":
                    original_intensity = new_sample[index_chemical].max_intensity
                    intensity = self._get_intensity(original_intensity, which_class, index_chemical)
                    adjusted_intensity = self._get_experimental_factor_effect(intensity, index_sample, index_chemical)
                    noisy_adjusted_intensity = self._get_noisy_intensity(adjusted_intensity)
                    new_sample[index_chemical].max_intensity = noisy_adjusted_intensity.tolist()[0]
            chemicals_to_keep = np.where((np.array(self.chemical_statuses)[which_class][0]) != "missing")
            new_sample = np.array(new_sample)[chemicals_to_keep].tolist()
            if self.save_location is not None:
                save_obj(new_sample, self.save_location + '\\sample_%d.p' % index_sample)
            self.samples.append(new_sample)

    def _get_chemical_statuses(self):
        chemical_statuses = [np.array(["unchanged" for i in range(len(self.original_dataset))])]
        chemical_statuses.extend([np.random.choice(["changed", "unchanged"], len(self.original_dataset),
                                                   p=[self.change_probabilities[i], 1 - self.change_probabilities[i]])
                                  for i in range(len(self.classes) - 1)])
        if self.dropout_probabilities is not None:
            for index_chemical in range(len(chemical_statuses)):
                missing = np.where(np.random.binomial(1, self.dropout_probabilities[index_chemical],
                                                      len(chemical_statuses[index_chemical])))
                chemical_statuses[index_chemical][missing] = "missing"
        return chemical_statuses

    def _get_experimental_statuses(self):
        experimental_statuses = []
        for i in range(len(self.experimental_classes)):
            class_allocation = np.random.choice(self.experimental_classes[i], sum(self.n_samples),
                                                p=self.experimental_probabilitities[i])
            experimental_statuses.append(class_allocation)
        return experimental_statuses

    def _get_experimental_effects(self):
        experimental_effects = []
        for i in range(len(self.experimental_classes)):
            coef = [np.random.normal(0, self.experimental_sds[i], len(self.experimental_classes[i])) for j in
                    range(len(self.original_dataset))]
            experimental_effects.append(coef)
        return experimental_effects

    def _get_chemical_differences_from_class1(self):
        chemical_differences_from_class1 = [np.array([0 for i in range(len(self.original_dataset))]) for j in
                                            range(len(self.classes))]
        for index_classes in range(1, len(self.classes)):
            coef_mean = self.change_differences_means[index_classes - 1]
            coef_sd = self.change_differences_sds[index_classes - 1]
            coef_len = sum(self.chemical_statuses[index_classes] == "changed")
            coef = np.random.normal(coef_mean, coef_sd, coef_len)
            chemical_differences_from_class1[index_classes][
                np.where(self.chemical_statuses[index_classes] == "changed")] = coef
        return chemical_differences_from_class1

    def _get_intensity(self, original_intensity, which_class, index_chemical):
        intensity = original_intensity + self.chemical_differences_from_class1[which_class[0][0]][index_chemical]
        return intensity

    def _get_experimental_factor_effect(self, intensity, index_sample, index_chemical):
        experimental_factor_effect = 0.0
        if self.experimental_classes == None:
            return intensity
        else:
            for index_factor in range(len(self.experimental_classes)):
                which_experimental_status = self.sample_experimental_statuses[index_factor][index_sample]
                which_experimental_class = np.where(
                    np.array(self.experimental_classes[index_factor]) == which_experimental_status)
                experimental_factor_effect += self.experimental_effects[index_factor][index_chemical][
                    which_experimental_class]
        return intensity + experimental_factor_effect

    def _get_noisy_intensity(self, adjusted_intensity):
        noisy_intensity = adjusted_intensity + np.random.normal(0, self.intensity_noise_sd[0], 1)
        if noisy_intensity < 0:
            print("Warning: Negative Intensities have been created")
        return noisy_intensity
