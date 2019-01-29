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

class Formula(object):

    def __init__(self, formula_string):
        self.formula_string = formula_string

    def get_isotope_distribution(self, proportion):
        raise NotImplementedError()


class Chemical(object):

    def get_mz_peaks(self, rt, ms_level, isolation_windows):
        raise NotImplementedError()

class UnknownChemical(Chemical):

    def __init__(self, mz, rt, max_intensity, chromatogram):
        self.mz = mz
        self.rt = rt
        self.max_intensity = max_intensity
        self.chromatogram = chromatogram

    def get_mz_peaks(self, query_rt, ms_level, isolation_windows):
        intensity = self.max_intensity * self.chromatogram.get_relative_intensity(query_rt - self.rt)
        mz = self.mz + self.chromatogram.get_mz(query_rt - self.rt)
        return [(mz, intensity)]

class KnownChemical(Chemical):

    def __init__(self, name, formula):
        self.name = name
        self.formula = formula

    def get_mz_peaks(self, rt, ms_level, isolation_windows):


class MassSpectrometer:
    # Make a generator

    def __next__(self):
        raise NotImplementedError()

class IndependentMassSpectrometer(MassSpectrometer):

class ThermoFusionMassSpectrometer:

    def __next__(self):
        raise NotImplementedError()


class Chromatogram:

    def get_relative_intensity(self, rt):
        raise NotImplementedError()

    def get_mz(self, rt):
        raise NotImplementedError()

class Controller:

class DIAController(Controller):

class TopNController(Controller):

class WindowGenerator:

class KaufmannTree:


class PeakGenerator:

class RestrictedCRP(PeakGenerator):






class Peak(object):
    def __init__(self, mz, rt,
                 intensity,                         # the maximum intensity value, 'maxo' in xcms
                 ms_level,                          # the ms-level: 1, 2, ...
                 parent=None,                       # the parent peak, it's another peak object
                 filename=None,                     # the name of the originating file
                 scan_number=None):                 # the scan number where this peak comes from
        self.mz = mz
        self.rt = rt
        self.intensity = intensity
        self.ms_level = ms_level
        self.filename = filename
        self.scan_number = scan_number
        self.parent = parent

    def __repr__(self):
        return 'Peak mz=%.4f rt=%.2f intensity=%.2f ms_level=%d' % (self.mz, self.rt, self.intensity, self.ms_level)

    def get(self, ms_level, rt, isolation_windows):
        if not ms_level == self.ms_level:
            return None
        if not self._rt_match(rt):
            return None
        if ms_level == 1:
            if self._isolation_match(isolation_windows):
                return (self._get_mz(rt), self._get_intensity())
            else:
                return None
        # if we get here, it's ms_level >1 and rt is ok
        if (self.parent._isolation_match(isolation_windows)):
            return (self._get_mz(rt), self._get_intensity(rt))

    def _get_mz(self, rt):
        return self.mz

    def _get_intensity(self, rt):
        return self.intensity

    def _rt_match(self, rt):
        if not rt:
            return True
        if self.ms_level == 1:
            return True  # eventually, some checking here and possibly returning false
        else:
            return self.parent._rt_match(rt)

    def _isolation_match(self, isolation_windows):
        # assumes list is formated like:
        # [(min_1,max_1),(min_2,max_2),...],
        if not isolation_windows:
            return True
        for window in isolation_windows:
            if (self.mz > window[0] and self.mz <= window[1]):
                return True
        return False
        # for i in range(0,len(isolation_window[0])):
        #         if self.parent.mz > isolation_window[0][i] and self.parent.mz <= isolation_window[1][i]: 


class ChromatographicPeak(Peak):
    def __init__(self, mz, rt,
                 intensity,                         # the maximum intensity value, 'maxo' in xcms
                 ms_level,                          # the ms-level: 1, 2, ...
                 parent=None,                       # the parent peak, it's another peak object
                 filename=None,                     # the filename where this peak comes from
                 scan_number=None,                  # the scan number where this peak comes from
                 pid=None,                          # the peak identifier, some string or number
                 sample_idx=None,                   # the index of the originating file in the dataset
                 sn=None,                           # the signal-to-noise ratio
                 mz_range=None,                     # the mz range of the signal
                 rt_range=None,                     # the rt range of the signal
                 integrated_intensity=None,         # the integrated intensity value, 'into' in xcms
                 mz_values=np.array([]),            # the mz values of the signal
                 rt_values=np.array([]),            # the rt values of the signal
                 intensity_values=np.array([])):    # the intensity values of the signal
        super(ChromatographicPeak, self).__init__(mz, rt, intensity, ms_level, parent=parent, filename=filename,
                                                  scan_number=scan_number)
        self.pid = pid
        self.sample_idx = sample_idx
        self.sn = sn
        self.mz_range = mz_range
        self.rt_range = rt_range
        self.integrated_intensity = integrated_intensity
        self.mz_values = mz_values
        self.rt_values = rt_values
        self.intensity_values = intensity_values

    def __repr__(self):
        return 'ChromatographicPeak mz=%.4f rt=%.2f intensity=%.2f ms_level=%d' % (self.mz, self.rt, self.intensity, self.ms_level)

    def get(self, ms_level, rt, isolation_windows):
        raise NotImplementedError

    def _get_mz(self, rt):
        raise NotImplementedError

    def _get_intensity(self, rt):
        raise NotImplementedError

    def _rt_match(self, rt):
        raise NotImplementedError

    def _isolation_match(self, isolation_windows):
        raise NotImplementedError


class NoisyPeak(Peak):
    def __init__(self, ms2_mz_noise_sd, ms2_intensity_noise_sd, mz, rt, intensity, ms_level, parent=None, filename=None, scan_number=None):
        super().__init__(mz, rt, intensity, ms_level, parent=None, filename=None, scan_number=None)
        self.ms2_mz_noise_sd = ms2_mz_noise_sd
        self.ms2_intensity_noise_sd = ms2_intensity_noise_sd
    def get(self, ms_level, rt, isolation_windows):
        if not ms_level == self.ms_level:
            return None
        if not self._rt_match(rt):
            return None
        if ms_level == 1:
            if self._isolation_match(isolation_windows):
                return (self._get_mz(rt), self._get_intensity())
            else:
                return None
        # if we get here, it's ms_level >1 and rt is ok
        if (self.parent._isolation_match(isolation_windows)):
            if self.ms_level==2:
                noisy_mz = self._get_mz(rt) + np.random.normal(0,self.ms2_mz_noise_sd,1)
                noisy_intensity = self._get_intensity(rt) + np.random.normal(0,self.ms2_intensity_noise_sd,1)
                return (noisy_mz, noisy_intensity)   
            else:
                return(self._get_mz(rt), self._get_intensity(rt))
            
class Compound(object):
    def __init__(self, name, chemical_formula, monisotopic_molecular_weight, smiles, inchi, inchikey):
        self.name = name
        self.chemical_formula = chemical_formula
        self.monisotopic_molecular_weight = monisotopic_molecular_weight
        self.smiles = smiles
        self.inchi = inchi
        self.inchikey = inchikey