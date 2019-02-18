import pickle
from sklearn.externals import joblib

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


def chromatogramDensityNormalisation(rts, intensities):
    """
    Definition to standardise the area under a chromatogram to 1. Returns updated intensities
    """
    area = 0.0
    for rt_index in range(len(rts) - 1):
        area += ((intensities[rt_index] + intensities[rt_index + 1]) / 2) / (rts[rt_index + 1] - rts[rt_index])
    new_intensities = [x * (1 / area) for x in intensities]
    return new_intensities