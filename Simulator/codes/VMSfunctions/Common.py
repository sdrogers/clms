import pickle
import logging

from sklearn.externals import joblib
from bisect import bisect_left

# some useful constants
MZ = 'mz'
INTENSITY = 'intensity'
MZ_INTENSITY = MZ + '_' + INTENSITY
RT = 'rt'
N_PEAKS = 'n_peaks'
SCAN_DURATION = 'scan_duration'
POSITIVE = 'positive'
NEGATIVE = 'negative'


def save_obj(obj, filename, use_joblib=True):
    """
    Save object to file
    :param obj: the object to save
    :param filename: the output file
    :param use_joblib: if true, use joblib to dump the model
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
    :param use_joblib: If true, use joblib to load the model
    :return: the loaded object
    """
    with open(filename, 'rb') as f:
        if use_joblib:
            return joblib.load(f)
        else:
            return pickle.load(f)


def chromatogramDensityNormalisation(rts, intensities):
    """
    Definition to standardise the area under a chromatogram to 1. Returns updated intensities
    """
    area = 0.0
    for rt_index in range(len(rts) - 1):
        area += ((intensities[rt_index] + intensities[rt_index + 1]) / 2) / (rts[rt_index + 1] - rts[rt_index])
    new_intensities = [x * (1 / area) for x in intensities]
    return new_intensities


def adductTransformation(mz, adduct):
    if adduct == "M+H":
        return (mz - 1.007276)
    elif adduct == "[M+ACN]+H":
        return (mz - 42.03383)
    elif adduct == "[M+CH3OH]+H":
        return (mz - 33.03349)
    elif adduct == "[M+NH3]+H":
        return (mz - 18.03383)    
    else:
        return None
    # turn this into a proper function


def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after
    else:
        return before


def get_logger(name, level=logging.DEBUG):
    # turn off annoying matplotlib messages
    if level == logging.DEBUG:
        mpl_logger = logging.getLogger('matplotlib')
        mpl_logger.setLevel(logging.WARNING)
    # initalise basic config for all loggers
    logging.basicConfig(level=level)
    logger = logging.getLogger(name)
    return logger


def set_log_level_info():
    logging.getLogger().setLevel(logging.INFO)


def set_log_level_debug():
    logging.getLogger().setLevel(logging.DEBUG)
