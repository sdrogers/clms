import collections
import gzip
import logging
import pickle

from bisect import bisect_left
from sklearn.externals import joblib

# some useful constants
MZ = 'mz'
INTENSITY = 'intensity'
MZ_INTENSITY = MZ + '_' + INTENSITY
RT = 'rt'
N_PEAKS = 'n_peaks'
SCAN_DURATION = 'scan_duration'
POSITIVE = 'positive'
NEGATIVE = 'negative'
DEFAULT_MS1_SCAN_WINDOW = (0, 1e3)

PROTON_MASS = 1.00727645199076


def save_obj(obj, filename):
    """
    Save object to file
    :param obj: the object to save
    :param filename: the output file
    :return: None
    """
    try:
        with gzip.GzipFile(filename, 'w') as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    except OSError:
        logging.getLogger().warning('Old, invalid or missing pickle in %s. Please regenerate this file.' % filename)
        with open(filename, 'wb') as f:
            joblib.dump(obj, f, compress=3)


def load_obj(filename):
    """
    Load saved object from file
    :param filename: The file to load
    :return: the loaded object
    """
    try:
        with gzip.GzipFile(filename, 'rb') as f:
            return pickle.load(f)
    except OSError:
        logging.getLogger().warning('Old, invalid or missing pickle in %s. Please regenerate this file.' % filename)
        with open(filename, 'rb') as f:
            return joblib.load(f)


def chromatogramDensityNormalisation(rts, intensities):
    """
    Definition to standardise the area under a chromatogram to 1. Returns updated intensities
    """
    area = 0.0
    for rt_index in range(len(rts) - 1):
        area += ((intensities[rt_index] + intensities[rt_index + 1]) / 2) / (rts[rt_index + 1] - rts[rt_index])
    new_intensities = [x * (1 / area) for x in intensities]
    return new_intensities

# TODO: add other options
# Note: M+H should come first in this dict because of the prior specification
POS_TRANSFORMATIONS = collections.OrderedDict()
POS_TRANSFORMATIONS['M+H'] = lambda mz : (mz + PROTON_MASS)
POS_TRANSFORMATIONS['[M+ACN]+H'] = lambda mz : (mz + 42.033823)
POS_TRANSFORMATIONS['[M+CH3OH]+H'] = lambda mz : (mz + 33.033489)
POS_TRANSFORMATIONS['[M+NH3]+H'] = lambda mz : (mz + 18.033823)
POS_TRANSFORMATIONS['M+Na'] = lambda mz : (mz + 22.989218)
POS_TRANSFORMATIONS['M+K'] = lambda mz : (mz + 38.963158)
POS_TRANSFORMATIONS['M+2Na-H'] = lambda mz : (mz + 44.971160)
POS_TRANSFORMATIONS['M+ACN+Na'] = lambda mz : (mz + 64.015765)
POS_TRANSFORMATIONS['M+2Na-H'] = lambda mz : (mz + 44.971160)
POS_TRANSFORMATIONS['M+2K+H'] = lambda mz : (mz + 76.919040)
POS_TRANSFORMATIONS['[M+DMSO]+H'] = lambda mz : (mz + 79.02122)
POS_TRANSFORMATIONS['[M+2ACN]+H'] = lambda mz : (mz + 83.060370)
POS_TRANSFORMATIONS['2M+H'] = lambda mz : (mz * 2) + 1.007276
POS_TRANSFORMATIONS['M+ACN+Na'] = lambda mz : (mz + 64.015765)
POS_TRANSFORMATIONS['2M+NH4'] = lambda mz : (mz * 2) + 18.033823


def adduct_transformation(mz, adduct):
    f = POS_TRANSFORMATIONS[adduct]
    return f(mz)


def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return 0
    if pos == len(myList):
        return -1
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return pos
    else:
        return pos - 1


def set_log_level_warning():
    logging.getLogger().setLevel(logging.WARNING)


def set_log_level_info():
    logging.getLogger().setLevel(logging.INFO)


def set_log_level_debug():
    logging.getLogger().setLevel(logging.DEBUG)


# see https://stackoverflow.com/questions/3375443/how-to-pickle-loggers
class LoggerMixin():
    @property
    def logger(self):
        # turn off annoying matplotlib messages
        mpl_logger = logging.getLogger('matplotlib')
        mpl_logger.setLevel(logging.WARNING)
        # initalise basic config for all loggers
        name = "{}".format(type(self).__name__)
        logging.basicConfig(level=logging.getLogger().level)
        logger = logging.getLogger(name)
        return logger