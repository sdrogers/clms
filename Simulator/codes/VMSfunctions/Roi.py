# Simon's code to extract all regions of interest from an mzML file
# TODO: this code is written for Python 2, need to check that it's actually Python 3 compatible

import collections

import bisect
import numpy as np
import pymzml
from scipy.stats import pearsonr

PROTON_MASS = 1.00727645199076

POS_TRANSFORMATIONS = collections.OrderedDict()
POS_TRANSFORMATIONS['M+H'] = lambda mz: (mz + PROTON_MASS)
POS_TRANSFORMATIONS['[M+ACN]+H'] = lambda mz: (mz + 42.033823)
POS_TRANSFORMATIONS['[M+CH3OH]+H'] = lambda mz: (mz + 33.033489)
POS_TRANSFORMATIONS['[M+NH3]+H'] = lambda mz: (mz + 18.033823)
POS_TRANSFORMATIONS['M+Na'] = lambda mz: (mz + 22.989218)
POS_TRANSFORMATIONS['M+K'] = lambda mz: (mz + 38.963158)
POS_TRANSFORMATIONS['M+2Na-H'] = lambda mz: (mz + 44.971160)
POS_TRANSFORMATIONS['M+ACN+Na'] = lambda mz: (mz + 64.015765)
POS_TRANSFORMATIONS['M+2Na-H'] = lambda mz: (mz + 44.971160)
POS_TRANSFORMATIONS['M+2K+H'] = lambda mz: (mz + 76.919040)
POS_TRANSFORMATIONS['[M+DMSO]+H'] = lambda mz: (mz + 79.02122)
POS_TRANSFORMATIONS['[M+2ACN]+H'] = lambda mz: (mz + 83.060370)
POS_TRANSFORMATIONS['2M+H'] = lambda mz: (mz * 2) + 1.007276
POS_TRANSFORMATIONS['M+ACN+Na'] = lambda mz: (mz + 64.015765)
POS_TRANSFORMATIONS['2M+NH4'] = lambda mz: (mz * 2) + 18.033823


# Object to store a RoI
# Maintains 3 lists -- mz, rt and intensity
# When a new point (mz,rt,intensity) is added, it updates the 
# list and the mean mz which is required.
class Roi(object):
    def __init__(self, mz, rt, intensity):
        self.mz_list = [mz]
        self.rt_list = [rt]
        self.intensity_list = [intensity]
        self.n = 1
        self.mz_sum = mz

    def get_mean_mz(self):
        return self.mz_sum / self.n

    def get_max_intensity(self):
        return max(self.intensity_list)

    def add(self, mz, rt, intensity):
        self.mz_list.append(mz)
        self.rt_list.append(rt)
        self.intensity_list.append(intensity)
        self.mz_sum += mz
        self.n += 1

    def __lt__(self, other):
        # if type(other) == 'Roi':
        return self.get_mean_mz() <= other.get_mean_mz()
    # else:
    # return self.get_mean_mz() <= other


# Find the RoI that a particular mz falls into
# If it falls into nothing, return None
# mz_tol is the window above and below the 
# mean_mz of the RoI. E.g. if mz_tol = 1 Da, then it looks
# plus and minus 1Da
def match(mz, roi_list, mz_tol, mz_units='Da'):
    if len(roi_list) == 0:
        return None
    pos = bisect.bisect_right(roi_list, mz)
    if pos == len(roi_list):
        return None
    if pos == 0:
        return None

    if mz_units == 'Da':
        dist_left = mz.get_mean_mz() - roi_list[pos - 1].get_mean_mz()
        dist_right = roi_list[pos].get_mean_mz() - mz.get_mean_mz()
    else:  # ppm
        dist_left = 1e6 * (mz.get_mean_mz() - roi_list[pos - 1].get_mean_mz()) / mz.get_mean_mz()
        dist_right = 1e6 * (roi_list[pos].get_mean_mz() - mz.get_mean_mz()) / mz.get_mean_mz()

    if dist_left < mz_tol and dist_right > mz_tol:
        return roi_list[pos - 1]
    elif dist_left > mz_tol and dist_right < mz_tol:
        return roi_list[pos]
    elif dist_left < mz_tol and dist_right < mz_tol:
        if dist_left <= dist_right:
            return roi_list[pos - 1]
        else:
            return roi_list[pos]
    else:
        return None


def roi_correlation(roi1, roi2, min_rt_point_overlap=5, method='pearson'):
    # flip around so that roi1 starts earlier (or equal)
    if roi2.rt_list[0] < roi1.rt_list[0]:
        temp = roi2
        roi2 = roi1
        roi1 = temp

    # check that they meet the min_rt_point overlap
    if roi1.rt_list[-1] < roi2.rt_list[0]:
        # no overlap at all
        return 0.0

    # find the position of the first element in roi2 in roi1
    pos = roi1.rt_list.index(roi2.rt_list[0])

    # print roi1.rt_list
    # print roi2.rt_list
    # print pos

    total_length = max([len(roi1.rt_list), len(roi2.rt_list) + pos])
    # print total_length

    r1 = np.zeros((total_length), np.double)
    r2 = np.zeros_like(r1)

    r1[:len(roi1.rt_list)] = roi1.intensity_list
    r2[pos:pos + len(roi2.rt_list)] = roi2.intensity_list

    # print 
    # for i,a in enumerate(r1):
    #     print "{:10.4f}\t{:10.4f}".format(a,r2[i])
    if method == 'pearson':
        r, _ = pearsonr(r1, r2)
    else:
        r = cosine_score(r1, r2)

    return r


def cosine_score(u, v):
    numerator = (u * v).sum()
    denominator = np.sqrt((u * u).sum()) * np.sqrt((v * v).sum())
    return numerator / denominator


# Make the RoI from an input file
# mz_units = Da for Daltons
# mz_units = ppm for ppm
def make_roi(input_file, mz_tol=0.001, mz_units='Da', min_length=10, min_intensity=50000, start_rt=0, stop_rt=10000000):
    # input_file = 'Beer_multibeers_1_fullscan1.mzML'

    if not mz_units == 'Da' and not mz_units == 'ppm':
        print("Unknown mz units, use Da or ppm")
        return None, None

    run = pymzml.run.Reader(input_file, MS1_Precision=5e-6,
                            extraAccessions=[('MS:1000016', ['value', 'unitName'])],
                            obo_version='4.0.1')

    live_roi = []
    dead_roi = []
    junk_roi = []

    for spectrum in run:
        # print spectrum['centroid_peaks']
        if spectrum['ms level'] == 1:
            live_roi.sort()
            current_ms1_scan_rt, units = spectrum['scan start time']
            if units == 'minute':
                current_ms1_scan_rt *= 60.0

            if current_ms1_scan_rt < start_rt:
                continue
            if current_ms1_scan_rt > stop_rt:
                break

            # print current_ms1_scan_rt
            # print spectrum.peaks
            not_grew = set(live_roi)
            for mz, intensity in spectrum.peaks:
                if intensity >= min_intensity:
                    match_roi = match(Roi(mz, 0, 0), live_roi, mz_tol, mz_units=mz_units)
                    if match_roi:
                        match_roi.add(mz, current_ms1_scan_rt, intensity)
                        if match_roi in not_grew:
                            not_grew.remove(match_roi)
                    else:
                        bisect.insort_right(live_roi, Roi(mz, current_ms1_scan_rt, intensity))

            for roi in not_grew:
                if roi.n >= min_length:
                    dead_roi.append(roi)
                else:
                    junk_roi.append(roi)
                pos = live_roi.index(roi)
                del live_roi[pos]

            print("Scan @ {}, {} live ROIs".format(current_ms1_scan_rt, len(live_roi)))

    # process all the live ones - keeping only those that 
    # are longer than the minimum length
    good_roi = dead_roi
    for roi in live_roi:
        if roi.n >= min_length:
            good_roi.append(roi)
        else:
            junk_roi.append(roi)
    return good_roi, junk_roi


def greedy_roi_cluster(roi_list, corr_thresh=0.75, corr_type='cosine'):
    # sort in descending intensity
    roi_list_copy = [r for r in roi_list]
    roi_list_copy.sort(key=lambda x: max(x.intensity_list), reverse=True)
    roi_clusters = []
    while len(roi_list_copy) > 0:
        roi_clusters.append([roi_list_copy[0]])
        remove_idx = [0]
        if len(roi_list_copy) > 1:
            for i, r in enumerate(roi_list_copy[1:]):
                corr = roi_correlation(roi_list_copy[0], r)
                if corr > corr_thresh:
                    roi_clusters[-1].append(r)
                    remove_idx.append(i + 1)
        remove_idx.sort(reverse=True)
        for r in remove_idx:
            del roi_list_copy[r]

    return roi_clusters


if __name__ == '__main__':
    make_roi()
