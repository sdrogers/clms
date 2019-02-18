import numpy as np
import math
from VMSfunctions.Common import load_obj, Peak, NoisyPeak
from VMSfunctions.DataGenerator import PeakSampler
from VMSfunctions.ChineseRestaurantProcess import *


class Sample_Dataset(object):
    """
    Creates a dataset of ms1 and ms2 peaks

    Args:
        file_location: location of PeakSampler object
        n_ms1: Number of ms1 peaks. Will generate from PeakSampler object if set to 'None'
        n_ms2: Number of ms2 peaks. Will generate from PeakSampler object if set to 'None'
        ms2_mz_noise_sd: standard deviation of mz noise in ms2 peaks
        ms2_intensity_noise_sd: standard deviation of intensity noise in ms2 peaks
        alpha: Chinese Restaurant Process parameter. 'alpha  = math.inf' gives completely independent ms2 peaks
        noise_peaks: 'True' or 'False'. Indicated whether to include noisey peaks
        noise_peaks_intensity: True intensity of noise peaks
        noise_peaks_intensity_sd: standard deviation of intensity of noise peaks
        noise_peaks_mz_sd: standard deviation of mz of noise peaks
        n_noise_peaks: number of noise peaks

    Attributes:
        peaks: list of ms1 and ms2 peaks
        ms2_mz_noise_sd: standard deviation of mz noise in ms2 peaks
        ms2_intensity_noise_sd: standard deviation of intensity noise in ms2 peaks
    """

    def __init__(self, file_location, rt=None, n_ms1=None, n_ms2=None, ms2_mz_noise_sd=0, ms2_intensity_noise_sd=0,
                 alpha=math.inf, noise_peaks=False, noise_peaks_intensity=None, noise_peaks_intensity_sd=None,
                 noise_peaks_mz_sd=None, rt_length=None, n_noise_peaks=None):
        ps = load_obj(file_location)
        self.peaks = []
        n_ms1_internal = ps.sample_n(1, n_ms1)
        crp_samples = []
        crp_index = []
        for i in range(n_ms1_internal):
            ms1 = ps.sample(ms_level=1, n_peaks=1, ms2_mz_noise_sd=ms2_mz_noise_sd,
                            ms2_intensity_noise_sd=ms2_intensity_noise_sd)
            if rt != None:
                ms1[0].rt = rt
            self.peaks.extend(ms1)
            n_ms2_internal = ps.sample_n(2, n_ms2)
            len_crp_current = 0
            for j in range(0, n_ms2_internal):
                next_crp = Restricted_Crp(alpha, crp_index, len_crp_current)
                len_crp_current += 1
                crp_index.append(next_crp)
                if next_crp == max(crp_index):
                    ms2 = ps.sample(ms_level=2, n_peaks=1, ms2_mz_noise_sd=ms2_mz_noise_sd,
                                    ms2_intensity_noise_sd=ms2_intensity_noise_sd)
                    ms2[0].parent = ms1[0]
                    if rt != None:
                        ms2[0].rt = rt
                    crp_samples.append(ms2)
                self.peaks.extend(crp_samples[next_crp])
        self.ms2_mz_noise_sd = ms2_mz_noise_sd
        self.ms2_intensity_noise_sd = ms2_intensity_noise_sd
        if noise_peaks == True:
            ms1_noise_peak = ps.sample_noise_peak(ms_level=1, noise_peaks_intensity=noise_peaks_intensity,
                                                  noise_peaks_intensity_sd=noise_peaks_intensity_sd,
                                                  noise_peaks_mz_sd=noise_peaks_mz_sd,
                                                  noise_peaks_rt_length=noise_peaks_rt_length,
                                                  n_noise_peaks=n_noise_peaks)
            if rt != None:
                for j in range(len(ms1_noise_peak)):
                    ms1[j].rt = rt
            self.peaks.extend(ms1_noise_peak)


class Dataset_Scan(object):
    """
    Scans dataset of ms1 and peaks, returning the appropriate output
    """

    def __init__(self, dataset, ms_level, rt, isolation_windows):
        self.mz_in_scan = []
        self.scan_intensities = []  # record intensities here
        true_mz = []
        true_mz_location = []
        for peak in range(0, len(dataset.peaks)):
            mz_internal = dataset.peaks[peak].get(ms_level, rt, isolation_windows)
            if mz_internal != None:
                if dataset.peaks[peak].mz in true_mz:
                    self.scan_intensities[true_mz_location[true_mz.index(dataset.peaks[peak].mz)]] += mz_internal[1]
                else:
                    self.mz_in_scan.append(mz_internal[0])
                    self.scan_intensities.append(mz_internal[1])
                    true_mz.append(dataset.peaks[peak].mz)
                    true_mz_location.append(len(self.scan_intensities) - 1)
        if self.mz_in_scan != []:
            self.mz_in_scan = np.concatenate(self.mz_in_scan)
            self.scan_intensities = np.concatenate(self.scan_intensities)


class Dia_Windows(object):
    """
    Create DIA window design
    """

    def __init__(self, dataset, dia_design, window_type, kaufmann_design, extra_bins, range_slack=0.01,
                 ms1_range=[(None, None)], num_windows=None):
        # set the ms1_range if not given
        self.ms1_values = []
        if ms1_range == [(None, None)]:
            ms1_range = [(math.inf, -math.inf)]
            for dataset_index in range(0, len(dataset.peaks)):
                if dataset.peaks[dataset_index].ms_level == 1:
                    ms1_range = [
                        (min(ms1_range[0][0], dataset.peaks[dataset_index].mz),
                         max(ms1_range[0][1], dataset.peaks[dataset_index].mz))
                    ]
                    self.ms1_values.append(dataset.peaks[dataset_index].mz)
            ms1_range_difference = ms1_range[0][1] - ms1_range[0][0]
        else:
            ms1_range_difference = ms1_range[0][1] - ms1_range[0][0]
            # set the number of windows for kaufmann method
        if dia_design == "kaufmann":
            num_windows = 64
        # dont allow extra bins for basic method    
        if dia_design == "basic" and extra_bins > 0:
            sys.exit("Cannot have extra bins with 'basic' dia design.")
        # find bin walls and extra bin walls
        if window_type == "even":
            internal_bin_walls = [ms1_range[0][0]]
            for window_index in range(0, num_windows):
                internal_bin_walls.append(ms1_range[0][0] + ((window_index + 1) / num_windows) * ms1_range_difference)
            internal_bin_walls[0] = internal_bin_walls[0] - range_slack * ms1_range_difference
            internal_bin_walls[-1] = internal_bin_walls[-1] + range_slack * ms1_range_difference
            internal_bin_walls_extra = None
            if extra_bins > 0:
                internal_bin_walls_extra = [min(self.ms1_values)]
                for window_index in range(0, num_windows * (2 ** extra_bins)):
                    internal_bin_walls_extra.append(ms1_range[0][0] + (
                                (window_index + 1) / (num_windows * (2 ** extra_bins))) * ms1_range_difference)
                internal_bin_walls_extra[0] = internal_bin_walls_extra[0] - range_slack * ms1_range_difference
                internal_bin_walls_extra[-1] = internal_bin_walls_extra[-1] + range_slack * ms1_range_difference
        elif window_type == "percentile":
            internal_bin_walls = np.percentile(self.ms1_values,
                                               np.arange(0, 100 + 100 / num_windows, 100 / num_windows)).tolist()
            internal_bin_walls[0] = internal_bin_walls[0] - range_slack * ms1_range_difference
            internal_bin_walls[-1] = internal_bin_walls[-1] + range_slack * ms1_range_difference
            internal_bin_walls_extra = None
            if extra_bins > 0:
                internal_bin_walls_extra = np.percentile(self.ms1_values,
                                                         np.arange(0, 100 + 100 / (num_windows * (2 ** extra_bins)),
                                                                   100 / (num_windows * (2 ** extra_bins)))).tolist()
                internal_bin_walls_extra[0] = internal_bin_walls_extra[0] - range_slack * ms1_range_difference
                internal_bin_walls_extra[-1] = internal_bin_walls_extra[-1] + range_slack * ms1_range_difference
        else:
            sys.exit("Incorrect window_type selected. Must be 'even' or 'percentile'.")
            # convert bin walls and extra bin walls into locations to scan
        if dia_design == "basic":
            self.locations = []
            for window_index in range(0, num_windows):
                self.locations.append([(internal_bin_walls[window_index], internal_bin_walls[window_index + 1])])
        elif dia_design == "kaufmann":
            self.locations = Kaufmann_Windows(internal_bin_walls, internal_bin_walls_extra, kaufmann_design,
                                              extra_bins).locations
        else:
            sys.ext("Incorrect dia_design selected. Must be 'basic' or 'kaufmann'.")
        # calculate bin walls
        self.bin_walls = list(set(np.array(sum(self.locations, [])).flatten()))
        self.bin_walls.sort()


class Dia_Methods(object):
    """
    Method for doing DIA on a dataset of ms1 and ms2 peaks. Creates windows and then return attributes of scan results
    """

    def __init__(self, dataset, ms_level, rt, dia_design, window_type, kaufmann_design=None, extra_bins=0,
                 range_slack=0.01, ms1_range=[(None, None)], num_windows=None):
        dia_windows = Dia_Windows(dataset, dia_design, window_type, kaufmann_design, extra_bins, range_slack, ms1_range,
                                  num_windows)
        self.bin_walls = dia_windows.bin_walls
        self.locations = dia_windows.locations
        self.ms1_values = dia_windows.ms1_values
        self.mz_in_scans = []
        self.intensities_in_scans = []
        for window_index in range(0, len(self.locations)):
            data_scan = Dataset_Scan(dataset, ms_level, rt, self.locations[window_index])
            self.mz_in_scans.append(np.array(data_scan.mz_in_scan))
            self.intensities_in_scans.append(np.array(data_scan.scan_intensities))


class Dia_Methods_Subsample(object):
    """
    Method for taking a sumsample of DIA results. Helpful for visualising results as scans progress
    """

    def __init__(self, dia_methods, num_scans):
        self.bin_walls = list(set(np.array(sum(dia_methods.locations[0:num_scans], [])).flatten()))
        self.bin_walls.sort()
        self.locations = dia_methods.locations[0:num_scans]
        self.ms1_values = dia_methods.ms1_values
        self.mz_in_scans = dia_methods.mz_in_scans[0:num_scans]
        self.intensities_in_scans = dia_methods.intensities_in_scans[0:num_scans]


class Kaufmann_Windows(object):
    """
    Method for creating window desigs based on Kaufmann paper - https://www.ncbi.nlm.nih.gov/pubmed/27188447 
    """

    def __init__(self, bin_walls, bin_walls_extra, kaufmann_design, extra_bins=0):
        self.locations = []
        if kaufmann_design == "nested":
            n_locations_internal = 4
            for i in range(0, 8):
                self.locations.append([(bin_walls[(0 + i * 8)], bin_walls[(8 + i * 8)])])
        elif kaufmann_design == "tree":
            n_locations_internal = 3
            self.locations.append([(bin_walls[0], bin_walls[32])])
            self.locations.append([(bin_walls[32], bin_walls[64])])
            self.locations.append([(bin_walls[16], bin_walls[48])])
            self.locations.append([(bin_walls[8], bin_walls[24]), (bin_walls[40], bin_walls[56])])
        else:
            print("not a valid design")
        locations_internal = [[] for i in range(n_locations_internal + extra_bins)]
        for i in range(0, 4):
            locations_internal[0].append((bin_walls[(4 + i * 16)], bin_walls[(12 + i * 16)]))
            locations_internal[1].append((bin_walls[(2 + i * 16)], bin_walls[(6 + i * 16)]))
            locations_internal[1].append((bin_walls[(10 + i * 16)], bin_walls[(14 + i * 16)]))
            locations_internal[2].append((bin_walls[(1 + i * 16)], bin_walls[(3 + i * 16)]))
            locations_internal[2].append((bin_walls[(9 + i * 16)], bin_walls[(11 + i * 16)]))
            if kaufmann_design == "nested":
                locations_internal[3].append((bin_walls[(5 + i * 16)], bin_walls[(7 + i * 16)]))
                locations_internal[3].append((bin_walls[(13 + i * 16)], bin_walls[(15 + i * 16)]))
            else:
                locations_internal[2].append((bin_walls[(5 + i * 16)], bin_walls[(7 + i * 16)]))
                locations_internal[2].append((bin_walls[(13 + i * 16)], bin_walls[(15 + i * 16)]))
        if extra_bins > 0:
            for j in range(extra_bins):
                for i in range(64 * (2 ** j)):
                    locations_internal[n_locations_internal + j].append((bin_walls_extra[int(
                        0 + i * ((2 ** extra_bins) / (2 ** j)))], bin_walls_extra[int(
                        ((2 ** extra_bins) / (2 ** j)) / 2 + i * ((2 ** extra_bins) / (2 ** j)))]))
        self.locations.extend(locations_internal)


class Scan_Results_Calculator(object):
    """
    Method for taking raw results, grouping ms2 fragments and determining in which scans they were found
    """

    def __init__(self, dia_results, ms2_mz_slack=0.00001, ms2_intensity_slack=0.1):
        self.intensities_in_scans = dia_results.intensities_in_scans
        self.mz_in_scans = dia_results.mz_in_scans
        self.locations = dia_results.locations
        self.bin_walls = dia_results.bin_walls
        self.ms1_values = dia_results.ms1_values
        self.results = [[] for i in range(len(dia_results.locations))]
        unlisted_mz_in_scans = np.concatenate(dia_results.mz_in_scans)
        unlisted_intensities_in_scans = np.concatenate(dia_results.intensities_in_scans)
        # find unique mz
        unique_mz = [[unlisted_mz_in_scans[0]]]
        unique_intensities = [[unlisted_intensities_in_scans[0]]]
        for unlisted_mz_index in range(1, len(unlisted_mz_in_scans)):
            unique_mz_min = math.inf
            for unique_mz_index in range(len(unique_mz)):
                unique_dist = abs(
                    sum(unique_mz[unique_mz_index]) / len(unique_mz[unique_mz_index]) - unlisted_mz_in_scans[
                        unlisted_mz_index])
                if (unique_dist < unique_mz_min):
                    unique_mz_min = unique_dist
                    unique_mz_which = unique_mz_index
            if unique_mz_min < ms2_mz_slack:
                unique_mz[unique_mz_which].append(unlisted_mz_in_scans[unlisted_mz_index])
                unique_intensities[unique_mz_which].append(unlisted_intensities_in_scans[unlisted_mz_index])
            else:
                unique_mz.append([unlisted_mz_in_scans[unlisted_mz_index]])
                unique_intensities.append([unlisted_intensities_in_scans[unlisted_mz_index]])
        self.ms2_intensities = unique_intensities
        self.ms2_mz = unique_mz
        # find where intensities are unique and assign them a scan result
        for unique_mz_index in range(len(unique_mz)):
            if max(abs(unique_intensities[0] - sum(unique_intensities[0]) / len(
                    unique_intensities[0]))) > ms2_intensity_slack:
                print("not ready yet")
            else:
                for location_index in range(len(dia_results.locations)):
                    TF_in_location = []
                    for unique_index in range(len(unique_mz[unique_mz_index])):
                        TF_in_location.append(
                            unique_mz[unique_mz_index][unique_index] in dia_results.mz_in_scans[location_index] and
                            unique_intensities[unique_mz_index][unique_index] in dia_results.intensities_in_scans[
                                location_index])
                    if any(TF_in_location):
                        self.results[location_index].append(1)
                    else:
                        self.results[location_index].append(0)


class Dia_Location_Finder(object):
    """
    Method for finding location of ms2 fragments based on which DIA scans they are seen in
    """

    def __init__(self, scan_results):
        self.locations = scan_results.locations
        self.results = scan_results.results
        self.bin_walls = scan_results.bin_walls
        self.ms1_values = scan_results.ms1_values
        self.ms2_intensities = scan_results.ms2_intensities
        self.ms2_mz = scan_results.ms2_mz
        self.location_all = []
        bin_mid_points = list((np.array(self.bin_walls[1:]) + np.array(self.bin_walls[:(len(self.bin_walls) - 1)])) / 2)
        for sample_index in range(0, len(self.results[0])):
            mid_point_TF = []
            for mid_points_index in range(0, len(bin_mid_points)):
                mid_point_TF.append(self._mid_point_in_location(bin_mid_points[mid_points_index], sample_index))
            self.location_all.append([(list(np.array(self.bin_walls)[np.where(np.array(mid_point_TF) == True)])[0],
                                       list(np.array(self.bin_walls[1:])[np.where(np.array(mid_point_TF) == True)])[
                                           0])])

    def _mid_point_in_location(self, mid_point, sample_index):
        for locations_index in range(0, len(self.locations)):
            if self._in_window(mid_point, self.locations[locations_index]) == True and self.results[locations_index][
                sample_index] == 0:
                return False
            if self._in_window(mid_point, self.locations[locations_index]) == False and self.results[locations_index][
                sample_index] == 1:
                return False
        else:
            return True

    def _in_window(self, mid_point, locations):
        for window in locations:
            if (mid_point > window[0] and mid_point <= window[1]):
                return True
        return False


class Entropy(object):
    """
    Method for calculating entropy based on locations of ms2 components
    """

    def __init__(self, dia_location_finder):
        self.entropy = []
        self.components_determined = []
        ms1_vec = []
        ms2_vec = []
        for i in range(0, len(dia_location_finder.bin_walls) - 1):
            ms2_vec.extend([0])
            ms1_vec.extend([len(np.where(
                np.logical_and(np.array(dia_location_finder.ms1_values) > dia_location_finder.bin_walls[i],
                               np.array(dia_location_finder.ms1_values) <= dia_location_finder.bin_walls[i + 1]))[0])])
            # fix this
            for j in range(0, len(dia_location_finder.location_all)):
                if [(dia_location_finder.bin_walls[i], dia_location_finder.bin_walls[i + 1])] == \
                        dia_location_finder.location_all[j]:
                    ms2_vec[i] += 1
        ms1_vec_nozero = [value for value in ms1_vec if value != 0]
        ms2_vec_nozero = [value for value in ms1_vec if value != 0]
        entropy_vec = []
        for j in range(0, len(ms2_vec_nozero)):
            entropy_vec.append(-ms2_vec_nozero[j] * ms1_vec_nozero[j] * math.log(1 / ms1_vec_nozero[j]))
        self.entropy = sum(entropy_vec)
        self.components_determined = sum(np.extract(np.array(ms1_vec_nozero) == 1, ms2_vec_nozero))
        self.components = sum(ms2_vec_nozero)


class Entropy_List(object):
    """
    Method for calculating entropy on multiple subsets of the DIA results. Useful for creating plots and monitoring performance over multiple scans
    """

    def __init__(self, dataset, ms_level, rt, dia_design, window_type, kaufmann_design, extra_bins=0, range_slack=0.01,
                 ms1_range=[(None, None)], num_windows=None, ms2_mz_slack=0.00001, ms2_intensity_slack=0.1):
        self.entropy = []
        self.components_determined = []
        if (dia_design != "kaufmann"):
            sys.exit("Only the 'kaufmann' method can be used with Entropy_List")
        if (kaufmann_design == "tree"):
            self.start_subsample_scans = 2
            self.end_subsample_scans = 7 + extra_bins
        elif (kaufmann_design == "nested"):
            self.start_subsample_scans = 8
            self.end_subsample_scans = 12 + extra_bins
        else:
            sys.exit("Cannot use Entropy_List with this design. Kaufmann 'nested' or 'tree' only.")
        for i in range(self.start_subsample_scans, self.end_subsample_scans):
            dia = Dia_Methods_Subsample(
                Dia_Methods(dataset, ms_level, rt, dia_design, window_type, kaufmann_design, extra_bins, range_slack,
                            ms1_range, num_windows), i)
            results = Entropy(Dia_Location_Finder(Scan_Results_Calculator(dia, ms2_mz_slack, ms2_intensity_slack)))
            self.entropy.append(results.entropy)
            self.components_determined.append(results.components_determined)
            self.components = results.components
