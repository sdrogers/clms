import glob
import pandas as pd
import os
import numpy as np

from VMSfunctions.Common import *


def get_schedule(n, schedule_dir):
    while True:
        files = sorted(glob.glob(schedule_dir + '\\*.csv'))
        if len(files) == n:
            last_file = files[-1]
            try:
                schedule = pd.read_csv(last_file)
                if schedule.shape[0] == 11951:
                    return last_file
            except:
                pass



class FragmentationPerformance(object):

    def __init__(self, controller_directory, min_acceptable_intensity, controller_file_spec = "*.p", aligned_chemicals=None):
        if aligned_chemicals is not None:
            os.chdir(controller_directory)
            file_names = glob.glob(controller_file_spec)
            n_chemicals_aligned = len(aligned_chemicals["mzmed"])
            self.chemicals_found = [[0 for i in range(n_chemicals_aligned)] for j in range(len(file_names))]
            for controller_index in range(len(file_names)):
                controller = load_obj(file_names[controller_index])
                chems = controller.mass_spec.chemicals
                min_rt = np.array([chems[i].rt for i in range(len(chems))])
                max_rt = np.array([chems[i].rt + chems[i].chromatogram.max_rt for i in range(len(chems))])
                ms2_scan_times = [controller.scans[2][i].rt for i in range(len(controller.scans[2]))]
                ms2_scan_mzs = [controller.scans[2][i].rt for i in range(len(controller.scans[2]))]
                chems_min_mz = [chems[i].isotopes[0][0] + min(chems[i].chromatogram.mzs) for i in range(len(chems))]
                chems_max_mz = [chems[i].isotopes[0][0] + max(chems[i].chromatogram.mzs) for i in range(len(chems))]
                for index_aligned in range(n_chemicals_aligned):
                    min_rt_aligned = aligned_chemicals.iloc[index_aligned, 4]
                    max_rt_aligned = aligned_chemicals.iloc[index_aligned, 5]
                    min_mz_aligned = aligned_chemicals.iloc[index_aligned, 1]
                    max_mz_aligned = aligned_chemicals.iloc[index_aligned, 2]
                    tf1 = np.logical_and(np.greater(max_rt_aligned, min_rt), np.less(min_rt_aligned, max_rt))
                    tf2 = np.logical_and(np.greater(max_mz_aligned, chems_min_mz), np.less(min_mz_aligned, chems_max_mz))
                    idx = np.where(np.logical_and(tf1, tf2))
                    possible_chems = np.array(chems)[idx]
                    found_chemical = False
                    for chem_index in range(len(possible_chems)):
                        for scan_index in range(len(ms2_scan_times)):
                            scan_mz_range = [[(ms2_scan_mzs[scan_index] - controller.isolation_window,
                                               ms2_scan_mzs[scan_index] + controller.isolation_window)]]
                            chem_scan_result = controller.mass_spec._get_all_mz_peaks(possible_chems[chem_index],
                                                                                      ms2_scan_times[scan_index], 1,
                                                                                      scan_mz_range)
                            if chem_scan_result is not None:
                                for scan_result_index in range(len(chem_scan_result)):
                                    if chem_scan_result[scan_result_index][1] > min_acceptable_intensity:
                                        found_chemical = True
                                        self.chemicals_found[controller_index][index_aligned] = 1
                            if found_chemical is True:
                                pass
                        if found_chemical is True:
                            pass
                print("Completed Controller", controller_index + 1)
            self.chemicals_found_total = [sum(sum(np.array(self.chemicals_found)[0:i, :]) > 0) for i in range(1, len(file_names)+1)]
            self.total_matched_chemicals = len(aligned_chemicals["mzmed"])
        else:
            print("Not Implemented")
            self.chemicals_found_total = None
            self.total_matched_chemicals = None
            self.chemicals_found = None