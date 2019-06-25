import glob
import pandas as pd
import os
import numpy as np
import tqdm

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


def fragmentation_performance_chemicals(controller_directory, min_acceptable_intensity, controller_file_spec = "*.p"):
    global total_matched_chemicals
    os.chdir(controller_directory)
    file_names = glob.glob(controller_file_spec)
    n_samples = len(file_names)
    controllers = []
    all_chemicals = []
    for controller_index in range(n_samples):
        controller = load_obj(file_names[controller_index])
        controllers.append(controller)
        all_chemicals.extend(controller.mass_spec.chemicals)
    all_rts = [chem.rt for chem in all_chemicals]
    chemicals_found_total = np.unique(all_rts)
    sample_chemical_start_rts = [[] for i in range(n_samples)]
    sample_chemical_start_rts_total = []
    for i in range(n_samples):
        for event in controllers[i].mass_spec.fragmentation_events:
            if event.ms_level == 2:
                if controllers[i].mass_spec._get_intensity(event.chem, event.query_rt, 0,
                                                           0) > min_acceptable_intensity:
                    sample_chemical_start_rts[i].append(event.chem.rt)
        sample_chemical_start_rts[i] = np.unique(np.array(sample_chemical_start_rts[i])).tolist()
        # at this point we have collected the RTs of the all the chemicals that
        # have been fragmented above the min_intensity threshold
        flatten_rts = []
        for l in sample_chemical_start_rts[0:(i + 1)]:
            flatten_rts.extend(l)
        sample_chemical_start_rts_total.append(len(np.unique(np.array(flatten_rts))))
        total_matched_chemicals = sample_chemical_start_rts_total
        print("Completed Controller", i + 1)
    return chemicals_found_total, total_matched_chemicals


def fragmentation_performance_aligned(param_dict):
    controller = load_obj(param_dict["controller_directory"])
    min_acceptable_intensity = param_dict["min_acceptable_intensity"]
    aligned_chemicals = pd.read_csv(param_dict["aligned_chemicals_location"])
    n_chemicals_aligned = len(aligned_chemicals["mzmed"])
    chemicals_found = 0
    for aligned_index in range(n_chemicals_aligned):
        all_relevant_events = []
        for event in controller.mass_spec.fragmentation_events:
            if event.ms_level == 2 and aligned_chemicals["rtmin"][
                aligned_index] < event.query_rt < aligned_chemicals["rtmax"][aligned_index]:
                all_relevant_events.append(event)
        for relevant_event_index in range(len(all_relevant_events)):
            event = controller.mass_spec.fragmentation_events[relevant_event_index]
            mz = controller.mass_spec._get_mz(event.chem, event.query_rt, 0, 0)
            if aligned_chemicals["mzmin"][aligned_index] < mz < aligned_chemicals["mzmax"][
                aligned_index]:
                inten = controller.mass_spec._get_intensity(event.chem, event.query_rt, 0, 0)
                if inten > min_acceptable_intensity:
                    chemicals_found += 1
                    break
    return chemicals_found


def create_frag_dicts(controller_directory, aligned_chemicals_location, min_acceptable_intensity, controller_file_spec="*.p"):
    os.chdir(controller_directory)
    file_names = glob.glob(controller_file_spec)
    params = []
    for controller_index in range(len(file_names)):
        params.append({
            'controller_directory': controller_directory + file_names[controller_index],
            'min_acceptable_intensity': min_acceptable_intensity,
            'aligned_chemicals_location': aligned_chemicals_location
        })
    return params