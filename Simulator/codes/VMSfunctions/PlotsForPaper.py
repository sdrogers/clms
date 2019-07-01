import os
from collections import defaultdict

import numpy as np
import pandas as pd
import pylab as plt
import seaborn as sns

from VMSfunctions.Chemicals import UnknownChemical, get_absolute_intensity, get_key
from VMSfunctions.Common import load_obj, PROTON_MASS
from VMSfunctions.TopNExperiment import get_chemicals, get_precursor_info, get_chem_to_frag_events


def get_N(row):
    if 'T10' in row['filename']:
        return 10
    else:
        return row['filename'].split('_')[3]


def get_dew(row):
    if 'T10' in row['filename']:
        return 15
    else:
        tok = row['filename'].split('_')[5]  # get the dew value in the filename
        return tok.split('.')[0]  # get the part before '.mzML'


def experiment_group(row):
    if 'experiment' in row:
        col_to_check = 'experiment'
    else:
        col_to_check = 'filename'

    if 'beer' in row[col_to_check]:
        return 'beer'
    else:
        return 'urine'


def add_group_column(df):
    df['group'] = df.apply(lambda row: experiment_group(row), axis=1)


def get_df(csv_file, min_ms1_intensity, rt_range, mz_range):
    df = pd.read_csv(csv_file)
    return filter_df(df, min_ms1_intensity, rt_range, mz_range)


def filter_df(df, min_ms1_intensity, rt_range, mz_range):
    # filter by rt range
    if rt_range is not None:
        df = df[(df['rt'] > rt_range[0][0]) & (df['rt'] < rt_range[0][1])]

    # filter by mz range
    if mz_range is not None:
        df = df[(df['rt'] > mz_range[0][0]) & (df['rt'] < mz_range[0][1])]

    # filter by min intensity
    intensity_col = 'maxo'
    if min_ms1_intensity is not None:
        df = df[(df[intensity_col] > min_ms1_intensity)]

    # add log intensity column
    df['log_intensity'] = df.apply(lambda row: np.log(row[intensity_col]), axis=1)

    # add N column
    try:
        df['N'] = df.apply(lambda row: get_N(row), axis=1)
        df[['N']] = df[['N']].astype('int')
    except IndexError:
        pass
    except ValueError:
        df['N'] = df.apply(lambda row: np.nan, axis=1)

    # add group column
    df['group'] = df.apply(lambda row: experiment_group(row), axis=1)
    return df


def make_boxplot(df, x, y, xticklabels, title, outfile=None):
    g = sns.catplot(x=x, y=y, kind='box', data=df)
    g.fig.set_size_inches(10, 3)
    if xticklabels is not None:
        g.set_xticklabels(xticklabels, rotation=90)
    else:
        g.set_xticklabels(rotation=90)
    plt.title(title)
    plt.tight_layout()
    if outfile is not None:
        plt.savefig(outfile, dpi=300)
    plt.show()


def make_hist(df, col_name, file_name, title):
    gb = df.groupby('filename')
    group_df = gb.get_group(file_name)
    vals = group_df[col_name].values
    print(vals, len(vals))
    _ = plt.hist(vals, bins=100)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def to_chemical(row):
    mz = row['mz'] - PROTON_MASS
    rt = row['rt']
    max_intensity = row['maxo']
    chrom = None
    chem = UnknownChemical(mz, rt, max_intensity, chrom, children=None)
    return chem


def df_to_chemicals(df, filename):
    filtered_df = df.loc[df['filename'] == filename]
    chems = filtered_df.apply(lambda row: to_chemical(row), axis=1).values
    return chems


def find_chem(to_find, min_rts, max_rts, min_mzs, max_mzs, chem_list):
    query_mz = to_find.isotopes[0][0]
    query_rt = to_find.rt
    min_rt_check = min_rts <= query_rt
    max_rt_check = query_rt <= max_rts
    min_mz_check = min_mzs <= query_mz
    max_mz_check = query_mz <= max_mzs
    idx = np.nonzero(min_rt_check & max_rt_check & min_mz_check & max_mz_check)[0]
    matches = chem_list[idx]

    # pick a match
    if len(matches) == 0:
        return None
    elif len(matches) == 1:
        return matches[0]
    else:  # multiple matches, take the closest in rt
        diffs = [np.abs(chem.rt - to_find.rt) for chem in matches]
        idx = np.argmin(diffs)
        return matches[idx]


def match(chemical_list_1, chemical_list_2, mz_tol, rt_tol, verbose=False):
    matches = {}
    chem_list = np.array(chemical_list_2)
    min_rts = np.array([chem.rt - rt_tol for chem in chem_list])
    max_rts = np.array([chem.rt + rt_tol for chem in chem_list])
    min_mzs = np.array([chem.isotopes[0][0] * (1 - mz_tol / 1e6) for chem in chem_list])
    max_mzs = np.array([chem.isotopes[0][0] * (1 + mz_tol / 1e6) for chem in chem_list])
    for i in range(len(chemical_list_1)):
        to_find = chemical_list_1[i]
        if i % 1000 == 0 and verbose:
            print('%d/%d found %d' % (i, len(chemical_list_1), len(matches)))
        match = find_chem(to_find, min_rts, max_rts, min_mzs, max_mzs, chem_list)
        if match:
            matches[to_find] = match
    return matches


def check_found_matches(matches, left_label, right_label, N=20):
    found = [key for key in matches if matches[key] is not None]
    print('Found %d/%d (%f)' % (len(found), len(matches), len(found) / len(matches)))

    print('%s\t\t\t\t\t\t%s' % (left_label, right_label))
    for key, value in list(matches.items())[0:N]:
        if value is not None:
            print('mz %.2f rt %.4f intensity %.4f\tmz %.2f rt %.4f intensity %.4f' % (
                key[0], key[1], key[2], value[0], value[1], value[2]))


def load_controller(results_dir, experiment_name, N, rt_tol):
    analysis_name = 'experiment_%s_N_%d_rttol_%d' % (experiment_name, N, rt_tol)
    pickle_in = '%s/%s.p' % (results_dir, analysis_name)
    print('Loading %s' % analysis_name)
    try:
        controller = load_obj(pickle_in)
    except FileNotFoundError:
        controller = None
    return controller


def load_controllers(results_dir, Ns, rt_tols):
    controllers = []
    for N in Ns:
        for rt_tol in rt_tols:
            controller = load_controller(results_dir, N, rt_tol)
            if controller is not None:
                controllers.append(controller)
    return controllers


def compute_performance_scenario_1(controller, dataset, min_ms1_intensity,
                                   fullscan_filename, P_peaks_df,
                                   matching_mz_tol, matching_rt_tol,
                                   chem_to_frag_events=None):
    if chem_to_frag_events is None:  # read MS2 fragmentation events from pickled controller
        chem_to_frag_events = get_frag_events(controller, 2)

    # match with xcms peak-picked ms1 data
    detected_ms1 = df_to_chemicals(P_peaks_df, fullscan_filename)
    matches_fullscan = match(dataset, detected_ms1, matching_mz_tol, matching_rt_tol, verbose=False)

    # check if matched and set a flag to indicate that
    update_matched_status(dataset, matches_fullscan, None)

    # positive instances are ground truth MS1 peaks found by XCMS
    # negative instances are chemicals that cannot be matched to XCMS output
    positives = list(filter(lambda x: x.found_in_fullscan, dataset))
    negatives = list(filter(lambda x: not x.found_in_fullscan, dataset))

    # for both positive and negative instances, count how many frag events they have
    # and whether it's above (good) or below (bad) the minimum ms1 intensity at the time of fragmentation.
    positives_count = get_chem_frag_counts(positives, chem_to_frag_events, min_ms1_intensity)
    negatives_count = get_chem_frag_counts(negatives, chem_to_frag_events, min_ms1_intensity)

    # TP = positive instances that are good only
    tp = [chem for chem in positives if positives_count[chem]['good'] > 0 and positives_count[chem]['bad'] == 0]

    # FP = negative instances that are fragmented (both good + bad)
    fp = [chem for chem in negatives if negatives_count[chem]['good'] > 0 or negatives_count[chem]['bad'] > 0]

    # FN = positive instances that are not fragmented at all + positive instances that are bad only
    fn = [chem for chem in positives if \
          (positives_count[chem]['good'] == 0 and positives_count[chem]['bad'] == 0) or \
          (positives_count[chem]['good'] == 0 and positives_count[chem]['bad'] > 0)]

    tp = len(tp)
    fp = len(fp)
    fn = len(fn)
    prec, rec, f1 = compute_pref_rec_f1(tp, fp, fn)
    return tp, fp, fn, prec, rec, f1


# def compute_performance_scenario_1(controller, chemicals, min_ms1_intensity,
#                                    fullscan_filename, P_peaks_df,
#                                    matching_mz_tol, matching_rt_tol,
#                                    chem_to_frag_events=None):

#     if chem_to_frag_events is None: # read MS2 fragmentation events from pickled controller
#         chem_to_frag_events = get_frag_events(controller, 2)

#     # match xcms picked ms1 peaks to fragmentation peaks
#     detected_ms1 = df_to_chemicals(P_peaks_df, fullscan_filename)
#     matches_fullscan = match(detected_ms1, chemicals, matching_mz_tol, matching_rt_tol, verbose=False)
#     matched_frags = set(matches_fullscan.values())
#     print('%d/%d %d/%d' % (len(matches_fullscan), len(detected_ms1), len(matched_frags), len(chemicals)))

#     # ms1 peaks that are also fragmented
#     positives = []
#     for ms1_peak in matches_fullscan:
#         frag_peak = matches_fullscan[ms1_peak]
#         frag_events = chem_to_frag_events[frag_peak]
#         if len(frag_events) > 0:
#             positives.append(frag_peak)

#     # fragmentation peaks that are not in ms1 peaks
#     negatives = []
#     for frag_peak in chemicals:
#         if frag_peak not in matched_frags:
#             frag_events = chem_to_frag_events[frag_peak]
#             if len(frag_events) > 0:
#                 negatives.append(frag_peak)

#     positives_count = get_chem_frag_counts(positives, chem_to_frag_events, min_ms1_intensity)
#     negatives_count = get_chem_frag_counts(negatives, chem_to_frag_events, min_ms1_intensity)

#     # peaks from ground truth (found in full-scan files) that are fragmented above the minimum intensity threshold
#     tp = [chem for chem in positives if positives_count[chem]['good'] > 0 and positives_count[chem]['bad'] == 0]
#     tp = len(tp)

#     # peaks from ground truth that are not fragmented + peaks from ground truth that are fragmented below the minimum intensity threshold.
#     fp = len(detected_ms1) - tp

#     # peaks not from ground truth that are fragmented above the minimum intensity threshold.
#     fn = [chem for chem in negatives if negatives_count[chem]['good'] > 0 and negatives_count[chem]['bad'] == 0]
#     fn = len(fn)

#     prec, rec, f1 = compute_pref_rec_f1(tp, fp, fn)
#     return tp, fp, fn, prec, rec, f1


def compute_performance_scenario_2(controller, dataset, min_ms1_intensity,
                                   fullscan_filename, fragfile_filename,
                                   P_peaks_df, Q_peaks_df,
                                   matching_mz_tol, matching_rt_tol,
                                   chem_to_frag_events=None):
    if chem_to_frag_events is None:  # read MS2 fragmentation events from pickled controller
        chem_to_frag_events = get_frag_events(controller, 2)

    # load the list of xcms-picked peaks
    detected_from_fullscan = df_to_chemicals(P_peaks_df, fullscan_filename)
    detected_from_fragfile = df_to_chemicals(Q_peaks_df, fragfile_filename)

    # match with xcms peak-picked ms1 data from fullscan file
    matches_fullscan = match(dataset, detected_from_fullscan, matching_mz_tol, matching_rt_tol, verbose=False)

    # match with xcms peak-picked ms1 data from fragmentation file
    matches_fragfile = match(dataset, detected_from_fragfile, matching_mz_tol, matching_rt_tol, verbose=False)

    # check if matched and set a flag to indicate that
    update_matched_status(dataset, matches_fullscan, matches_fragfile)

    # True positive: a peak that is fragmented above the minimum MS1 intensity and is picked by XCMS from
    # the MS1 information in the DDA file and is picked in the fullscan file.
    found_in_both = list(filter(lambda x: x.found_in_fullscan and x.found_in_fragfile, dataset))
    frag_count = get_chem_frag_counts(found_in_both, chem_to_frag_events, min_ms1_intensity)
    tp = [chem for chem in found_in_both if frag_count[chem]['good'] > 0 and frag_count[chem]['bad'] == 0]
    tp = len(tp)

    # False positive: any peak that is above minimum intensity and is picked by XCMS
    # from the DDA file but is not picked from the fullscan.
    found_in_dda_only = list(filter(lambda x: not x.found_in_fullscan and x.found_in_fragfile, dataset))
    frag_count = get_chem_frag_counts(found_in_dda_only, chem_to_frag_events, min_ms1_intensity)
    fp = [chem for chem in found_in_dda_only if frag_count[chem]['good'] > 0 and frag_count[chem]['bad'] == 0]
    fp = len(fp)

    # False negative: any peak that is picked from fullscan data, and is not fragmented, or
    # is fragmented below the minimum intensity.
    found_in_fullscan = list(filter(lambda x: x.found_in_fullscan, dataset))
    fn = len(found_in_fullscan) - tp

    prec, rec, f1 = compute_pref_rec_f1(tp, fp, fn)
    return tp, fp, fn, prec, rec, f1


def get_frag_events(controller, ms_level):
    '''
    Gets the fragmentation events for all chemicals for an ms level from the controller
    :param controller: A Top-N controller object
    :param ms_level: The MS-level (usually 2)
    :return: A dictionary where keys are chemicals and values are a list of fragmentation events
    '''
    filtered_frag_events = list(filter(lambda x: x.ms_level == ms_level, controller.mass_spec.fragmentation_events))
    chem_to_frag_events = defaultdict(list)
    for frag_event in filtered_frag_events:
        key = frag_event.chem
        chem_to_frag_events[key].append(frag_event)
    return dict(chem_to_frag_events)


def count_frag_events(chem, chem_to_frag_events, min_ms1_intensity):
    '''
    Counts how many good and bad fragmentation events for each chemical (key).
    Good fragmentation events are defined as fragmentation events that occur when at the time of fragmentation,
    the chemical MS1 intensity is above the min_ms1_intensity threshold.
    :param chem: the chemical to count
    :param chem_to_frag_events: a dictionary of chemicals to frag events (from get_frag_events above())
    :return: a tuple of good and bad fragmentation event counts
    '''
    frag_events = chem_to_frag_events[chem]
    good_count = 0
    bad_count = 0
    for frag_event in frag_events:
        chem = frag_event.chem
        query_rt = frag_event.query_rt
        if get_absolute_intensity(chem, query_rt) < min_ms1_intensity:
            bad_count += 1
        else:
            good_count += 1
    return good_count, bad_count


def get_chem_frag_counts(chem_list, chem_to_frag_events, min_ms1_intensity):
    # get the count of good/bad fragmentation events for all chemicals in chem_list
    results = {}
    for i in range(len(chem_list)):
        chem = chem_list[i]
        try:
            good_count, bad_count = count_frag_events(chem, chem_to_frag_events, min_ms1_intensity)
        except KeyError:
            good_count = 0
            bad_count = 0
        results[chem] = {
            'good': good_count,
            'bad': bad_count
        }
    return results


def update_matched_status(dataset, matches_fullscan, matches_fragfile):
    '''
    Update a boolean flag in the Chemical object that tells us if it is found in fullscan or fragmentation data
    :param dataset: a list of Chemicals
    :param matches_fullscan: the result of matching Chemicals in dataset to fullscan file
    :param matches_fragfile: the result of matching Chemicals in dataset to fragmentation file
    :return: None, but the Chemical objects in dataset is modified
    '''
    found_in_fullscan = 0
    found_in_fragfile = 0
    for chem in dataset:
        if matches_fullscan is not None:  # check if a match is found in fullscan mzML
            if chem in matches_fullscan:
                chem.found_in_fullscan = True
                found_in_fullscan += 1
            else:
                chem.found_in_fullscan = False

        if matches_fragfile is not None:  # check if a match is found in fragmentation mzML
            if chem in matches_fragfile:
                chem.found_in_fragfile = True
                found_in_fragfile += 1
            else:
                chem.found_in_fragfile = False

    print('Matched %d/%d in fullscan data, %d/%d in fragmentation data' % (found_in_fullscan, len(dataset),
                                                                           found_in_fragfile, len(dataset)))


def compute_pref_rec_f1(tp, fp, fn):
    prec = tp / (tp + fp)
    rec = tp / (tp + fn)
    f1 = (2 * prec * rec) / (prec + rec)
    return prec, rec, f1


def calculate_performance(params):
    # get parameters
    fragfile = params['fragfile']
    N = params['N']
    rt_tol = params['rt_tol']
    roi_mz_tol = params['roi_mz_tol']
    roi_min_ms1_intensity = params['roi_min_ms1_intensity']
    fragmentation_min_ms1_intensity = params['fragmentation_min_ms1_intensity']
    min_rt = params['min_rt']
    max_rt = params['max_rt']
    roi_min_length = params['roi_min_length']
    fullscan_filename = params['fullscan_filename']
    P_peaks_df = params['P_peaks_df']
    Q_peaks_df = params['Q_peaks_df']
    matching_mz_tol = params['matching_mz_tol']
    matching_rt_tol = params['matching_rt_tol']
    scenario = params['scenario']

    controller_file = params['controller_file']
    chemicals_file = params['chemicals_file']

    if chemicals_file.endswith('.p'):
        print('Loading chemicals')
        chemicals = load_obj(chemicals_file)
    else:
        print('Extracting chemicals')
        chemicals = get_chemicals(chemicals_file, roi_mz_tol, roi_min_ms1_intensity, min_rt, max_rt,
                                  min_length=roi_min_length)

    if type(chemicals) == list:
        chemicals = np.array(chemicals)

    if controller_file.endswith('.p'):
        print('Loading fragmentation events')
        controller = load_obj(controller_file)
        chem_to_frag_events = None
    else:
        print('Extracting fragmentation events')
        controller = None
        precursor_df = get_precursor_info(controller_file)
        chem_to_frag_events = get_chem_to_frag_events(chemicals, precursor_df)

    # compute performance under each scenario
    print('Computing performance under scenario %d' % scenario)
    tp, fp, fn, prec, rec, f1 = 0, 0, 0, 0, 0, 0
    if scenario == 1:
        tp, fp, fn, prec, rec, f1 = compute_performance_scenario_1(controller, chemicals,
                                                                   fragmentation_min_ms1_intensity,
                                                                   fullscan_filename, P_peaks_df,
                                                                   matching_mz_tol, matching_rt_tol,
                                                                   chem_to_frag_events=chem_to_frag_events)
    elif scenario == 2:
        fragfile_filename = os.path.basename(fragfile)
        tp, fp, fn, prec, rec, f1 = compute_performance_scenario_2(controller, chemicals,
                                                                   fragmentation_min_ms1_intensity,
                                                                   fullscan_filename, fragfile_filename,
                                                                   P_peaks_df, Q_peaks_df, matching_mz_tol,
                                                                   matching_rt_tol,
                                                                   chem_to_frag_events=chem_to_frag_events)

    return N, rt_tol, scenario, tp, fp, fn, prec, rec, f1


def evaluate_serial(all_params):
    results = []
    for params in all_params:
        res = calculate_performance(params)
        results.append(res)
        print('N=%d rt_tol=%d scenario=%d tp=%d fp=%d fn=%d prec=%.3f rec=%.3f f1=%.3f\n' % res)
    result_df = pd.DataFrame(results, columns=['N', 'rt_tol', 'scenario', 'TP', 'FP', 'FN', 'Prec', 'Rec', 'F1'])
    return result_df


def evaluate_parallel(all_params, pushed_dict=None):
    import ipyparallel as ipp
    rc = ipp.Client()
    dview = rc[:]  # use all engines​
    with dview.sync_imports():
        import os
        from VMSfunctions.TopNExperiment import get_chemicals, get_precursor_info, get_chem_to_frag_events
        from VMSfunctions.Common import load_obj

    if pushed_dict is not None:
        dview.push(pushed_dict)

    results = dview.map_sync(calculate_performance, all_params)
    result_df = pd.DataFrame(results, columns=['N', 'rt_tol', 'scenario', 'TP', 'FP', 'FN', 'Prec', 'Rec', 'F1'])
    return result_df
