import os
import pathlib

from VMSfunctions.Controller import *


def run_experiment(param):
    analysis_name = param['analysis_name']
    mzml_out = param['mzml_out']
    pickle_out = param['pickle_out']
    if os.path.isfile(mzml_out) and os.path.isfile(pickle_out):
        print('Skipping %s' % (analysis_name))
    else:
        print('Processing %s' % (analysis_name))
        mass_spec = IndependentMassSpectrometer(param['ionisation_mode'], param['data'], density=param['density'])
        controller = TopNController(mass_spec, param['N'], param['isolation_window'],
                                    param['mz_tol'], param['rt_tol'], param['min_ms1_intensity'])
        controller.run(param['min_rt'], param['max_rt'], progress_bar=param['pbar'])
        controller.write_mzML(analysis_name, mzml_out)
        save_obj(controller, pickle_out)
        return analysis_name


def run_parallel_experiment(params):
    import ipyparallel as ipp
    rc = ipp.Client()
    dview = rc[:]  # use all engines​
    with dview.sync_imports():
        pass

    analysis_names = dview.map_sync(run_experiment, params)
    for analysis_name in analysis_names:
        print(analysis_name)


def run_serial_experiment(param, i, total):
    print('Processing \t%d/%d\t%s' % (i + 1, total, param['analysis_name']))
    run_experiment(param)


def get_params(experiment_name, Ns, rt_tols, mz_tol, isolation_window, ionisation_mode, data, density,
               min_ms1_intensity, min_rt, max_rt,
               out_dir, pbar):
    create_if_not_exist(out_dir)

    print('N =', Ns)
    print('rt_tol =', rt_tols)
    params = []
    for N in Ns:
        for rt_tol in rt_tols:
            analysis_name = 'experiment_%s_N_%d_rttol_%d' % (experiment_name, N, rt_tol)
            mzml_out = os.path.join(out_dir, '%s.mzML' % analysis_name)
            pickle_out = os.path.join(out_dir, '%s.p' % analysis_name)
            params.append({
                'N': N,
                'mz_tol': mz_tol,
                'rt_tol': rt_tol,
                'min_ms1_intensity': min_ms1_intensity,
                'isolation_window': isolation_window,
                'ionisation_mode': ionisation_mode,
                'data': data,
                'density': density,
                'min_rt': min_rt,
                'max_rt': max_rt,
                'analysis_name': analysis_name,
                'mzml_out': mzml_out,
                'pickle_out': pickle_out,
                'pbar': pbar
            })
    print('len(params) =', len(params))
    return params

def get_key(chem):
    # turn a chem into (mz, rt, intensity) for equal comparison
    return (tuple(chem.isotopes), chem.rt, chem.max_intensity)

def get_frag_events(controller, ms_level):
    # get the fragmentation events for all chemicals for an ms level
    filtered_frag_events = list(filter(lambda x: x.ms_level == ms_level, controller.mass_spec.fragmentation_events))
    chem_to_frag_events = defaultdict(list)
    for frag_event in filtered_frag_events:
        key = get_key(frag_event.chem)
        chem_to_frag_events[key].append(frag_event)
    return dict(chem_to_frag_events)

def count_frag_events(key, chem_to_frag_events):
    # count how many good and bad fragmentation events for each chemical (key)
    frag_events = chem_to_frag_events[key]
    good_count = 0
    bad_count = 0
    for frag_event in frag_events:
        chem = frag_event.chem
        rt_match = chem.chromatogram._rt_match(frag_event.query_rt - chem.rt)
        if rt_match:
            good_count += 1
        else:
            bad_count += 1
    return good_count, bad_count

def get_chem_frag_counts(chem_list, chem_to_frag_events):
    # get the count of good/bad fragmentation events for all chemicals in chem_list
    results = {}
    for i in range(len(chem_list)):
        chem = chem_list[i]
        key = get_key(chem)
        try:
            good_count, bad_count = count_frag_events(key, chem_to_frag_events)
        except KeyError:
            good_count = 0
            bad_count = 0
        results[chem] = {
            'good': good_count,
            'bad': bad_count
        }
    return results

def compute_performance(controller, dataset):
    ms_level = 2
    chem_to_frag_events = get_frag_events(controller, ms_level)
    positives = list(filter(lambda x: x.type == 'data', dataset))
    negatives = list(filter(lambda x: x.type == 'noise', dataset))
    positives_count = get_chem_frag_counts(positives, chem_to_frag_events)
    negatives_count = get_chem_frag_counts(negatives, chem_to_frag_events)

    # count the following:
    # true positive = is an xcms peak (positive) and is fragmented within the chemical's elution time
    # false positive = is not an xcms peak (negative) and is fragmented within the chemical's elution time
    # false negative = is an xcms peak (positive) and is not fragmented within the chemical's elution time
    tp = len([chem for chem in positives if positives_count[chem]['good'] > 0])
    fp = len([chem for chem in negatives if negatives_count[chem]['good'] > 0])
    fn = len([chem for chem in positives if positives_count[chem]['good'] == 0])

    prec = tp / (tp + fp)
    rec = tp / (tp + fn)
    f1 = ( 2 * prec * rec) / (prec + rec)
    prec, rec, f1
    return tp, fp, fn, prec, rec, f1

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

def make_plot(df, X, Y, title, ylabel):
    df.plot.line(x=X, y=Y)
    plt.title(title)
    plt.ylabel(ylabel)