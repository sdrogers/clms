import seaborn as sns

from VMSfunctions.Controller import *


def get_N(row):
    if 'T10' in row['filename']:
        return 10
    else:
        return row['filename'].split('_')[3]


def get_dew(row):
    if 'T10' in row['filename']:
        return 15
    else:
        tok = row['filename'].split('_')[5] # get the dew value in the filename
        return tok.split('.')[0] # get the part before '.mzML'


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
    intensity_col = 'maxo'
    df = df[(df['rt'] > rt_range[0][0]) & (df['rt'] < rt_range[0][1])]
    df = df[(df['rt'] > mz_range[0][0]) & (df['rt'] < mz_range[0][1])]
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


def make_boxplot(df, x, y, xticklabels, title):
    g = sns.catplot(x=x, y=y, kind='box', data=df)
    g.fig.set_size_inches(10, 3)
    if xticklabels is not None:
        g.set_xticklabels(xticklabels, rotation=90)
    else:
        g.set_xticklabels(rotation=90)
    plt.title(title)
    plt.tight_layout()
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


def compute_performance(controller, dataset, min_ms1_intensity):
    ms_level = 2
    chem_to_frag_events = get_frag_events(controller, ms_level)

    # positive instances are ground truth MS1 peaks found by XCMS
    # negative instances are chemicals that cannot be matched to XCMS output
    positives = list(filter(lambda x: x.found_in_ms1, dataset))
    negatives = list(filter(lambda x: not x.found_in_ms1, dataset))

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
    prec = tp / (tp + fp)
    rec = tp / (tp + fn)
    f1 = (2 * prec * rec) / (prec + rec)
    prec, rec, f1
    return tp, fp, fn, prec, rec, f1


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


def get_absolute_intensity(chem, query_rt):
    return chem.max_intensity * chem.chromatogram.get_relative_intensity(query_rt - chem.rt)


def count_frag_events(key, chem_to_frag_events, min_ms1_intensity):
    # count how many good and bad fragmentation events for each chemical (key)
    frag_events = chem_to_frag_events[key]
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
        key = get_key(chem)
        try:
            good_count, bad_count = count_frag_events(key, chem_to_frag_events, min_ms1_intensity)
        except KeyError:
            good_count = 0
            bad_count = 0
        results[chem] = {
            'good': good_count,
            'bad': bad_count
        }
    return results
