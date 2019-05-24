import sys

sys.path.append('../codes')

from VMSfunctions.Controller import *
from VMSfunctions.DataGenerator import *

set_log_level_info()

base_dir = 'C:\\Users\\joewa\\University of Glasgow\\Vinny Davies - CLDS Metabolomics Project\\'
mzml_path = os.path.join(base_dir, 'Data\\Fusion_1578_Ronan_Daly_CLP_pHILIC_22May19\\Positive\\QCB_MS2')
file_name = 'QCB_N10_DEW015.mzML'

experiment_name = 'qcb_compare_debug'
experiment_out_dir = os.path.join(base_dir, 'Simulator Outputs', experiment_name)

min_rt = 0
max_rt = 1600

kde_min_ms1_intensity = 0 # min intensity to be selected for kdes
kde_min_ms2_intensity = 0

roi_mz_tol = 10
roi_min_length = 2
roi_min_intensity = 0
roi_start_rt = min_rt
roi_stop_rt = max_rt

isolation_window = 1   # the isolation window in Dalton around a selected precursor ion
ionisation_mode = POSITIVE
N = 10
rt_tol = 15
mz_tol = 10
min_ms1_intensity = 0 # minimum ms1 intensity to fragment

ds = DataSource()
ds.load_data(mzml_path, file_name=file_name)

densities = PeakDensityEstimator(kde_min_ms1_intensity, kde_min_ms2_intensity, min_rt, max_rt, plot=True)
densities.kde(ds, file_name, 2, bandwidth_mz_intensity=1.0, bandwidth_rt=5.0, bandwidth_n_peaks=1.0)
ps = PeakSampler(densities)

mzml_file = os.path.join(mzml_path, file_name)
good_roi, junk = make_roi(mzml_file, mz_tol=roi_mz_tol, mz_units='ppm', min_length=roi_min_length,
                          min_intensity=roi_min_intensity, start_rt=roi_start_rt, stop_rt=roi_stop_rt)
all_roi = good_roi

set_log_level_debug()
rtcc = RoiToChemicalCreator(ps, all_roi)

data = rtcc.chemicals
save_obj(data, os.path.join(experiment_out_dir, 'dataset.p'))

density = ps.density_estimator
set_log_level_warning()
pbar = True

mass_spec = IndependentMassSpectrometer(ionisation_mode, data, density=density)
controller = TopNController(mass_spec, N, isolation_window, mz_tol,
                            rt_tol, min_ms1_intensity)
controller.run(min_rt, max_rt, pbar)

mzml_out = os.path.join(experiment_out_dir, 'QCB_N10_DEW015_simulated.mzML')
controller.write_mzML('my_analysis', mzml_out)