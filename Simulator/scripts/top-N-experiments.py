import sys

sys.path.append('../codes')

from VMSfunctions.Chemicals import *
from VMSfunctions.Chromatograms import *
from VMSfunctions.MassSpec import *
from VMSfunctions.Controller import *
from VMSfunctions.Common import *
from VMSfunctions.DataGenerator import *

set_log_level_debug()

base_dir = 'C:\\Users\\joewa\\University of Glasgow\\Vinny Davies - CLDS Metabolomics Project\\'
ps = load_obj(base_dir + 'Trained Models\\peak_sampler_beer1pos_fragmentation.p')

base_dir = 'C:\\Users\\joewa\\Work\\git\\clms\\Simulator'
dataset = load_obj(os.path.join(base_dir, 'experiments\\beer1pos_mzml_compare\\dataset.p'))

isolation_window = 1  # the isolation window in Dalton around a selected precursor ion
ionisation_mode = POSITIVE
data = dataset
density = ps.density_estimator
min_ms1_intensity = 1.75E5 # minimum ms1 intensity to fragment
min_rt = 3*60
max_rt = 21*60
N = 10
mz_tol = 10
rt_tol = 15
pbar = False

set_log_level_warning()
set_log_level_debug()
mass_spec = IndependentMassSpectrometer(ionisation_mode, data, density=density)
controller = TopNController(mass_spec, N, isolation_window, mz_tol,
                            rt_tol, min_ms1_intensity)
controller.run(min_rt, max_rt, pbar)
# controller.write_mzML('my_analysis', '../experiments/beer1pos_mzml_compare/out.mzML')