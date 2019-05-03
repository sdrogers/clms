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
dataset = load_obj(os.path.join(base_dir, 'models\\dda_results\\dataset.p'))
noisy_dataset = load_obj(os.path.join(base_dir, 'models\\dda_results\\noisy_dataset.p'))

isolation_window = 0.5   # the isolation window in Dalton around a selected precursor ion
ionisation_mode = POSITIVE
data = noisy_dataset
density = ps.density_estimator
min_ms1_intensity = 2.5E5 # minimum ms1 intensity to fragment
min_rt = 3*60
max_rt = 21*60
N = 10
rt_tol = 15

# isolation_window = 0.5   # the isolation window in Dalton around a selected precursor ion
# ionisation_mode = POSITIVE
# data = dataset[0:100]
# density = ps.density_estimator
# min_ms1_intensity = 2.5E5 # minimum ms1 intensity to fragment
# min_rt = 200
# max_rt = 600
# N = 10
# rt_tol = 15

set_log_level_warning()
mass_spec = IndependentMassSpectrometer(ionisation_mode, data, density=density)
controller = TopNController(mass_spec, N, isolation_window,
                            rt_tol, min_ms1_intensity)
controller.run(min_rt, max_rt, True)
controller.write_mzML('my_analysis', '../models/test/mzML/out.mzML')