import sys

sys.path.append('../codes')

from VMSfunctions.Chemicals import *
from VMSfunctions.Chromatograms import *
from VMSfunctions.MassSpec import *
from VMSfunctions.Controller import *
from VMSfunctions.Common import *
from VMSfunctions.DataGenerator import *
from VMSfunctions.DsDA import *

base_dir = 'C:\\Users\\Vinny\\OneDrive - University of Glasgow\\CLDS Metabolomics Project\\'
schedule_dir = 'C:\\Users\\Vinny\\OneDrive - University of Glasgow\\CLDS Metabolomics Project\\DsDA\\DsDA_Test\\settings'
data_dir = 'C:\\Users\\Vinny\\OneDrive - University of Glasgow\\CLDS Metabolomics Project\\DsDA\\DsDA_Test\\Data'

ps = load_obj(os.path.join(base_dir, 'Trained Models\\peak_sampler_beer1pos_fragmentation.p'))
cc = ChromatogramCreator(os.path.join(base_dir, 'Trained Models\\chromatogram_beer1pos.csv.gz'))
chemicals = ChemicalCreator(ps)

min_ms1_intensity = 2E5
rt_range = [(180, 1260)]
mz_range = [(0, 1050)]
chemicals = ChemicalCreator(ps)
n_peaks = 1000
dataset = chemicals.sample(cc, mz_range, rt_range, min_ms1_intensity, n_peaks, 2, "Unknown", None, None, fixed_mz = True)

# set_log_level_warning()
# last_schedule = get_schedule(29, schedule_dir)
# mass_spec = DsDAMassSpec(POSITIVE, dataset, density=ps.density_estimator)
# controller = DsDAController(mass_spec, 1, 0.5, 15, 2E5)
# controller.run(last_schedule)
#
# controller.write_mzML('my_analysis', data_dir + '\\hello.mzML')

isolation_window = 1   # the isolation window in Dalton around a selected precursor ion
N = 4
rt_tol = 15
mz_tol = 10
set_log_level_warning()
last_schedule = get_schedule(29, schedule_dir)
mass_spec = IndependentMassSpectrometer(POSITIVE, dataset, density=ps.density_estimator, schedule_file=last_schedule)
controller = TopNController(mass_spec, N, isolation_window, mz_tol, rt_tol, min_ms1_intensity)
controller.run(rt_range[0][0], rt_range[0][1])