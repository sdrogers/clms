import sys

sys.path.append('../codes')

from VMSfunctions.Controller import *
from VMSfunctions.Common import *

set_log_level_debug()
ps = load_obj('../models/peak_sampler_4_beers.p')
xcms_output = '../models/beer_ms1_peaks.csv.gz'
cc = ChromatogramCreator(xcms_output)

max_rt = 10  # the maximum retention time of scans to generate
N = 5  # top-5 DDA fragmentation
mz_tol = 5  # the mz isolation window around a selected precursor ion
rt_tol = 15  # the rt window around a selected precursor ion to prevent it from fragmented multiple times
min_ms2_intensity = 5000  # the minimum ms2 peak intensity

set_log_level_warning()
mass_spec = IndependentMassSpectrometer(POSITIVE, cc.chemicals, density=ps.density_estimator)
controller = TopNController(mass_spec, N, mz_tol, rt_tol, min_ms2_intensity=min_ms2_intensity)
controller.run(0, max_rt)
