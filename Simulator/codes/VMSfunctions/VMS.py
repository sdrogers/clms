import numpy as np
import math
from VMSfunctions.Common import load_obj
from VMSfunctions.DataGenerator import PeakSampler

class Sample_Dataset:
    def __init__ (self,file_location,n):
        ps = load_obj(file_location)
        ms1 = ps.sample(ms_level=1, n_peaks=n)
        ms2 = ps.sample(ms_level=2, n_peaks=n)
        for i in range(0,n):
            ms2[i].parent = ms1[i]
            ms1[i].rt = 0
            ms2[i].rt = 0
        self.peaks = ms1 + ms2
        self.ms1_names = []
        self.ms2_names = []
        self.ms1_range = [math.inf,-math.inf]
        for dataset_index in range(0,len(self.peaks)):
            if self.peaks[dataset_index].ms_level == 2:
                self.ms1_names.append(self.peaks[dataset_index].parent)
                self.ms2_names.append(self.peaks[dataset_index])
                self.ms1_range = np.array([min(self.ms1_range[0],self.peaks[dataset_index].parent.mz), max(self.ms1_range[1],self.peaks[dataset_index].parent.mz)]).tolist()
                
    def ms1_values(self):
        self.ms1_values = []
        for dataset_index in range(0,len(self.peaks)):
            if self.peaks[dataset_index].ms_level == 1:
                self.ms1_values.append(self.peaks[dataset_index].mz)
        self.ms1_values = sum(np.array(self.ms1_values).tolist(),[])