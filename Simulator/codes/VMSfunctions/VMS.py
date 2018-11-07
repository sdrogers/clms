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
        self.ms1_range = [(self.ms1_range[0][0],self.ms1_range[1][0])]
                
    def ms1_values(self):
        ms1_values = []
        for dataset_index in range(0,len(self.peaks)):
            if self.peaks[dataset_index].ms_level == 1:
                ms1_values.append(self.peaks[dataset_index].mz)
        ms1_values = sum(np.array(ms1_values).tolist(),[])
        return(ms1_values)
    
class Dataset_Scan(object):
    def __init__(self,dataset,ms_level,rt,isolation_windows):
        self.scan_result = []
        for peak in range(0,len(dataset.peaks)):
            mz_internal = dataset.peaks[peak].get(ms_level,rt,isolation_windows)
            if ms_level == dataset.peaks[peak].ms_level:
                if mz_internal == None:
                    self.scan_result.append(0)
                else:
                    self.scan_result.append(1)  
                    
                    
import math
class Dia_Basic(object):
    def __init__(self,dataset,num_windows,window_type,rt=0,ms_level=2,ms1_range=[(None,None)],range_slack=0.01):
        self.locations = []
        self.scan_results = []
        # if no ms1_range provided, initialise the scan range from dataset
        if ms1_range==[(None,None)]:
            ms1_range_difference = dataset.ms1_range[0][1] - dataset.ms1_range[0][0]
            ms1_range = [[dataset.ms1_range[0][0] - range_slack*ms1_range_difference],[dataset.ms1_range[0][1] + range_slack*ms1_range_difference]]
        else:
            ms1_range_difference = ms1_range[0][1] - ms1_range[0][0] 
        # Look at whether components are within the different isolation windows
        if window_type=="even":
            for window_index in range(0,num_windows):
                ms1_window_range = [(ms1_range[0][0] + window_index*((ms1_range_difference * (1 + 2*range_slack))/num_windows),ms1_range[0][0] + (window_index+1)*((ms1_range_difference * (1 + 2*range_slack))/num_windows))]
                self.locations.append(ms1_window_range)
                self.scan_results.append(Dataset_Scan(dataset,ms_level,rt,ms1_window_range).scan_result)
        elif window_type=="percentile":
            adjusted_percentiles = np.percentile(dataset.ms1_values(),range(0,100 + int(100/num_windows),int(100/num_windows))).tolist()
            adjusted_percentiles[0] = adjusted_percentiles[0] - range_slack*ms1_range_difference
            adjusted_percentiles[-1] = adjusted_percentiles[-1] + range_slack*ms1_range_difference
            for window_index in range(0,num_windows):
                ms1_window_range = [(adjusted_percentiles[window_index],adjusted_percentiles[window_index+1])]
                self.locations.append(ms1_window_range)
                self.scan_results.append(Dataset_Scan(dataset,ms_level,rt,ms1_window_range).scan_result)
        else:
            sys.exit("Incorrect window_type specified")
        self.ms1_names = dataset.ms1_names
        self.ms2_names = dataset.ms2_names
        
    def bin_walls(self):
        bin_walls = list(set(np.array(sum(d1.locations,[])).flatten()))
        return(bin_walls)
    
    def location_all(self):
        location_all = []
        for i in range(0,len(self.ms1_names)):
            found=False
            j=0
            while found==False:
                if self.scan_results[j][i]==1:
                    location_all.append(self.locations[j])
                    found=True
                else:
                    j+=1  
        return(location_all)
            
    def entropy(self):
        entropy_vec = []
        for i in range(0,len(self.locations)):
            entropy_vec.extend([0])
            for j in range(0,len(self.ms1_names)):
                if self.locations[i]==self.location_all()[j]:
                    entropy_vec[i] +=1
        entropy_vec_nozero = [value for value in entropy_vec if value != 0]
        entropy = sum([math.log(1/y) for y in entropy_vec_nozero])
        return(entropy)