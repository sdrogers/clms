import numpy as np
import math
from VMSfunctions.Common import load_obj
from VMSfunctions.DataGenerator import PeakSampler

class Sample_Dataset:
    def __init__ (self,file_location,n_ms1=None,n_ms2=None):
        ps = load_obj(file_location)
        self.peaks = []
        if n_ms1==None:
            n_ms1 = len(ps.sample(ms_level=1))
        for i in range(0,n_ms1):
            ms1 = ps.sample(ms_level=1,n_peaks=1)
            ms1[0].rt=0
            if n_ms2 == None:
                ms2 = ps.sample(ms_level=2)
                if len(ms2)<0:
                    ms2 = ps.sample(ms_level=2,n_peaks=0)
            else:
                ms2 = ps.sample(ms_level=2,n_peaks=n_ms2)
            for j in range(0,len(ms2)):
                ms2[j].parent = ms1[0]
                ms2[j].rt = 0
            self.peaks.extend(ms1)
            self.peaks.extend(ms2)
        self.ms1_names = []
        self.ms2_names = []
        self.ms1_range = [math.inf,-math.inf]
        for dataset_index in range(0,len(self.peaks)):
            if self.peaks[dataset_index].ms_level == 2:
                self.ms2_names.append(self.peaks[dataset_index])
            elif self.peaks[dataset_index].ms_level == 1:
                self.ms1_names.append(self.peaks[dataset_index])
                self.ms1_range = np.array([min(self.ms1_range[0],self.peaks[dataset_index].mz), max(self.ms1_range[1],self.peaks[dataset_index].mz)]).tolist()
            else:
                sys.exit("MS level beyond MS1 / MS2. This has not been implemented yet")
        self.ms1_range = [(self.ms1_range[0],self.ms1_range[1])]
    def ms1_values(self):
        ms1_values = []
        for dataset_index in range(0,len(self.peaks)):
            if self.peaks[dataset_index].ms_level == 1:
                ms1_values.append(self.peaks[dataset_index].mz)
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
        ms2_vec = []
        locations_all=self.location_all()
        for i in range(0,len(self.locations)):
            ms2_vec.extend([0])
            for j in range(0,len(self.ms1_names)):
                if self.locations[i]==locations_all[j]:
                    ms2_vec[i] +=1
        ms2_vec_nozero = [value for value in ms2_vec if value != 0]
        entropy = sum([-y*math.log(1/y) for y in ms2_vec_nozero])
        return(entropy)
    
import math
class Dia_Kaufmann(object):
    def __init__(self,dataset,num_windows,window_type,kaufmann_design,rt=0,ms_level=2,ms1_range=[(None,None)], range_slack=0.01,extra_bins=0):
        if num_windows!=64:
            print("num_windows must equal 64. Has not been set up for any other values yet.")
            return None
        self.num_windows=num_windows
        self.extra_bins=extra_bins
        self.ms1_values = dataset.ms1_values()
        self.locations = []
        self.scan_results = []
        # if no ms1_range provided, initialise the scan range from dataset
        if ms1_range==[(None,None)]:
            ms1_range_difference = dataset.ms1_range[0][1] - dataset.ms1_range[0][0]
            ms1_range=dataset.ms1_range
        else:
            ms1_range_difference = ms1_range[0][1] - ms1_range[0][0] 
        # get names
        self.ms1_names = dataset.ms1_names
        self.ms2_names = dataset.ms2_names
        # Create bin_walls        
        if window_type=="even":
            self.bin_walls=[min(dataset.ms1_values())]
            for window_index in range(0,num_windows):
                self.bin_walls.append(ms1_range[0][0] + ((window_index+1)/num_windows)*ms1_range_difference)
            self.bin_walls[0] = self.bin_walls[0] - range_slack*ms1_range_difference
            self.bin_walls[-1] = self.bin_walls[-1] + range_slack*ms1_range_difference
            self.bin_walls_extra = None
            if extra_bins > 0:
                self.bin_walls_extra=[min(dataset.ms1_values())]
                for window_index in range(0,num_windows*(2**extra_bins)):
                    self.bin_walls_extra.append(ms1_range[0][0] + ((window_index+1)/(num_windows*(2**extra_bins)))*ms1_range_difference)
                self.bin_walls_extra[0] = self.bin_walls_extra[0] - range_slack*ms1_range_difference
                self.bin_walls_extra[-1] = self.bin_walls_extra[-1] + range_slack*ms1_range_difference
        elif window_type=="percentile":
            self.bin_walls = np.percentile(dataset.ms1_values(),np.arange(0,100 + 100/num_windows,100/num_windows)).tolist()
            self.bin_walls[0] = self.bin_walls[0] - range_slack*ms1_range_difference
            self.bin_walls[-1] = self.bin_walls[-1] + range_slack*ms1_range_difference
            self.bin_walls_extra = None
            if extra_bins > 0:
                self.bin_walls_extra = np.percentile(dataset.ms1_values(),np.arange(0,100 + 100/(num_windows*(2**extra_bins)), 100/(num_windows*(2**extra_bins)))).tolist()
                self.bin_walls_extra[0] = self.bin_walls_extra[0] - range_slack*ms1_range_difference
                self.bin_walls_extra[-1] = self.bin_walls_extra[-1] + range_slack*ms1_range_difference
        else:
            sys.exit("Incorrect window_type specified")

        # use separate class which creates windows here
        self.locations=Kaufmann_Windows(self.bin_walls,self.bin_walls_extra,kaufmann_design,extra_bins).locations
        # scan locations
        self.scan_results = []
        for window_index in range(0,len(self.locations)):
            self.scan_results.append(Dataset_Scan(dataset,ms_level,rt,self.locations[window_index]).scan_result)

    def location_finder(self,num_scans):
        return(Dia_Location_Finder(self.scan_results[0:num_scans],self.locations[0:num_scans],"kaufmann",self.extra_bins))
    
    def entropy(self,start):
        entropy = [None]*start
        components_determined = [None]*start
        for i in range(start,len(self.locations)+1):
            ms1_vec = []
            ms2_vec = []
            locations_finder=self.location_finder(i)
            locations_all=locations_finder.location_all
            bin_walls=locations_finder.bin_walls
            for i in range(0,len(bin_walls)-1):
                ms2_vec.extend([0])
                ms1_vec.extend([len(np.where(np.logical_and(np.array(self.ms1_values)>bin_walls[i], np.array(self.ms1_values)<=bin_walls[i+1]))[0])])
                for j in range(0,len(self.ms2_names)):
                    if [(bin_walls[i],bin_walls[i+1])]==locations_all[j]:
                        ms2_vec[i] +=1    
            ms1_vec_nozero = [value for value in ms1_vec if value != 0]
            ms2_vec_nozero = [value for value in ms1_vec if value != 0]
            entropy_vec = []
            for j in range(0,len(ms2_vec_nozero)):
                entropy_vec.append(-ms2_vec_nozero[j] *ms1_vec_nozero[j]*math.log(1/ms1_vec_nozero[j])) 
            entropy.append(sum(entropy_vec))
            components_determined.append(sum(np.extract(np.array(ms1_vec_nozero)==1,ms2_vec_nozero)))
            components = sum(ms2_vec_nozero)
        return(entropy,components_determined,components)

class Kaufmann_Windows(object):
    def __init__(self,bin_walls,bin_walls_extra,kaufmann_design,extra_bins):
        if kaufmann_design == "nested":
            self.locations = []
            for i in range(0,8):
                self.locations.append([(bin_walls[(0+i*8)],bin_walls[(8+i*8)])])
            if extra_bins > 0:
                locations_internal =  [ [] for i in range(4 + extra_bins) ]
            else:
                locations_internal = [[],[],[],[]]
            for i in range(0,4):
                locations_internal[0].append((bin_walls[(4 + i*16)],bin_walls[(12 + i*16)]))
                locations_internal[1].append((bin_walls[(2 + i*16)],bin_walls[(6 + i*16)]))
                locations_internal[1].append((bin_walls[(10 + i*16)],bin_walls[(14 + i*16)]))
                locations_internal[2].append((bin_walls[(1 + i*16)],bin_walls[(3 + i*16)]))
                locations_internal[2].append((bin_walls[(9 + i*16)],bin_walls[(11 + i*16)]))
                locations_internal[3].append((bin_walls[(5 + i*16)],bin_walls[(7 + i*16)]))
                locations_internal[3].append((bin_walls[(13 + i*16)],bin_walls[(15 + i*16)]))
            if extra_bins > 0:
                for j in range(extra_bins):
                    for i in range(64*(2**j)):
                        locations_internal[4+j].append((bin_walls_extra[int(0  + i*((2**extra_bins)/(2**j)))], bin_walls_extra[int(((2**extra_bins)/(2**j))/2  + i*((2**extra_bins)/(2**j)))]))
            self.locations.extend(locations_internal)    
        elif kaufmann_design == "tree":
            self.locations = []
            self.locations.append([(bin_walls[0],bin_walls[32])])
            self.locations.append([(bin_walls[32],bin_walls[64])])
            self.locations.append([(bin_walls[16],bin_walls[48])])
            self.locations.append([(bin_walls[8],bin_walls[24]),(bin_walls[40],bin_walls[56])])
            if extra_bins > 0:
                locations_internal =  [ [] for i in range(3 + extra_bins) ]
            else:
                locations_internal = [[],[],[]]
            for i in range(0,4):
                locations_internal[0].append((bin_walls[(4 + i*16)],bin_walls[(12 + i*16)]))
                locations_internal[1].append((bin_walls[(2 + i*16)],bin_walls[(6 + i*16)]))
                locations_internal[1].append((bin_walls[(10 + i*16)],bin_walls[(14 + i*16)]))
                locations_internal[2].append((bin_walls[(1 + i*16)],bin_walls[(3 + i*16)]))
                locations_internal[2].append((bin_walls[(9 + i*16)],bin_walls[(11 + i*16)]))
                locations_internal[2].append((bin_walls[(5 + i*16)],bin_walls[(7 + i*16)]))
                locations_internal[2].append((bin_walls[(13 + i*16)],bin_walls[(15 + i*16)]))
            if extra_bins > 0:
                for j in range(extra_bins):
                    for i in range(64*(2**j)):
                        locations_internal[3+j].append((bin_walls_extra[int(0  + i*((2**extra_bins)/(2**j)))], bin_walls_extra[int(((2**extra_bins)/(2**j))/2  + i*((2**extra_bins)/(2**j)))]))
            self.locations.extend(locations_internal)
        else:
            print("not a valid design")
            
class Dia_Location_Finder(object):
    def __init__(self,scan_results,locations,dia_method,extra_bins):
        bin_walls = []
        for i in range(0,len(locations)):
            bin_walls.extend(list(set(np.array(locations[i]).flatten().tolist())))
        bin_walls=list(set(bin_walls))
        bin_walls.sort()
        self.bin_walls=bin_walls   
        self.scan_results=scan_results
        self.locations=locations
        if dia_method=="kaufmann":
            self.location_all = []
            bin_mid_points = list((np.array(bin_walls[1:]) + np.array(bin_walls[:(len(bin_walls)-1)])) / 2)
            for sample_index in range(0,len(scan_results[0])):
                mid_point_TF = []
                for mid_points_index in range(0,len(bin_mid_points)):
                    mid_point_TF.append(self._mid_point_in_location(bin_mid_points[mid_points_index],sample_index))
                self.location_all.append([(list(np.array(bin_walls)[np.where(np.array(mid_point_TF)==True)])[0],list(np.array(bin_walls[1:])[np.where(np.array(mid_point_TF)==True)])[0])])
                          
    def _mid_point_in_location(self,mid_point,sample_index):
        for locations_index in range(0,len(self.locations)):
            if self._in_window(mid_point,self.locations[locations_index])==True and self.scan_results[locations_index][sample_index]==0:   
                #print(sample_index)
                #print('hello1')
                #print(self.scan_results[locations_index][sample_index])
                #print(self._in_window(mid_point,self.locations[locations_index]))
                return False
            if self._in_window(mid_point,self.locations[locations_index])==False and self.scan_results[locations_index][sample_index]==1:
                return False
        else:
            return True
            
    def _in_window(self,mid_point,locations):
        for window in locations:
            #print(mid_point)
            #print(window[0])
            #print(window[1])
            if(mid_point > window[0] and mid_point <= window[1]):
                return True
        return False
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    