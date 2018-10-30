import pickle

def save_obj(obj, filename):
    with open(filename, 'wb') as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)    

def load_obj(filename):
    with open(filename, 'rb') as f:        
        return pickle.load(f)            

class Peak(object):
    def __init__(self, mz, rt, intensity, ms_level, parent=None, filename=None, scan_number=None):
        self.mz = mz
        self.rt = rt
        self.intensity = intensity
        self.ms_level = ms_level
        self.filename = filename
        self.scan_number = scan_number
        self.parent = parent
        
    def __repr__(self):
        return 'Peak mz=%.4f rt=%.2f intensity=%.2f ms_level=%d' % (self.mz, self.rt, self.intensity, self.ms_level)
    
    def get(self,ms_level,rt,isolation_window):
        if self.ms_level != ms_level:
            self.get_parent_intensity = None
            self.get_parent_mz = None
            self.get_intensity = None
            self.get_mz = None
            self.in_range = None
        elif ms_level==1: # ms1 scan of ms1 Peak
            self.get_intensity = self.intensity
            self.get_mz = self.mz
            self.get_parent_intensity = None
            self.get_parent_mz = None
            self.in_range = None
        else: # ms2 scan of ms2 Peak
            self.in_range = 0
            for i in range(0,len(isolation_window[0])):
                if self.parent.mz > isolation_window[0][i] and self.parent.mz <= isolation_window[1][i]:
                    self.in_range=1
                    self.get_parent_intensity = self.parent.intensity
                    self.get_parent_mz = self.parent.mz
                    self.get_intensity = self.intensity
                    self.get_mz = self.mz
            if self.in_range == 0:
                self.get_parent_intensity = None
                self.get_parent_mz = None
                self.get_intensity = None
                self.get_mz = None
            