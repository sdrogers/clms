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
    
    
    def get(self,ms_level,rt,isolation_windows):
        if not ms_level == self.ms_level:
            return None
        if not self._rt_match(rt):
            return None
        if ms_level == 1:
            if self._isolation_match(isolation_windows):
                return (self._get_mz(rt),self._get_intensity())
            else:
                return None
        # if we get here, it's ms_level >1 and rt is ok
        if(self.parent._isolation_match(isolation_windows)):
            return (self._get_mz(rt),self._get_intensity(rt))

    def _get_mz(self,rt):
        return self.mz

    def _get_intensity(self,rt):
        return self.intensity

    def _rt_match(self,rt):
        if not rt:
            return True
        if self.ms_level == 1:
            return True # eventually, some checking here and possibly returning false
        else:
            return self.parent._rt_match(rt)
    def _isolation_match(self,isolation_windows):
        # assumes list is formated like:
        # [(min_1,max_1),(min_2,max_2),...],
        if not isolation_windows:
            return True
        for window in isolation_windows:
            if(self.mz > window[0] and self.mz < window[1]):
                return True
        return False
        # for i in range(0,len(isolation_window[0])):
        #         if self.parent.mz > isolation_window[0][i] and self.parent.mz <= isolation_window[1][i]: