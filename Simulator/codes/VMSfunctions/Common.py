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