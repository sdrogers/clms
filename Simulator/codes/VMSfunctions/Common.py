class Peak(object):
    def __init__(self, mz, rt, intensity, ms_level, parent, filename=None, scan_number=None):
        self.mz = mz
        self.rt = rt
        self.intensity = intensity
        self.ms_level = ms_level
        self.filename = filename
        self.scan_number = scan_number
        self.parent = parent