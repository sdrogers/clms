import glob
import os

import numpy as np
import pylab as plt
import pymzml
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV

from .Common import Peak

class DataSource(object):
    def __init__(self, min_ms1_intensity=0, min_ms2_intensity=0):
        self.data = {}
        self.peaks = {
            1: [],
            2: []
        }
        self.ms1_precision = 5e-6
        self.obo_version = '4.0.1'
        self.min_ms1_intensity = min_ms1_intensity
        self.min_ms2_intensity = min_ms2_intensity
        
    def load_data(self, data_path):
        for filename in glob.glob(os.path.join(data_path, '*.mzML')):
            run = pymzml.run.Reader(filename, obo_version = self.obo_version,
                                    MS1_Precision=self.ms1_precision,
                                    extraAccessions=[('MS:1000016', ['value', 'unitName'])])
            ms1 = []
            ms2 = []
            for scan_number, spectrum in enumerate(run):
                # if isinstance(spectrum, pymzml.spec.Chromatogram):
                #    print('chromatogram')
                for mz, intensity in spectrum.peaks('centroided'):
                    rt, units = spectrum.scan_time
                    if units == 'minute':
                        rt *= 60.0
                    ms_level = spectrum.ms_level
                    p = Peak(mz, rt, intensity, ms_level)
                    if ms_level == 1 and intensity > self.min_ms1_intensity:
                        ms1.append(p)
                    elif ms_level == 2 and intensity > self.min_ms2_intensity:
                        ms2.append(p)
            self.peaks[1].extend(ms1)
            self.peaks[2].extend(ms2)
            print('%s (ms1=%d, ms2=%d)' % (filename, len(ms1), len(ms2)))
                            
    def plot_histogram(self, data_type, ms_level, log=False, bins=100):
        X = self.get_data(data_type, ms_level, log)
        plt.figure()
        _ = plt.hist(X, bins=100)
        plt.plot(X[:, 0], np.full(X.shape[0], -0.01), '|k')
        plt.title('%s histogram' % data_type)
        plt.show()
                
    def get_data(self, data_type, ms_level, log=False): # data_type is 'mz', 'rt' or 'intensity'
        peaks = self.peaks[ms_level]
        values = list(getattr(x, data_type) for x in peaks)
        X = np.array(values)[:, np.newaxis]     
        if log:
            X = np.log(X)
        return X
        
class DensityEstimator(object):
    def __init__(self):
        self.bandwidths = 10.0**np.arange(-4, 1)
        self.kdes = {}
        
    def cv(self, data_source, data_type, ms_level, log=False, bandwidths=None, cv=2, plot=False):
        X = data_source.get_data(data_type, ms_level, log=log)
        bw = self.bandwidths if bandwidths is None else bandwidths
        grid = GridSearchCV(KernelDensity(kernel='gaussian', rtol=1e-2), 
                            {'bandwidth': bw},
                            cv=cv, verbose=True, n_jobs=-1)
        grid.fit(X)
        
        bandwidth_cv = grid.best_params_['bandwidth']
        print('Best bandwidth: %.3f' % bandwidth_cv)        
        title = '%s density estimation - bandwidth %.3f' % (data_type, bandwidth_cv)
        kde = self._train(bandwidth_cv, X)
        if plot:
            self._plot(kde, X, title)
        self.kdes[(data_type, ms_level)] = kde
        
    def kde(self, data_source, data_type, ms_level, log=False, bandwidth=1.0, plot=False):
        X = data_source.get_data(data_type, ms_level, log=log)
        title = '%s density estimation - bandwidth %.3f' % (data_type, bandwidth)
        kde = self._train(bandwidth, X)
        if plot:
            self._plot(kde, X, title)
        self.kdes[(data_type, ms_level)] = kde        
        
    def sample(self, data_type, ms_level, n_sample):
        kde = self.kdes[(data_type, ms_level)]
        return kde.sample(n_sample)
                
    def _train(self, bandwidth_cv, X):
        kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth_cv).fit(X)
        return kde
    
    def _plot(self, kde, X, title):
        X_plot = np.linspace(np.min(X), np.max(X), 1000)[:, np.newaxis]
        log_dens = kde.score_samples(X_plot)                
        plt.figure()        
        plt.fill_between(X_plot[:, 0], np.exp(log_dens), alpha=0.5)
        plt.plot(X[:, 0], np.full(X.shape[0], -0.01), '|k')
        plt.title(title)
        plt.show()
        
class PeakSampler(object):
    def __init__(self, densities):
        self.densities = densities
        
    def sample(self, ms_level, n_peaks):
        peaks = []
        mzs = self.densities.sample('mz', ms_level, n_peaks)
        rts = self.densities.sample('rt', ms_level, n_peaks)
        intensities = np.exp(self.densities.sample('intensity', ms_level, n_peaks))
        for i in range(n_peaks):
            p = Peak(mzs[i], rts[i], intensities[i], ms_level)
            peaks.append(p)
        return peaks