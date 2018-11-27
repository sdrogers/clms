import glob
import os

import numpy as np
import pylab as plt
import pymzml
from sklearn.neighbors import KernelDensity

from .Common import Peak, MZ, RT, INTENSITY, N_PEAKS, MZ_INTENSITY


class DataSource(object):
    """A class to load and extract centroided peaks from mzML files.
    """
    def __init__(self, min_ms1_intensity=0, min_ms2_intensity=0):
        self.data = {}
        self.peaks = { # all peaks across all scans
            1: [],
            2: []
        }
        self.num_peaks = { # number of peaks for each ms_level per scan
            1: [],
            2: []
        }
        self.ms1_precision = 5e-6
        self.obo_version = '4.0.1'
        self.min_ms1_intensity = min_ms1_intensity
        self.min_ms2_intensity = min_ms2_intensity

    def load_data(self, data_path):
        for filename in glob.glob(os.path.join(data_path, '*.mzML')):
            run = pymzml.run.Reader(filename, obo_version=self.obo_version,
                                    MS1_Precision=self.ms1_precision,
                                    extraAccessions=[('MS:1000016', ['value', 'unitName'])])

            total_peaks = { 1: 0, 2: 0 }
            for scan_number, spectrum in enumerate(run):
                spectrum_peaks = { 1: [], 2: [] }

                # loop through spectrum and get all peaks above threshold
                for mz, intensity in spectrum.peaks('centroided'):
                    ms_level = spectrum.ms_level
                    rt, units = spectrum.scan_time
                    if units == 'minute':
                        rt *= 60.0

                    p = Peak(mz, rt, intensity, ms_level)
                    min_intensity = self.min_ms1_intensity if ms_level == 1 else self.min_ms2_intensity
                    if intensity > min_intensity:
                        spectrum_peaks[ms_level].append(p)
                        # self.peaks[ms_level].append(p)
                        # count[ms_level] += 1

                # store the results from each spectrum
                for ms_level in spectrum_peaks:
                    n = len(spectrum_peaks[ms_level])
                    if n > 0:
                        self.peaks[ms_level].extend(spectrum_peaks[ms_level])
                        self.num_peaks[ms_level].append(n)
                        total_peaks[ms_level] += n
            print('%s (ms1=%d, ms2=%d)' % (filename, total_peaks[1], total_peaks[2]))
                            
    def plot_histogram(self, data_type, ms_level, log=False, bins=100):
        X = self.get_data(data_type, ms_level, log)
        plt.figure()
        _ = plt.hist(X, bins=bins)
        plt.plot(X[:, 0], np.full(X.shape[0], -0.01), '|k')
        plt.title('%s histogram' % data_type)
        plt.show()
                
    def get_data(self, data_type, ms_level, log=False): # data_type is 'mz', 'rt', 'intensity' or 'n_peaks'
        if data_type == N_PEAKS: # get observed peak counts across all spectra in the data
            values = self.num_peaks[ms_level]
        else: # get observed mz, rt or intensity values of all peaks in the data
            peaks = self.peaks[ms_level]
            values = list(getattr(x, data_type) for x in peaks)

        # convert into Nx1 array and also log-transform if necessary
        X = np.array(values)[:, np.newaxis]
        if log:
            X = np.log(X)
        return X


class DensityEstimator(object):
    """A class to perform kernel density estimation. Takes as input a DataSource class."""
    def __init__(self):
        self.kdes = {}
        self.kernel = 'gaussian'

    def kde(self, data_source, data_type, ms_level, log=False, bandwidth=1.0, plot=False):
        X = data_source.get_data(data_type, ms_level, log=log)
        title = '%s density estimation - bandwidth %.3f' % (data_type, bandwidth)
        kde = KernelDensity(kernel=self.kernel, bandwidth=bandwidth).fit(X)
        if plot:
            self._plot(kde, X, title)
        self.kdes[(data_type, ms_level)] = kde

    def sample(self, ms_level, n_sample):
        mz_vals = self.kdes[(MZ, ms_level)].sample(n_sample)
        intensity_vals = self.kdes[(INTENSITY, ms_level)].sample(n_sample)
        rt_vals = self.kdes[(RT, ms_level)].sample(n_sample)
        vals = np.concatenate((mz_vals, intensity_vals, rt_vals), axis=1)
        return vals

    def n_peaks(self, ms_level, n_sample):
        return self.kdes[(N_PEAKS, ms_level)].sample(n_sample)

    def _plot(self, kde, X, title):
        X_plot = np.linspace(np.min(X), np.max(X), 1000)[:, np.newaxis]
        log_dens = kde.score_samples(X_plot)                
        plt.figure()        
        plt.fill_between(X_plot[:, 0], np.exp(log_dens), alpha=0.5)
        plt.plot(X[:, 0], np.full(X.shape[0], -0.01), '|k')
        plt.title(title)
        plt.show()


class PeakDensityEstimator(object):
    """A class to perform better kernel density estimation for peak data. Takes as input a DataSource class."""

    def __init__(self):
        self.kdes = {}
        self.kernel = 'gaussian'

    def kde(self, data_source, ms_level, bandwidth_mz_int=1.0, bandwidth_rt=5.0, bandwidth_n_peaks=1.0):
        # train kde on mz-intensity values
        mz = data_source.get_data(MZ, ms_level)
        intensity = data_source.get_data(INTENSITY, ms_level, log=True)
        X = np.concatenate((mz, intensity), axis=1)
        self.kdes[(MZ_INTENSITY, ms_level)] = KernelDensity(kernel=self.kernel, bandwidth=bandwidth_mz_int).fit(X)

        # train kde on rt values
        X = data_source.get_data(RT, ms_level)
        self.kdes[(RT, ms_level)] = KernelDensity(kernel=self.kernel, bandwidth=bandwidth_rt).fit(X)

        # train kde on number of peaks
        X = data_source.get_data(N_PEAKS, ms_level)
        self.kdes[(N_PEAKS, ms_level)] = KernelDensity(kernel=self.kernel, bandwidth=bandwidth_n_peaks).fit(X)

    def sample(self, ms_level, n_sample):
        mz_intensity_vals = self.kdes[(MZ_INTENSITY, ms_level)].sample(n_sample)
        rt_vals = self.kdes[(RT, ms_level)].sample(n_sample)
        vals = np.concatenate((mz_intensity_vals, rt_vals), axis=1)
        return vals

    def n_peaks(self, ms_level, n_sample):
        return self.kdes[(N_PEAKS, ms_level)].sample(n_sample)


class PeakSampler(object):
    """A class to sample peaks from a trained density estimator"""
    def __init__(self, density_estimator):
        self.density_estimator = density_estimator
        
    def sample(self, ms_level, n_peaks=None):
        if n_peaks is None:
            n_peaks = self.density_estimator.n_peaks(ms_level, 1).astype(int)[0][0]
        vals = self.density_estimator.sample(ms_level, n_peaks)
        mzs = vals[:, 0]
        intensities = np.exp(vals[:, 1])
        rts = vals[:, 2]
        peaks = []
        for i in range(n_peaks):
            p = Peak(mzs[i], rts[i], intensities[i], ms_level)
            peaks.append(p)
        return peaks
