import glob
import math
from collections import defaultdict

import pandas as pd
import numpy as np
import os
import pylab as plt
import pymzml
from sklearn.neighbors import KernelDensity

from .Common import MZ, RT, INTENSITY, N_PEAKS, MZ_INTENSITY, SCAN_DURATION


class PeakSample(object):
    """
    A simple class to represent an empirical or sampled scan-level peak object
    """

    def __init__(self, mz, rt, intensity, ms_level):
        self.mz = mz
        self.rt = rt
        self.intensity = intensity
        self.ms_level = ms_level

    def __repr__(self):
        return 'PeakSample mz=%.4f rt=%.2f intensity=%.2f ms_level=%d' % (self.mz, self.rt, self.intensity, self.ms_level)



class DataSource(object):
    """
    A class to load and extract centroided peaks from CSV and mzML files.
    :param min_ms1_intensity: minimum ms1 intensity for filtering
    :param min_ms2_intensity: maximum ms2 intensity for filtering
    :param min_rt: minimum RT for filtering
    :param max_rt: maximum RT for filtering
    """

    def __init__(self, min_ms1_intensity=0, min_ms2_intensity=0, min_rt=0, max_rt=math.inf):
        # A dictionary that stores scan-level ms1 and ms2 peaks
        self.peaks = defaultdict(list)
        self.noisy_peaks = defaultdict(list)

        # A dictionary to store the distribution on the number of peaks per scan for each ms_level
        self.num_peaks = defaultdict(list)

        # pymzml parameters
        self.ms1_precision = 5e-6
        self.obo_version = '4.0.1'

        # filtering parameters
        self.min_ms1_intensity = min_ms1_intensity
        self.min_ms2_intensity = min_ms2_intensity
        self.min_rt = min_rt
        self.max_rt = max_rt

    def load_data(self, ms1_df, mzml_path):
        """
        Loads data and generate peaks from mzML files. The resulting peak objects will not have chromatographic peak
        shapes, because no peak picking has been performed yet.
        :param data_path: the input folder containing the mzML files
        :param xcms_df: the path to a dataframe (in CSV) format containing chromatographic peak information
        :return: nothing, but the instance variables self.peaks and self.num_peaks will be populated by the peak objects
        at different ms level.
        """
        scan_peaks, scan_num_peaks, scan_noisy_peaks, scan_durations = self._load_mzml(mzml_path)
        self.peaks = scan_peaks
        self.noisy_peaks = scan_noisy_peaks
        self.num_peaks = scan_num_peaks
        self.scan_durations = scan_durations

    def plot_histogram(self, data_type, ms_level, log=False, bins=100):
        """
        Makes a histogram plot on the distribution of the item of interest
        :param data_type: data_type is 'mz', 'rt', 'intensity' or 'n_peaks'
        :param ms_level: level 1 or 2
        :param log: whether to log-transform the data
        :param bins: number of histogram bins
        :return: nothing. A plot is shown.
        """
        X = self.get_data(data_type, ms_level, log)
        plt.figure()
        _ = plt.hist(X, bins=bins)
        plt.plot(X[:, 0], np.full(X.shape[0], -0.01), '|k')
        plt.title('%s histogram' % data_type)
        plt.show()

    def plot_peak(self, peak):
        f, axarr = plt.subplots(2, sharex=True)
        axarr[0].plot(peak.rt_values, peak.intensity_values)
        axarr[1].plot(peak.rt_values, peak.mz_values, linestyle='None', marker='o', markersize=1.0, color='b')

    def get_data(self, data_type, ms_level, log=False):
        """
        Retrieves values as numpy array
        :param data_type: data_type is 'mz', 'rt', 'intensity' or 'n_peaks'
        :param ms_level: level 1 or 2
        :return: an Nx1 numpy array of all the values from the data
        """
        if data_type == N_PEAKS:  # get observed peak counts across all spectra in the data
            values = self.num_peaks[ms_level]
        elif data_type == SCAN_DURATION:
            values = self.scan_durations[ms_level]
        else:  # get observed mz, rt or intensity values of all peaks in the data
            peaks = self.peaks[ms_level]
            values = list(getattr(x, data_type) for x in peaks)

        # convert into Nx1 array and also log-transform if necessary
        X = np.array(values)[:, np.newaxis]
        if log:
            X = np.log(X)
        return X

    def _load_mzml(self, data_path):
        """
        Load MS1 and MS2 peaks from scan data from all mzML files in data_path
        :param data_path: the path containing mzML files to load
        :return: the tuple of two dictionaries, where key corresponds to the ms level (1 or 2):
           - scan_peaks is a dictionary of all peaks across all scans
           - scan_num_peaks is a dictionary of the number of peaks for each scan
        """
        scan_peaks = defaultdict(list)  # good quality peaks across all scans
        noisy_peaks = defaultdict(list)  # noisy peaks across all scans (outside threshold)
        scan_num_peaks = defaultdict(list)  # number of peaks per scan for each ms_level
        scan_durations = defaultdict(list) # durations per scan for each ms level
        for filename in glob.glob(os.path.join(data_path, '*.mzML')):
            run = pymzml.run.Reader(filename, obo_version=self.obo_version,
                                    MS1_Precision=self.ms1_precision,
                                    extraAccessions=[('MS:1000016', ['value', 'unitName'])])

            total_peaks = defaultdict(int)
            start_time = 0
            for scan_number, spectrum in enumerate(run):
                rt, units = spectrum.scan_time
                if units == 'minute':
                    rt *= 60.0

                # loop through spectrum and get all peaks above threshold
                spectrum_peaks = defaultdict(list)
                for mz, intensity in spectrum.peaks('centroided'):
                    ms_level = spectrum.ms_level
                    p = PeakSample(mz, rt, intensity, ms_level)
                    if self._valid_peak(p):
                        spectrum_peaks[ms_level].append(p)
                    else:
                        noisy_peaks[ms_level].append(p)

                # store the results from each spectrum
                end_time = rt
                duration = end_time - start_time
                scan_durations[ms_level].append(duration)
                start_time = end_time

                for ms_level in spectrum_peaks:
                    n = len(spectrum_peaks[ms_level])
                    if n > 0:
                        scan_peaks[ms_level].extend(spectrum_peaks[ms_level])
                        scan_num_peaks[ms_level].append(n)
                        total_peaks[ms_level] += n
            print('%s (ms1=%d, ms2=%d)' % (filename, total_peaks[1], total_peaks[2]))
        return scan_peaks, scan_num_peaks, noisy_peaks, scan_durations

    def _valid_peak(self, peak):
        ms_level = peak.ms_level
        min_intensity = self.min_ms1_intensity if ms_level == 1 else self.min_ms2_intensity
        if peak.intensity < min_intensity:
            return False
        elif peak.rt < self.min_rt:
            return False
        elif peak.rt > self.max_rt:
            return False
        return True


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

    def scan_durations(self, ms_level, n_sample):
        return self.kdes[(SCAN_DURATION, ms_level)].sample(n_sample)

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

    def kde(self, data_source, ms_level, bandwidth_mz_int=1.0, bandwidth_rt=5.0, bandwidth_n_peaks=1.0,
            bandwidth_scan_durations=0.01):
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

        # train kde on scan durations
        X = data_source.get_data(SCAN_DURATION, ms_level)
        self.kdes[(SCAN_DURATION, ms_level)] = KernelDensity(kernel=self.kernel, bandwidth=bandwidth_scan_durations).fit(X)


    def sample(self, ms_level, n_sample):
        mz_intensity_vals = self.kdes[(MZ_INTENSITY, ms_level)].sample(n_sample)
        rt_vals = self.kdes[(RT, ms_level)].sample(n_sample)
        vals = np.concatenate((mz_intensity_vals, rt_vals), axis=1)
        return vals

    def n_peaks(self, ms_level, n_sample):
        return self.kdes[(N_PEAKS, ms_level)].sample(n_sample)

    def scan_durations(self, ms_level, n_sample):
        return self.kdes[(SCAN_DURATION, ms_level)].sample(n_sample)


class PeakSampler(object):
    """A class to sample peaks from a trained density estimator"""

    def __init__(self, density_estimator):
        self.density_estimator = density_estimator

    def sample(self, ms_level, n_peaks=None):
        if n_peaks is None:
            n_peaks = max(self.density_estimator.n_peaks(ms_level, 1).astype(int)[0][0], 0)

        peaks = []
        while (len(peaks) < n_peaks):
            vals = self.density_estimator.sample(ms_level, 1)
            intensity = np.exp(vals[0, 1])
            if intensity > 0:
                mz = vals[0, 0]
                rt = vals[0, 2]
                p = PeakSample(mz, rt, intensity, ms_level)
                peaks.append(p)
        return peaks