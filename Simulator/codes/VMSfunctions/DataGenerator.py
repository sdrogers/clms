import glob
import math
import os
from collections import defaultdict

import numpy as np
import pylab as plt
import pymzml
from sklearn.neighbors import KernelDensity
import pandas as pd
import math

from VMSfunctions.Chromatograms import EmpiricalChromatogram, UnknownChemical
from VMSfunctions.Common import *


class Peak(object):
    """
    A simple class to represent an empirical or sampled scan-level peak object
    """

    def __init__(self, mz, rt, intensity, ms_level):
        self.mz = mz
        self.rt = rt
        self.intensity = intensity
        self.ms_level = ms_level

    def __repr__(self):
        return 'Peak mz=%.4f rt=%.2f intensity=%.2f ms_level=%d' % (self.mz, self.rt, self.intensity, self.ms_level)

    def __eq__(self, other):
        if not isinstance(other, Peak):
            # don't attempt to compare against unrelated types
            return NotImplemented

        return math.isclose(self.mz, other.mz) and \
            math.isclose(self.rt, other.rt) and \
            math.isclose(self.intensity, other.intensity) and \
            self.ms_level == other.ms_level


class RegionOfInterest(object):
    def __init__(self, file_name, mode, pickedPeak, mzrange, rtrange, scrange):
        self.file_name = file_name
        self.mode = mode
        self.pickedPeak = pickedPeak
        self.mzrange = mzrange
        self.rtrange = rtrange
        self.scrange = scrange
        self.peaks = []

    def add(self, p):
        if (self.mzrange[0] <= p.mz <= self.mzrange[1]) and (self.rtrange[0] <= p.rt <= self.rtrange[1]):
            self.peaks.append(p)

    def num_scans(self):
        return self.scrange[1] - self.scrange[0] + 1

    def mzs(self):
        return np.array([p.mz for p in self.peaks])

    def rts(self):
        return np.array([p.rt for p in self.peaks])

    def intensities(self):
        return np.array([p.intensity for p in self.peaks])

    def to_chromatogram(self):
        if len(self.peaks) == 0 or len(self.peaks) == 1:
            return None
        chrom = EmpiricalChromatogram(self.rts(), self.mzs(), self.intensities(), roi=self)
        return chrom

    def __repr__(self):
        return 'ROI %s %s picked %s mz (%.4f-%.4f) rt (%.4f-%.4f) scans (%d-%d)' % (
        self.file_name, self.mode, self.pickedPeak,
        self.mzrange[0], self.mzrange[1],
        self.rtrange[0], self.rtrange[1],
        self.scrange[0], self.scrange[1])


class DataSource(LoggerMixin):
    """
    A class to load and extract centroided peaks from CSV and mzML files.
    :param min_ms1_intensity: minimum ms1 intensity for filtering
    :param min_ms2_intensity: maximum ms2 intensity for filtering
    :param min_rt: minimum RT for filtering
    :param max_rt: maximum RT for filtering
    """

    def __init__(self):
        # A dictionary that stores the actual pymzml spectra for each filename
        self.file_spectra = {} # key: filename, value: a dict where key is scan_number and value is spectrum

        # A dictionary to store the distribution on scan durations for each ms_level
        self.scan_durations = defaultdict(list)

        # A dictionary to stores region of interests
        self.all_rois = {}

        # Dataframe of exported ROI information
        self.roi_df = None

        # pymzml parameters
        self.ms1_precision = 5e-6
        self.obo_version = '4.0.1'

    def load_data(self, mzml_path):
        """
        Loads data and generate peaks from mzML files. The resulting peak objects will not have chromatographic peak
        shapes, because no peak picking has been performed yet.
        :param mzml_path: the input folder containing the mzML files
        :return: nothing, but the instance variable file_spectra and scan_durations are populated
        """
        file_scan_durations = {} # key: filename, value: a dict where key is ms level and value is scan durations
        file_spectra = {} # key: filename, value: a dict where key is scan_number and value is spectrum
        for filename in glob.glob(os.path.join(mzml_path, '*.mzML')):
            run = pymzml.run.Reader(filename, obo_version=self.obo_version,
                                    MS1_Precision=self.ms1_precision,
                                    extraAccessions=[('MS:1000016', ['value', 'unitName'])])

            fname = os.path.basename(filename)
            self.logger.info('Loading %s' % fname)
            file_spectra[fname] = {}
            file_scan_durations[fname] = {
                1: [],
                2: []
            }
            start_time = 0
            for scan_number, spectrum in enumerate(run):
                file_spectra[fname][scan_number] = spectrum

                # get spectrum ms level and retention time
                ms_level = spectrum.ms_level
                rt = self._get_rt(spectrum)

                # store the scan duration of each spectrum
                end_time = rt
                duration = end_time - start_time
                item = [start_time, duration]
                file_scan_durations[fname][ms_level].append(item)
                start_time = end_time

        self.file_scan_durations = file_scan_durations
        self.file_spectra = file_spectra

    def plot_histogram(self, X, data_type, bins=100):
        """
        Makes a histogram plot on the distribution of the item of interest
        :param X: a numpy array
        :param bins: number of histogram bins
        :return: nothing. A plot is shown.
        """
        plt.figure()
        _ = plt.hist(X, bins=bins)
        plt.plot(X[:, 0], np.full(X.shape[0], -0.01), '|k')
        plt.title('Histogram for %s -- shape %s' % (data_type, str(X.shape)))
        plt.show()

    def plot_boxplot(self, X, data_type):
        """
        Makes a boxplot on the distribution of the item of interest
        :param X: a numpy array
        :return: nothing. A plot is shown.
        """
        plt.figure()
        _ = plt.boxplot(X)
        plt.title('Boxplot for %s -- shape %s' % (data_type, str(X.shape)))
        plt.show()

    def plot_peak(self, peak):
        f, axarr = plt.subplots(2, sharex=True)
        axarr[0].plot(peak.rt_values, peak.intensity_values)
        axarr[1].plot(peak.rt_values, peak.mz_values, linestyle='None', marker='o', markersize=1.0, color='b')

    def get_data(self, data_type, filename, ms_level, min_intensity=None,
                 min_rt=None, max_rt=None, log=False, max_data=100000):
        """
        Retrieves values as numpy array
        :param data_type: data_type is 'mz', 'rt', 'intensity' or 'n_peaks'
        :param filename: the mzml filename or None for all files
        :param ms_level: level 1 or 2
        :param min_intensity: minimum ms2 intensity for thresholding
        :param min_rt: minimum RT value for thresholding
        :param max_rt: max RT value for thresholding
        :param log: if true, the returned values will be logged
        :return: an Nx1 numpy array of all the values requested
        """
        if data_type == SCAN_DURATION:
            if filename is None: # use the scan durations from all files
                values = []
                for f in self.file_scan_durations:
                    temp = self._get_filtered_scan_durations(f, ms_level, min_rt, max_rt)
                    values.extend(temp)
            else: # use the scan durations from that file only
                values = self._get_filtered_scan_durations(filename, ms_level, min_rt, max_rt)
        else:
            # get spectra from either one file or all files
            if filename is None: # use all spectra
                all_spectra = []
                for f in self.file_spectra:
                    spectra_for_f = list(self.file_spectra[f].values())
                    all_spectra.extend(spectra_for_f)
            else: # use spectra for that file only
                all_spectra = self.file_spectra[filename].values()

            # loop through spectrum and get all peaks above threshold
            values = []
            for spectrum in all_spectra:
                # if wrong ms level, skip this spectrum
                if spectrum.ms_level != ms_level:
                    continue

                # collect all valid Peak objects in a spectrum
                spectrum_peaks = []
                for mz, intensity in spectrum.peaks('raw'):
                    rt = self._get_rt(spectrum)
                    p = Peak(mz, rt, intensity, spectrum.ms_level)
                    if self._valid_peak(p, min_intensity, min_rt, max_rt):
                        spectrum_peaks.append(p)

                if data_type == N_PEAKS:
                    n_peaks = len(spectrum_peaks)
                    if n_peaks > 0:
                        values.append(n_peaks)
                elif data_type == MZ_INTENSITY:
                    mzs = list(getattr(x, MZ) for x in spectrum_peaks)
                    intensities = list(getattr(x, INTENSITY) for x in spectrum_peaks)
                    values.extend(list(zip(mzs, intensities)))
                else: # MZ, INTENSITY or RT
                    attrs = list(getattr(x, data_type) for x in spectrum_peaks)
                    values.extend(attrs)

        # log-transform if necessary
        X = np.array(values)
        if log:
            if data_type == MZ_INTENSITY: # just log the intensity part
                X[:, 1] = np.log(X[:, 1])
            else:
                X = np.log(X)

        # pick random samples
        try:
            idx = np.arange(len(X))
            rnd_idx = np.random.choice(idx, size=int(max_data), replace=False)
            sampled_X = X[rnd_idx]
        except ValueError:
            sampled_X = X

        # return values
        if data_type == MZ_INTENSITY:
            return sampled_X # it's already a Nx2 array
        else:
            # convert into Nx1 array
            return sampled_X[:, np.newaxis]

    def _get_filtered_scan_durations(self, filename, ms_level, min_rt, max_rt):
        temp = self.file_scan_durations[filename][ms_level]
        if min_rt is not None:
            temp = filter(lambda item: min_rt < item[0], temp)
        if max_rt is not None:
            temp = filter(lambda item: item[0] < max_rt, temp)
        values = list(temp.map(lambda item: item[1]))
        return values

    def extract_roi(self, roi_file, min_rt=None, max_rt=None, filename=None):
        self.roi_df = pd.read_csv(roi_file)
        if filename is None:
            roi_filenames = self.roi_df['file'].unique()
        else:
            roi_filenames = np.array([filename])

        # convert each row in the dataframe to a ROI object
        for fname in roi_filenames:
            self.logger.info('Creating ROI objects for %s' % fname)
            rois_data = {
                'rois': [],
                'mzmin': [],
                'mzmax': [],
                'rtmin': [],
                'rtmax': []
            }

            # convert each row of the dataframe to roi objects
            df = self.roi_df[self.roi_df['file'] == fname]
            for idx, row in df.iterrows(): # TODO: make this faster
                if (idx % 10000 == 0):
                    self.logger.debug('%6d/%6d' % (len(rois_data['rois']), df.shape[0]))
                file_name = row['file']
                mzmin = row['mzmin']
                mzmax = row['mzmax']
                rtmin = row['rtmin']
                rtmax = row['rtmax']
                scmin = row['scmin']
                scmax = row['scmax']
                pickedPeak = row['pickedPeak']
                mode = row['mode']
                roi = RegionOfInterest(file_name, mode, pickedPeak, (mzmin, mzmax), (rtmin, rtmax), (scmin, scmax))
                if self._valid_roi(roi, min_rt, max_rt):
                    rois_data['rois'].append(roi)
                    rois_data['mzmin'].append(mzmin)
                    rois_data['mzmax'].append(mzmax)
                    rois_data['rtmin'].append(rtmin)
                    rois_data['rtmax'].append(rtmax)

            # convert all values to numpy arrays
            rois_data['rois'] = np.array(rois_data['rois'])
            rois_data['mzmin'] = np.array(rois_data['mzmin'])
            rois_data['mzmax'] = np.array(rois_data['mzmax'])
            rois_data['rtmin'] = np.array(rois_data['rtmin'])
            rois_data['rtmax'] = np.array(rois_data['rtmax'])
            self.all_rois[fname] = rois_data
            self.logger.info('Extracted %d ROIs for %s' % (len(rois_data['rois']), fname))

    def _valid_roi(self, roi, min_rt=None, max_rt=None):
        if min_rt is not None and roi.rtrange[0] < min_rt:
            return False
        if max_rt is not None and roi.rtrange[1] > max_rt:
            return False
        return True

    def populate_roi(self, filename=None):
        if filename is None:
            roi_filenames = self.roi_df['file'].unique()
        else:
            roi_filenames = np.array([filename])

        # assign raw spectrum peaks to ROI
        for fname in roi_filenames:
            rois_data = self.all_rois[fname]
            self.logger.info('Populating ROI objects for %s' % fname)

            # get spectra for a file
            spectra = self.file_spectra[fname]
            for scan_id, spectrum in spectra.items(): # loop over all scans
                if scan_id % 100 == 0:
                    self.logger.debug('%6d/%6d processing spectrum %s' % (scan_id, len(spectra), spectrum))
                rt = self._get_rt(spectrum)
                for mz, intensity in spectrum.peaks('raw'):
                    # find the ROIs that contain this spectrum peak
                    p = Peak(mz, rt, intensity, spectrum.ms_level)
                    rois = self._get_containing_rois(p, rois_data)
                    for roi in rois: # if found, assign peaks to ROIs
                        roi.add(p)

    def save_roi(self, out_file):
        save_obj(self.all_rois, out_file)

    def load_roi(self, in_file):
        self.all_rois = load_obj(in_file)

    def _get_rt(self, spectrum):
        rt, units = spectrum.scan_time
        if units == 'minute':
            rt *= 60.0
        return rt

    def _valid_peak(self, peak, min_intensity, min_rt, max_rt):
        if min_intensity is not None and peak.intensity < min_intensity:
            return False
        elif min_rt is not None and peak.rt < min_rt:
            return False
        elif max_rt is not None and peak.rt > max_rt:
            return False
        else:
            return True

    def _get_containing_rois(self, peak, rois_data):
        mzmin_check = rois_data['mzmin'] <= peak.mz
        mzmax_check = peak.mz <= rois_data['mzmax']
        rtmin_check = rois_data['rtmin'] <= peak.rt
        rtmax_check = peak.rt <= rois_data['rtmax']
        idx = np.nonzero(mzmin_check & mzmax_check & rtmin_check & rtmax_check)[0]
        rois = rois_data['rois'][idx]
        return rois


class DensityEstimator(LoggerMixin):
    """A class to perform kernel density estimation. Takes as input a DataSource class."""

    def __init__(self, min_ms1_intensity, min_ms2_intensity, min_rt, max_rt, plot=False):
        self.kdes = {}
        self.kernel = 'gaussian'
        self.min_ms1_intensity = min_ms1_intensity
        self.min_ms2_intensity = min_ms2_intensity
        self.min_rt = min_rt
        self.max_rt = max_rt
        self.plot = plot

    def kde(self, data_source, filename, ms_level, bandwidth_mz=1.0, bandwidth_intensity=0.01,
            bandwidth_rt=5.0, bandwidth_n_peaks=1.0, bandwidth_scan_durations=0.01, max_data=100000):
        params = [
            {'data_type': MZ, 'bandwidth': bandwidth_mz},
            {'data_type': INTENSITY, 'bandwidth': bandwidth_intensity},
            {'data_type': RT, 'bandwidth': bandwidth_rt},
            {'data_type': N_PEAKS, 'bandwidth': bandwidth_n_peaks},
            {'data_type': SCAN_DURATION, 'bandwidth': bandwidth_scan_durations},
        ]
        for param in params:
            data_type = param['data_type']

            # get data
            self.logger.debug('Retrieving %s values from %s' % (data_type, data_source))
            min_intensity = self.min_ms1_intensity if ms_level == 1 else self.min_ms2_intensity
            log = True if data_type == INTENSITY else False
            X = data_source.get_data(data_type, filename, ms_level, min_intensity=min_intensity,
                                     min_rt=self.min_rt, max_rt=self.max_rt, log=log, max_data=max_data)
            self.logger.debug('Retrieved %s values' % str(X.shape))

            # fit kde
            self.logger.debug('Fitting kde on %s' % data_type)
            bandwidth = param['bandwidth']
            kde = KernelDensity(kernel=self.kernel, bandwidth=bandwidth).fit(X)
            self.kdes[(data_type, ms_level)] = kde

            # plot if necessary
            if self.plot:
                self._plot(kde, X, data_type, filename, bandwidth)

    def sample(self, ms_level, n_sample):
        """
        Samples mz, rt and intensity values
        :param ms_level: the ms level
        :param n_sample: the number of samples to draw
        :return:
        """
        mz_vals = self.kdes[(MZ, ms_level)].sample(n_sample)
        intensity_vals = self.kdes[(INTENSITY, ms_level)].sample(n_sample)
        rt_vals = self.kdes[(RT, ms_level)].sample(n_sample)
        vals = np.concatenate((mz_vals, intensity_vals, rt_vals), axis=1)
        return vals

    def n_peaks(self, ms_level, n_sample):
        return self.kdes[(N_PEAKS, ms_level)].sample(n_sample)

    def scan_durations(self, ms_level, n_sample):
        return self.kdes[(SCAN_DURATION, ms_level)].sample(n_sample)

    def _plot(self, kde, X, data_type, filename, bandwidth):
        if self.plot:
            fname = 'All' if filename is None else filename
            title = '%s density estimation for %s - bandwidth %.3f' % (data_type, fname, bandwidth)
            X_plot = np.linspace(np.min(X), np.max(X), 1000)[:, np.newaxis]
            log_dens = kde.score_samples(X_plot)
            plt.figure()
            plt.fill_between(X_plot[:, 0], np.exp(log_dens), alpha=0.5)
            plt.plot(X[:, 0], np.full(X.shape[0], -0.01), '|k')
            plt.title(title)
            plt.show()


class PeakDensityEstimator(LoggerMixin):
    """A class to perform better kernel density estimation for peak data. Takes as input a DataSource class."""

    def __init__(self, min_ms1_intensity, min_ms2_intensity, min_rt, max_rt, plot=False):
        self.kdes = {}
        self.kernel = 'gaussian'
        self.min_ms1_intensity = min_ms1_intensity
        self.min_ms2_intensity = min_ms2_intensity
        self.min_rt = min_rt
        self.max_rt = max_rt
        self.plot = plot

    def kde(self, data_source, filename, ms_level, bandwidth_mz_intensity=1.0,
            bandwidth_rt=5.0, bandwidth_n_peaks=1.0, bandwidth_scan_durations=0.01, max_data=100000):
        params = [
            {'data_type': MZ_INTENSITY, 'bandwidth': bandwidth_mz_intensity},
            {'data_type': RT, 'bandwidth': bandwidth_rt},
            {'data_type': N_PEAKS, 'bandwidth': bandwidth_n_peaks},
            {'data_type': SCAN_DURATION, 'bandwidth': bandwidth_scan_durations},
        ]
        for param in params:
            data_type = param['data_type']

            # get data
            self.logger.debug('Retrieving %s values from %s' % (data_type, data_source))
            min_intensity = self.min_ms1_intensity if ms_level == 1 else self.min_ms2_intensity
            log = True if data_type == MZ_INTENSITY else False
            X = data_source.get_data(data_type, filename, ms_level, min_intensity=min_intensity,
                                      min_rt=self.min_rt, max_rt=self.max_rt, log=log, max_data=max_data)

            # fit kde
            bandwidth = param['bandwidth']
            kde = KernelDensity(kernel=self.kernel, bandwidth=bandwidth).fit(X)
            self.kdes[(data_type, ms_level)] = kde

            # plot if necessary
            self._plot(kde, X, data_type, filename, bandwidth)

    def sample(self, ms_level, n_sample):
        mz_intensity_vals = self.kdes[(MZ_INTENSITY, ms_level)].sample(n_sample)
        rt_vals = self.kdes[(RT, ms_level)].sample(n_sample)
        vals = np.concatenate((mz_intensity_vals, rt_vals), axis=1)
        return vals

    def n_peaks(self, ms_level, n_sample):
        return self.kdes[(N_PEAKS, ms_level)].sample(n_sample)

    def scan_durations(self, ms_level, n_sample):
        return self.kdes[(SCAN_DURATION, ms_level)].sample(n_sample)

    def _plot(self, kde, X, data_type, filename, bandwidth):
        if self.plot:
            if data_type == MZ_INTENSITY:
                self.logger.debug('Plotting for %s not implemented' % MZ_INTENSITY)
            else:
                fname = 'All' if filename is None else filename
                title = '%s density estimation for %s - bandwidth %.3f' % (data_type, fname, bandwidth)
                X_plot = np.linspace(np.min(X), np.max(X), 1000)[:, np.newaxis]
                log_dens = kde.score_samples(X_plot)
                plt.figure()
                plt.fill_between(X_plot[:, 0], np.exp(log_dens), alpha=0.5)
                plt.plot(X[:, 0], np.full(X.shape[0], -0.01), '|k')
                plt.title(title)
                plt.show()


class PeakSampler(LoggerMixin):
    """A class to sample peaks from a trained density estimator"""

    def __init__(self, density_estimator):
        self.density_estimator = density_estimator

    def sample(self, ms_level, n_peaks=None, min_mz=None, max_mz=None, min_rt=None, max_rt=None, min_intensity=None):
        if n_peaks is None:
            n_peaks = max(self.density_estimator.n_peaks(ms_level, 1).astype(int)[0][0], 0)

        peaks = []
        while len(peaks) < n_peaks:
            vals = self.density_estimator.sample(ms_level, 1)
            intensity = np.exp(vals[0, 1])
            mz = vals[0, 0]
            rt = vals[0, 2]
            p = Peak(mz, rt, intensity, ms_level)
            if self._is_valid(p, min_mz, max_mz, min_rt, max_rt, min_intensity): # othwerise we just keep rejecting
                peaks.append(p)
        return peaks

    def _is_valid(self, peak, min_mz, max_mz, min_rt, max_rt, min_intensity):
        if peak.intensity < 0:
            return False
        if min_mz is not None and min_mz > peak.mz:
            return False
        if max_mz is not None and max_mz < peak.mz:
            return False
        if min_rt is not None and min_rt > peak.rt:
            return False
        if max_rt is not None and max_rt < peak.rt:
            return False
        if min_intensity is not None and min_intensity > peak.intensity:
            return False
        return True