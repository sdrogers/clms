import numpy as np
import pylab as plt
import bisect


class BlockData(object):
    def __init__(self, datasets, mz_step, rt_step, rt_range=[(0, 1450)], mz_range=[(50, 1070)]):
        self.datasets = datasets
        self.mz_step = mz_step
        self.rt_step = rt_step
        self.rt_range = rt_range
        self.mz_range = mz_range
        self.keys = []
        for j in range(len(self.datasets)):
            self.keys.append(list(self.datasets[j].file_spectra.keys()))
        self.mz_bin_lower = np.arange(mz_range[0][0], mz_range[0][1], mz_step)
        self.rt_bin_lower = np.arange(rt_range[0][0], rt_range[0][1], rt_step)
        self.n_mz_bin_lower = len(self.mz_bin_lower)
        self.n_rt_bin_lower = len(self.rt_bin_lower)

        self.intensity_mats = []
        for j in range(len(self.datasets)):
            for i in range(len(self.keys[j])):
                self.intensity_mats.append(self._block_file(i, j))
                print("Processed", self.keys[j][i])

    def _block_file(self, num, data_num):
        intensity_mat = np.zeros((self.n_mz_bin_lower, self.n_rt_bin_lower), np.double)
        spectra = self.datasets[data_num].file_spectra[self.keys[data_num][num]]
        c1 = 0
        for scan_num in spectra:
            scan_rt = spectra[scan_num].scan_time[0]
            if scan_rt < self.rt_bin_lower[0]:
                continue
            if scan_rt > self.rt_bin_lower[-1] + self.rt_step:
                continue
            else:
                rt_pos = bisect.bisect_right(self.rt_bin_lower, scan_rt)
                rt_pos -= 1
                for peak_num in range(len(spectra[scan_num].mz)):
                    mz = spectra[scan_num].mz[peak_num]
                    intensity = spectra[scan_num].i[peak_num]
                    if mz < self.mz_bin_lower[0]:
                        continue
                    if mz > self.mz_bin_lower[-1] + self.mz_step:
                        break
                    mz_pos = bisect.bisect_right(self.mz_bin_lower, mz)
                    mz_pos -= 1

                    intensity_mat[mz_pos, rt_pos] += intensity
        return intensity_mat

    def plot(self, data_num):
        print("Warning: Python prints plots in a stupid stupid way!")
        plt.imshow(np.log(self.intensity_mats[data_num][-1:0:-1] + 1), aspect='auto')

    def combine(self, plot=True):
        combined = []
        for mat in self.intensity_mats:
            combined.append(list(mat.flatten()))
        return combined


def gibbs_sampler(X, observed, R, prior_u, prec_u, prior_v, prec_v, alpha, n_its = 1000, burn_in = 100, true_V=[]):
    # initialise
    N, M = X.shape
    U = np.random.normal(size=(N, R))
    if len(true_V) == 0:
        V = np.random.normal(size=(M, R))
    else:
        V = true_V
    tot_U = np.zeros((N, R))
    tot_V = np.zeros((M, R))
    all_err = []
    for it in range(n_its):
        # loop over u, updating them
        # first compute the covariance - shared if all data observed
        prec_mat = prec_u + alpha * np.dot(V.T, V)
        cov_mat = np.linalg.inv(prec_mat)
        for n in range(N):
            if observed[n, :].sum() < M:
                # not all data observed, compute specific precision
                this_prec_mat = prec_u + alpha * np.dot(np.dot(V.T, np.diag(observed[n, :])), V)
                this_cov_mat = np.linalg.inv(this_prec_mat)
            else:
                this_prec_mat = prec_mat
                this_cov_mat = cov_mat
            s = np.zeros(R)
            for m in range(M):
                if observed[n, m]:
                    s += X[n, m] * V[m, :]
            s *= alpha
            s += np.dot(prec_u, prior_u)
            cond_mu = np.dot(this_cov_mat, s)
            U[n, :] = np.random.multivariate_normal(cond_mu, this_cov_mat)

        # loop over v updating them
        # first covariance
        if len(true_V) == 0:
            prec_mat = prec_v + alpha * np.dot(U.T, U)
            cov_mat = np.linalg.inv(prec_mat)
            for m in range(M):
                if observed[:, m].sum() < N:
                    this_prec_mat = prec_v + alpha * np.dot(np.dot(U.T, np.diag(observed[:, m])), U)
                    this_cov_mat = np.linalg.inv(this_prec_mat)
                else:
                    this_prec_mat = prec_mat
                    this_cov_mat = cov_mat

                s = np.zeros(R)
                for n in range(N):
                    if observed[n, m]:
                        s += X[n, m] * U[n, :]
                s *= alpha
                s += np.dot(prec_v, prior_v)
                cond_mu = np.dot(this_cov_mat, s)
                V[m, :] = np.random.multivariate_normal(cond_mu, this_cov_mat)
        if it >  burn_in:
            tot_U += U
            tot_V += V

        recon_error = np.sqrt(((X - np.dot(U, V.T)) ** 2).mean())
        all_err.append(recon_error)

    if len(true_V) == 0:
        return tot_U / (n_its - burn_in), tot_V / (n_its - burn_in)
    else:
        return tot_U / (n_its - burn_in), true_V

