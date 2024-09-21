import numpy as np
import math
import h5py
import json
import os
from pathlib import Path
from scipy.stats import norm
from lmfit.models import GaussianModel
from scipy.optimize import curve_fit
from scipy.optimize import minimize_scalar
from sklearn.decomposition import PCA
from quantify_core.data.handling import set_datadir, get_tuids_containing, extract_parameter_from_snapshot, load_snapshot, load_quantities_of_interest, load_dataset, to_gridded_dataset, load_processed_dataset

#get parameters from snapshots
def get_RO_freq(tuid):
    snapshot = load_snapshot(tuid)
    RO_freq = extract_parameter_from_snapshot(snapshot, "FX8.clock_freqs.readout")['value']
    return RO_freq

def get_flux_pulse_amp(tuid):
    snapshot = load_snapshot(tuid)
    fpa = extract_parameter_from_snapshot(snapshot, "FX8.measure_flux.flux_amp_1")['value']
    return fpa

def get_integration_time(tuid):
    snapshot = load_snapshot(tuid)
    int_time = extract_parameter_from_snapshot(snapshot, "FX8.measure_flux.integration_time")['value']
    return int_time

def get_fid_discr(tuid):
    qois = load_quantities_of_interest(tuid, 'SSROAnalysis')
    fid = qois['F_discrimination_fit_FX8']
    return fid

def get_ro_amp(tuid):
    snapshot = load_snapshot(tuid)
    ro_amp = extract_parameter_from_snapshot(snapshot, "FX8.measure_flux.pulse_amp")['value']
    return ro_amp

#functions for calculating assignment fidelity (fit) as defined by OQS
def double_norm_cdf(
    x,
    amplitude_0: float,
    center_0: float,
    sigma_0: float,
    amplitude_1: float,
    center_1: float,
    sigma_1: float,
):
    """
    The cumulative density function (cdf) of a double gaussian probability density
    function.

    Parameters
    ----------

    """
    y = amplitude_0 * norm.cdf(x, center_0, sigma_0) + amplitude_1 * norm.cdf(
        x, center_1, sigma_1
    )
    return y


def double_gauss_cdf_diff(
    x,
    amp_0_g: float,
    center_0_g: float,
    sigma_0_g: float,
    amp_1_g: float,
    center_1_g: float,
    sigma_1_g: float,
    amp_0_e: float,
    center_0_e: float,
    sigma_0_e: float,
    amp_1_e: float,
    center_1_e: float,
    sigma_1_e: float,
):
    """
    The difference between two double gaussian cumulative density functions.

    The first (0 or 1) subscript refers to which of the two peaks within the double
    gauss.
    The second (g or e) subscript refers to which of the two distinct curves, g is for
    the qubit prepared in the grounds state, e for the qubit prepared in the excited
    state.
    """
    curve_g = double_norm_cdf(
        x, amp_0_g, center_0_g, sigma_0_g, amp_1_g, center_1_g, sigma_1_g
    )
    curve_e = double_norm_cdf(
        x, amp_0_e, center_0_e, sigma_0_e, amp_1_e, center_1_e, sigma_1_e
    )
    return curve_g - curve_e

# only get raw assignment fidelity from counting shots in each bin
def get_raw_ass_fidelity(tuid):
    f = load_processed_dataset(tuid, analysis_name='SSROAnalysis')
    data_y0 = np.array(list(f["y0"]))
    data_y1 = np.array(list(f["y1"]))

    s0 = data_y0[0,:]+1j*data_y1[0,:]
    s1 = data_y0[1,:]+1j*data_y1[1,:]

    X = np.zeros((2*np.size(s0), 2))
    X[:,0] = np.concatenate((s0.real, s1.real))
    X[:,1] = np.concatenate((s0.imag, s1.imag))

    pca = PCA(n_components=2)
    X_new = pca.fit_transform(X)

    s0 = X_new[0:np.size(s0),0] + 1j*X_new[0:np.size(s0),1]
    s1 = X_new[np.size(s0):,0] + 1j*X_new[np.size(s0):,1]

    signal0 = np.real(s0)
    signal1 = np.real(s1)

    signal0_mean = np.mean(signal0)
    signal1_mean = np.mean(signal1)

    # Make sure rotation puts 0 shots on left and 1 shots on right

    if signal0_mean > signal1_mean:
        # 0 shots right of 1 shots
        signal0, signal1 = -signal0, -signal1

    else:
        signal0, signal1 = signal0, signal1


    signals = np.concatenate((signal0, signal1))
    bin_distance = 0.5
    num_bins = int(np.round(np.abs(np.max(signals)-np.min(signals))/bin_distance))
    num_bins = 150
    hist_range = float(np.min(signals)), float(np.max(signals))
    counts0, bin_edges0 = np.histogram(signal0, bins=num_bins, range=hist_range, density=True)
    counts1, bin_edges1 = np.histogram(signal1, bins=num_bins, range=hist_range, density=True)
    bin_centers0 = (bin_edges0[1:] + bin_edges0[:-1]) / 2
    bin_centers1 = (bin_edges1[1:] + bin_edges1[:-1]) / 2

    bins = bin_centers0
    counts_0 = counts0
    counts_1 = counts1

    # extract raw assignment fidelity (counts shots in each bin)
    cumulative_counts_0 = np.cumsum(counts_0) / np.sum(counts_0)
    cumulative_counts_1 = np.cumsum(counts_1) / np.sum(counts_1)

    max_diff_idx = np.argmax(np.array(cumulative_counts_0 - cumulative_counts_1))
    voltage_threshold = float(bin_centers0[max_diff_idx])
    max_diff = np.array(cumulative_counts_0 - cumulative_counts_1)[max_diff_idx]
    assign_fid = 1 - (1 - max_diff) / 2
    return assign_fid

# gets SNR as defined in previous work: Stefanski, T. V., & Andersen, C. K. (2023). Flux-pulse-assisted Readout of a Fluxonium Qubit. arXiv preprint arXiv:2309.17286.
def get_paperdef_snr_error(tuid):
    qois = load_quantities_of_interest(tuid, 'SSROAnalysis')
    mu0 = qois['fit_result']['g_single_fit_FX8']['params']['center']['value']
    mu1 = qois['fit_result']['e_single_fit_FX8']['params']['center']['value']
    sigma0 = qois['fit_result']['g_single_fit_FX8']['params']['sigma']['value']
    sigma1 = qois['fit_result']['e_single_fit_FX8']['params']['sigma']['value']
    signal = np.abs(mu0-mu1)
    noise = np.sqrt(sigma0**2+sigma1**2)
    snr = signal/noise
    error = (math.erfc(snr/2))/2
    return snr, error

# fitting functions for measurement efficiency

def gaussian(x, amp, mu, sigma):
    return (amp/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-mu)**2)/(2*sigma**2))

def lin_fit(e, a):
    return a*e

def coherence_decay(epsilon, b, sigma):
    return b*np.exp(-epsilon**2/(2*sigma**2))