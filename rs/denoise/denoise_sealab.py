import nibabel as nib
from nilearn import signal
from nilearn.signal import clean
import numpy as np
import pandas as pd
from scipy import stats
from scipy.interpolate import CubicSpline
import subprocess
import statsmodels.api as sm


"""
# Denoise functions being used in SEA Lab -- for rs fMRI data
# 1. Resample to 32K; 
# 2. Plot gray plots
# 3. Compute new global signal and new DVARS
# 4. For each participant/run with mean FD >0.2mm, z-score DVARS and identify FD of volumes with >2 SD change in DVARS
"""


def wb_resample_32k(ts, out_ts, atlas_dlabel='Gordon333_SeitzmanSubcortical.32k_fs_LR.dlabel.nii'):
    """
    This is a wrapper of wb_command -cifti-resample
    Input Parameters
    ----------
    ts: full path to time series data. i.e., dtseries data (should work for ptseries data as well)
        string
    out_ts: full path of output time series data. i.e., dtseries data (should work for ptseries data as well)
            string
    atlas_dlabel: full path to dlabel data; Default: Gordon333_SeitzmanSubcortical.32k_fs_LR.dlabel.nii
        string
    Returns
    ----------
    None
    """
    subprocess.call(f'wb_command -cifti-resample {ts} COLUMN {atlas_dlabel} COLUMN ADAP_BARY_AREA CUBIC {out_ts}',
                    shell=True)


def rescale_signal(ts, out_ts, method='mean', scale=1000):
    """
    Normalize signal. Default to normalizing signal to a mean of 1000.
    Input Parameters
    ----------
    ts: full path to time series data. i.e., dtseries data (should work for ptseries data as well)
        string
    out_ts: full path of output time series data. i.e., dtseries data (should work for ptseries data as well)
            string
    method: method used to rescale, either 'mean' or 'mode'. Default to 'mean'.
            string
    scale: used scale. Default to 1000.
    Returns
    ----------
    None
    """
    # pull timeseries
    img = nib.load(ts)
    img_pre_data = img.get_fdata()

    if method == 'mean':
        data = (img_pre_data / np.mean(img_pre_data)) * scale
    elif method == 'mode':
        mode_array = stats.mode(img_pre_data, axis=None)
        mode_value = mode_array.mode[0]
        data = (img_pre_data / mode_value) * scale
    else:
        print('=======ERROR======= def rescale_signal; if method == mean:')

    img_ss = nib.cifti2.cifti2.Cifti2Image(data, (img.header.get_axis(0), img.header.get_axis(1)))
    nib.save(img_ss, out_ts)


def make_noise_matrix(fmriprep_confounds_tsv_file, out_file, headsize=50, fd_thres=None):
    """
    36 nuisance regressors were selected from the preprocessing confounds, according to the ‘36P’ strategy.
    These nuisance regressors included six motion parameters, mean global signal, mean white matter signal,
    mean CSF signal with their temporal derivatives, and the quadratic expansion of six motion parameters,
    tissues signals and their temporal derivatives (Ciric et al. 2017; Satterthwaite et al. 2013).
    If needed, later can expand to other nuisance regressors, like 24 etc.
    Input Parameters
    ----------
    fmriprep_confounds_tsv_file: full path to confounds_tsv_file being generated by fmriprep
        string
    out_file: full path of output txt file. this funciton will save a txt file for the nuisance regressors
            string
    headsize: headsize. Default to 50 (adult). The recommended value for infants is 35.
              int
    fd_thres: fd threshold to despike. Default to None (will not despike).
              None or float.
    Returns
    ----------
    None
    """

    # load movement metrics
    nuissance = pd.read_csv(fmriprep_confounds_tsv_file, sep='\t')
    # for debug: nuissance = pd.read_csv('/Volumes/ppm/PPM/NCAP_bids/Data/derivatives/fmriprep/sub-10001scan1/
    # ses-1/func/sub-10001scan1_ses-1_task-restPA_run-2_desc-confounds_timeseries.tsv', sep='\t')
    nuissance = nuissance.fillna(0)

    regressors = nuissance.loc[:, ['global_signal', 'global_signal_derivative1',
                                   'global_signal_power2', 'global_signal_derivative1_power2',
                                   'csf', 'csf_derivative1',
                                   'csf_derivative1_power2', 'csf_power2',
                                   'white_matter', 'white_matter_derivative1',
                                   'white_matter_power2', 'white_matter_derivative1_power2',
                                   'trans_x', 'trans_x_derivative1',
                                   'trans_x_derivative1_power2', 'trans_x_power2',
                                   'trans_y', 'trans_y_derivative1',
                                   'trans_y_derivative1_power2', 'trans_y_power2',
                                   'trans_z', 'trans_z_derivative1',
                                   'trans_z_power2', 'trans_z_derivative1_power2',
                                   'rot_x', 'rot_x_derivative1',
                                   'rot_x_derivative1_power2', 'rot_x_power2',
                                   'rot_y', 'rot_y_derivative1',
                                   'rot_y_power2', 'rot_y_derivative1_power2',
                                   'rot_z', 'rot_z_derivative1',
                                   'rot_z_power2', 'rot_z_derivative1_power2']]

    if fd_thres is None:
        np.savetxt(out_file, regressors.to_numpy())
    else:
        # pull FD
        fd = nuissance['framewise_displacement'].to_numpy()

        # create timeseries of volumes to censor
        vols_to_censor = fd > fd_thres
        n_vols = np.sum(vols_to_censor)
        if n_vols > 0:
            spikes = np.zeros((len(fd), n_vols))
            b = 0
            for a in range(0, len(fd)):
                if vols_to_censor[a] == 1:
                    spikes[a, b] = 1
                    b = b + 1

            regressors_np = np.hstack((regressors.to_numpy(), spikes))
            np.savetxt(out_file, regressors_np)
        else:
            np.savetxt(out_file, regressors.to_numpy())


def regress_out_nuissance(noise_file, standard_ts, regressed_ts_filename):
    """
    regressed out the nuissance from the noise_file file. Save the resid into time series file.
    Input Parameters
    ----------
    noise_file: full path to noise_file.txt file (saved in 'make_noise_matrix' function)
                string
    standard_ts: full path of time series data that needs to be regressed out the nuissance
                 string
    regressed_ts_filename: full path of output time series data (the resid data)
                           string
    Returns
    ----------
    None
    """
    # load nuissance regressors
    noise_mat = np.loadtxt(noise_file)

    # intercept terms
    fit_term = sm.add_constant(noise_mat)

    # load data for debug: standard_ts =
    # '/Volumes/ppm/PPM/NCAP_bids_xcp/Data/derivatives/fmriprep/sub-10001scan1/ses-1/func/sub-10001scan1_ses-1_task-
    # restPA_run-1_space-fsLR_den-91k_bold.32k_fs_LR.dtseries.nii'
    func_data_rescaled = nib.load(standard_ts).get_fdata()

    # fit and take the resid
    model = sm.OLS(func_data_rescaled, fit_term, missing='drop')
    model_result = model.fit()
    resid_data = model_result.resid

    # make cifti header to save residuals and coefficients
    ax1 = nib.load(standard_ts).header.get_axis(0)
    ax2 = nib.load(standard_ts).header.get_axis(1)
    header = (ax1, ax2)

    # save outputs
    resid_image = nib.cifti2.cifti2.Cifti2Image(resid_data, header)
    resid_image.to_filename(regressed_ts_filename)


def compute_dvars(dat_4compute):
    """
    Compute DVARS from a functional timeseries CIFTI file.
    https://mriqc.readthedocs.io/en/latest/iqms/bold.html
    https://nipype.readthedocs.io/en/latest/api/generated/nipype.algorithms.confounds.html#computedvars
    Input Parameters
    ----------
    dat_4compute: data needed to compute dvar and std_dvar
                  np array
    Returns
    ----------
    dvars_nstd: dvar -- inserted the first row np.nan
                np array
    dvars_std: std_dvar -- inserted the first row np.nan
               np array
    """

    dat_func_diff = np.diff(dat_4compute, axis=0)
    # DVARS (no standardization)
    dvars_nstd = np.sqrt(np.square(dat_func_diff).mean(axis=1))
    dvars_std = stats.zscore(dvars_nstd)
    dvars_nstd = np.insert(dvars_nstd, 0, np.nan)
    dvars_std = np.insert(dvars_std, 0, np.nan)

    return dvars_nstd, dvars_std


def get_custom_thres(fmriprep_confounds_tsv_file, regressed_ts_filename, thres_std_dvar=2, min_thres_fd=0.2):
    """
    Estimate customized FD threshold. Based on the regressed time series data, recompute the dvar and std_dvar,
    and find the min FD among new std_dvars that are greater than 2. By default, the minumin FD is 0.2.
    Input Parameters
    ----------
    fmriprep_confounds_tsv_file: full path to confounds_tsv_file being generated by fmriprep
                                 string
    regressed_ts_filename: full path to regressed_ts_filename (generated by regress_out_nuissance)
    Returns
    ----------
    threshold: customized FD threshold
               float
    """

    func_regressed = nib.load(regressed_ts_filename).get_fdata()

    # FD
    nuissance = pd.read_csv(fmriprep_confounds_tsv_file, sep='\t')
    fd = nuissance['framewise_displacement'].to_numpy()

    # DVARS after
    dvars_nstd_after, dvars_std_after = compute_dvars(func_regressed)
    cur_above_2 = fd[dvars_std_after > thres_std_dvar]
    min_fd = np.min(cur_above_2)
    print(f'Current min FD is {min_fd}')

    if min_thres_fd < min_fd:
        threshold = min_fd
    else:
        threshold = min_thres_fd

    return threshold


def get_sample_mask(fmriprep_confounds_tsv_file, fd_thres=0.2):
    """
    Estimate a sample mask based on the threshold. return time point * 1 arrary bool (T,F).
    True -- less than fd_thres (need to remain);
    False -- greater than fd_thres (need to remove/cencor/interpolate);
    Input Parameters
    ----------
    fmriprep_confounds_tsv_file: full path to confounds_tsv_file being generated by fmriprep
                                 string
    fd_thres: FD threshold. Default to 0.2
    Returns
    ----------
    sample_mask: sample mask (True -- less than fd_thres (need to remain))
                 time point, 1D np array
    """

    confounds_tsv = pd.read_table(fmriprep_confounds_tsv_file)["framewise_displacement"]
    confounds_tsv = confounds_tsv.fillna(0)

    return (confounds_tsv <= fd_thres).to_numpy().astype(bool)


def signal_interpolate(ts, interpolated_out, sample_mask, t_r):
    """
    !!!incomplete -- should add save to tseries data!!!
    Estimate a sample mask based on the threshold. return time point * 1 arrary bool (T,F).
    True -- less than fd_thres (need to remain);
    False -- greater than fd_thres (need to remove/cencor/interpolate);
    Input Parameters
    ----------
    ts: full path to time series data that needs to be interpolated
        string
    interpolated_out: output full path to timer series data that is already interpolated
                      string
    sample_mask: sample mask (true -- need to remain), generated by get_sample_mask
    t_r: TR in second
    Returns
    ----------
    None
    """

    func_data = nib.load(ts).get_fdata()
    interpolated_data = interpolate_volumes(func_data, sample_mask=sample_mask, t_r=t_r)

    ax1 = nib.load(ts).header.get_axis(0)
    ax2 = nib.load(ts).header.get_axis(1)
    header = (ax1, ax2)

    # save outputs
    interped_image = nib.cifti2.cifti2.Cifti2Image(interpolated_data, header)
    interped_image.to_filename(interpolated_out)


def interpolate_volumes(volumes, sample_mask, t_r):
    """
    !!!incomplete -- should add save to tseries data!!!
    Interpolate censored volumes in signals/confounds.
    copy from nilearn. https://github.com/nilearn/nilearn/blob/c33cca7a/nilearn/signal.py#L510
    """
    frame_times = np.arange(volumes.shape[0]) * t_r
    remained_vol = frame_times[sample_mask]
    remained_x = volumes[sample_mask, :]
    cubic_spline_fitter = CubicSpline(remained_vol, remained_x)
    volumes_interpolated = cubic_spline_fitter(frame_times)
    volumes[~sample_mask, :] = volumes_interpolated[~sample_mask, :]
    return volumes


def bandpass_ts(denoised_ts, filtered_out, lowpass, highpass, t_r):
    """
    Simple filter by fourier-transform
    Input Parameters
    ----------
    denoised_ts: full path to denoised time series data that needs to be interpolated
                 string
    filtered_out: output full path to timer series data that is already filtered
                  string
    lowpass and highpass: lowpass and highpass in Hz
    t_r: TR in second
    Returns
    ----------
    None
    """
    # load data and preallocate output
    func_data = nib.load(denoised_ts).get_fdata()

    sampling_rate = 1/t_r
    n_timepoints = func_data.shape[1]
    F = np.zeros((n_timepoints))

    lowidx = int(np.round(lowpass / sampling_rate * n_timepoints))
    highidx = int(np.round(highpass / sampling_rate * n_timepoints))
    F[highidx:lowidx] = 1
    F = ((F + F[::-1]) > 0).astype(int)
    filt_data = np.real(np.fft.ifftn(np.fft.fftn(func_data) * F))

    # make cifti header to save filtered data
    ax1 = nib.load(denoised_ts).header.get_axis(0)
    ax2 = nib.load(denoised_ts).header.get_axis(1)
    header = (ax1,ax2)
    # make and save image
    filt_image = nib.cifti2.cifti2.Cifti2Image(filt_data, header)
    filt_image.to_filename(filtered_out)
    return


def drop_high_motion_vols(filtered_ts, dropped_out, outlier_vols):
    """
    Last step to drop high motion volums based on the threshold
    Input Parameters
    ----------
    filtered_ts: full path to denoised time series data that needs to be dropped
                 string
    dropped_out: output full path to timer series data that is already dropped
                  string
    outlier_vols: a 1D numpy arrary (time poin) to indicate which volumes should be dropped
                  (inverse of the results from get_sample_mask)
                 time point, 1D np array
    Returns
    ----------
    None
    """
    # remove high motion volumes
    func_data = nib.load(filtered_ts).get_fdata()
    reduced_data = func_data[outlier_vols]

    # make cifti file to save data with
    ax1 = nib.load(filtered_ts).header.get_axis(0)
    ax2 = nib.load(filtered_ts).header.get_axis(1)
    ax1.size = reduced_data.shape[0]
    header = (ax1, ax2)

    # make and save image
    red_image = nib.cifti2.cifti2.Cifti2Image(reduced_data, header)
    red_image.to_filename(dropped_out)
    return
