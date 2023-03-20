from glob import glob
import nibabel as nib
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import re
from scipy.stats import pearsonr
import seaborn as sns
import statistics
import sys

sys.path.append("/Users/yanbinniu/Projects/NCAP/scripts/rs/utils")
import denoise_sealab as den_sl


# change to current directory
os.chdir('/Users/yanbinniu/Projects/NCAP/scripts/rs/custom_thres')
"""
# 1. Compute the correlation between FD and FC.
# 2. Visualize the correlation matrix.
"""


def get_mean_FD_BeforeAfter_censor (fmriprep_confounds_tsv_file, idx_censor):
    nuissance = pd.read_csv(fmriprep_confounds_tsv_file, sep='\t')
    fd = nuissance['framewise_displacement'].to_numpy()
    return np.nanmean(fd), np.nanmean(fd[~idx_censor])


# target directory
target_dir = '/Volumes/ppm/PPM/NCAP_bids/Data/derivatives'
cus_thres = '/Volumes/ppm/PPM/NCAP_bids/Code/fmriprep/custom_fd_th_info.csv'
df_cus_thres = pd.read_csv(cus_thres)

# Step 1 -- compute the following lists
# mean FD (before censor) for all subject
mFD_beforeCensor_list = []

# mean FD (after censor) for all subject
mFD_afterCensor_list = []

# number of frames left for all subject
frame_size_list = []

# FC for all subject, 407*407*30
fc_all = np.empty((407, 407, 30))

# subject list
subj_list = []

for root, dirs, files in os.walk(target_dir):
    for directory in dirs:
        if directory == 'post_fmriprep':
            print(root)
            print(directory)

            cur_func_root = os.path.join(root, directory)

            # compute FD
            cur_conn = glob(f'{cur_func_root}/*gordon.32k_fs_LR.pconn.nii')

            # validity checking
            if 1 != len(cur_conn):
                print('ERROR ======= if 1 != len(cur_conn)')
                continue

            cur_idx = len(mFD_afterCensor_list)
            # save fc
            temp_fc = nib.load(cur_conn[0]).get_fdata()
            fc_all[:, :, cur_idx] = temp_fc
            cur_subj = root.split('/')[-3]
            subj_list.append(cur_subj.split('-')[-1])


            cur_mFD_before_Censor = []
            cur_mFD_after_Censor = []
            cur_frame_size = []
            cur_func_par = os.path.dirname(cur_func_root)
            cur_confounds = glob(f'{cur_func_par}/*desc-confounds_timeseries.tsv')

            for tsv in cur_confounds:
                # get the run_x for session column
                match = re.search(r'run-(\d)', tsv)
                if match:
                    one_char = match.group(1)
                else:
                    one_char = "1"

                tmp_scan = f'{cur_subj}_run-{one_char}'
                # customized fd threshold
                # use .iloc to select the row(s) where subj equals tmp_scan, and get the value of fd_thres
                tmp_fd_th = df_cus_thres.loc[df_cus_thres['session'] == tmp_scan, 'fd_thres'].values[0]

                # get sample mask
                cur_sample_mask = den_sl.get_sample_mask(tsv, tmp_fd_th)

                mFD1, mFD2 = get_mean_FD_BeforeAfter_censor(tsv, ~cur_sample_mask)
                cur_mFD_before_Censor.append(mFD1)
                cur_mFD_after_Censor.append(mFD2)
                cur_frame_size.append(cur_sample_mask.sum())

            # compute for this scan
            mFD_beforeCensor_list.append(statistics.mean(cur_mFD_before_Censor))
            mFD_afterCensor_list.append(statistics.mean(cur_mFD_after_Censor))
            frame_size_list.append(sum(cur_frame_size))

# save data into numpy files
tmp_arr1 = np.array(mFD_beforeCensor_list)
tmp_arr2 = np.array(mFD_afterCensor_list)
tmp_arr3 = np.array(frame_size_list)
subj_arr = np.array(subj_list)
np.save("mFD_beforeCensor_list.npy", tmp_arr1)
np.save("mFD_afterCensor_list.npy", tmp_arr2)
np.save("frame_size_list.npy", tmp_arr3)
np.save("subj_list.npy", subj_arr)
np.save("fc_all.npy", fc_all)


# tabulate frame size info
frame_size_list = np.load("frame_size_list.npy")
subject_list = np.load("subj_list.npy")
mFD_beforeCensor_list = np.load("mFD_beforeCensor_list.npy")
mFD_afterCensor_list = np.load("mFD_afterCensor_list.npy")
time_length_info = pd.DataFrame({"subject": subject_list,
                                 "frame_size": frame_size_list,
                                 "mFD_beforeCensor": mFD_beforeCensor_list,
                                 "mFD_afterCensor": mFD_afterCensor_list})
# tr = 2.0
time_length_info['time_length'] = time_length_info['frame_size'] * 2.0/60
time_length_info = time_length_info.sort_values(by=['subject'])
time_length_info.to_csv('mFD_frameSize_info.csv', index=False)


# tabulate mfd
df_mFD_frameSize = pd.DataFrame({"subject": subject_list,
                                 "frame_size": frame_size_list})


# compute correlation matrix
mFD_beforeCensor_list = np.load("mFD_beforeCensor_list.npy").tolist()
mFD_afterCensor_list = np.load("mFD_afterCensor_list.npy").tolist()
frame_size_list = np.load("frame_size_list.npy")
fc_all = np.load("fc_all.npy")
loop_size = fc_all.shape[0]
matrix_fd_beforeCensor = np.zeros((loop_size, loop_size))
matrix_fd_afterCensor = np.zeros((loop_size, loop_size))
matrix_fd_frameSize = np.zeros((loop_size, loop_size))

for i in range(0, loop_size):
    for j in range(0, loop_size):
        if i == j:
            matrix_fd_beforeCensor[i, j] = np.nan
            matrix_fd_afterCensor[i, j] = np.nan
            matrix_fd_frameSize[i, j] = np.nan
            continue

        tmp_fc_list = fc_all[i, j, :].tolist()
        if pd.isna(tmp_fc_list).sum() == 30:
            matrix_fd_beforeCensor[i, j] = np.nan
            matrix_fd_afterCensor[i, j] = np.nan
            matrix_fd_frameSize[i, j] = np.nan
        else:
            matrix_fd_beforeCensor[i, j] = pearsonr(tmp_fc_list, mFD_beforeCensor_list)[0]
            matrix_fd_afterCensor[i, j] = pearsonr(tmp_fc_list, mFD_afterCensor_list)[0]
            matrix_fd_frameSize[i, j] = pearsonr(tmp_fc_list, frame_size_list)[0]


# save
np.savetxt('matrix_fd_beforeCensor.txt', matrix_fd_beforeCensor)
np.savetxt('matrix_fd_afterCensor.txt', matrix_fd_afterCensor)
np.savetxt('matrix_fd_frameSize.txt', matrix_fd_frameSize)

# read in matrix
matrix_fd_beforeCensor = np.loadtxt('matrix_fd_beforeCensor.txt')
matrix_fd_afterCensor = np.loadtxt('matrix_fd_afterCensor.txt')
matrix_fd_frameSize = np.loadtxt('matrix_fd_frameSize.txt')

# plotting
# get network labels
parcel_labels = np.load("/Users/yanbinniu/Projects/NCAP/scripts/rs/utils/parcel_labels.npy")
network_labels = np.load("/Users/yanbinniu/Projects/NCAP/scripts/rs/utils/network_labels.npy")
network_names = np.load("/Users/yanbinniu/Projects/NCAP/scripts/rs/utils/network_names.npy")

network_labels_cleaned = []
for idx, ele in enumerate(network_labels):
    if 0 == idx:
        network_labels_cleaned.append(ele)
    else:
        if ele == network_labels[idx-1]:
            network_labels_cleaned.append('')
        else:
            network_labels_cleaned.append(ele)


# plotting matrix_fd_beforeCensor
sns.set_theme(style="white")
# Set up the matplotlib figure
plt.rcParams["figure.figsize"] = [10, 14]
plt.rcParams["figure.autolayout"] = True
fig, (ax1, ax2) = plt.subplots(nrows=2)
# Generate a custom diverging colormap
cmap = sns.diverging_palette(230, 20, as_cmap=True)
# Draw the heatmap with the mask and correct aspect ratio
sns.heatmap(matrix_fd_beforeCensor, cmap=cmap, square=True, linewidths=0, ax=ax1,
            vmin=-.8, vmax=.8, yticklabels=network_labels_cleaned, xticklabels=False)
ax2.hist(matrix_fd_beforeCensor.flatten(), bins=40, color='#607c8e')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
fig.savefig("matrix_fd_beforeCensor.png")
plt.close(fig)

# plotting matrix_fd_afterCensor
sns.set_theme(style="white")
# Set up the matplotlib figure
plt.rcParams["figure.figsize"] = [10, 14]
plt.rcParams["figure.autolayout"] = True
fig, (ax1, ax2) = plt.subplots(nrows=2)
# Generate a custom diverging colormap
cmap = sns.diverging_palette(230, 20, as_cmap=True)
# Draw the heatmap with the mask and correct aspect ratio
sns.heatmap(matrix_fd_afterCensor, cmap=cmap, square=True, linewidths=0, ax=ax1,
            vmin=-.8, vmax=.8, yticklabels=network_labels_cleaned, xticklabels=False)
ax2.hist(matrix_fd_afterCensor.flatten(), bins=40, color='#607c8e')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
fig.savefig("matrix_fd_afterCensor.png")
plt.close(fig)

# plotting matrix_fd_afterCensor
sns.set_theme(style="white")
# Set up the matplotlib figure
plt.rcParams["figure.figsize"] = [10, 14]
plt.rcParams["figure.autolayout"] = True
fig, (ax1, ax2) = plt.subplots(nrows=2)
# Generate a custom diverging colormap
cmap = sns.diverging_palette(230, 20, as_cmap=True)
# Draw the heatmap with the mask and correct aspect ratio
sns.heatmap(matrix_fd_frameSize, cmap=cmap, square=True, linewidths=0, ax=ax1,
            vmin=-.8, vmax=.8, yticklabels=network_labels_cleaned, xticklabels=False)
ax2.hist(matrix_fd_frameSize.flatten(), bins=40, color='#607c8e')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
fig.savefig("matrix_fd_frameSize.png")
plt.close(fig)
