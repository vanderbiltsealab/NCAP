from glob import glob
import os
import pandas as pd
import re
import sys

sys.path.append("/Users/yanbinniu/Projects/NCAP/scripts/rs/utils")
import denoise_sealab as den_sl

"""
# Estimate custom FD threshold for xcp_d pipeline
"""

# change to current directory
os.chdir('/Users/yanbinniu/Projects/NCAP/scripts/rs/custom_thres')

# first estimate a csv for custom FD threshold with two columns (one for scan, i.e., sub-10001scanx
# and one for custom FD threshold)

# target directory
target_dir = '/Volumes/ppm/PPM/NCAP_bids/Data/derivatives'
atlas_dlabel = '/Users/yanbinniu/Projects/NCAP/scripts/rs/utils/Gordon333_SeitzmanSubcortical.32k_fs_LR.dlabel.nii'
custom_fd_th_info = pd.DataFrame(columns=['session', 'fd_thres'])
t_r = 2.0
lowpass = 0.1  # in Hz
highpass = 0.008  # in Hz

for root, dirs, files in os.walk(target_dir):
    for directory in dirs:
        if (directory == 'func') and ('xcp_d' not in root):
            print(root)
            print(directory)

            cur_func_root = os.path.join(root, directory)
            cur_working_root = os.path.join(cur_func_root, 'post_fmriprep')
            cur_scan = root.split('/')[-2]

            # create tractseg directory under current subject
            if not os.path.exists(cur_working_root):
                os.makedirs(cur_working_root)

            # glob confounds and data preprocessed by fmriprep
            cur_confounds = glob(f'{cur_func_root}/*desc-confounds_timeseries.tsv')
            cur_cifti = glob(f'{cur_func_root}/*fsLR_den-91k_bold.dtseries.nii')

            if len(cur_confounds) != len(cur_cifti):
                print(f'======ERROR======if len(cur_confounds) != len(cur_cifti):')
                continue

            for tsv in cur_confounds:
                tmp_cifti = tsv.replace('desc-confounds_timeseries.tsv', 'space-fsLR_den-91k_bold.dtseries.nii')
                tsv_base = os.path.basename(tsv)
                tmp_cifti_base = os.path.basename(tmp_cifti)

                # resample
                resampled_out = os.path.join(cur_working_root,
                                            tmp_cifti_base.replace('.dtseries', '.32k_fs_LR.dtseries'))
                den_sl.wb_resample_32k(tmp_cifti, resampled_out, atlas_dlabel)

                # rescale
                rescaled_out = resampled_out.replace('.32k_fs_LR', '_rescaled.32k_fs_LR')
                den_sl.rescale_signal(resampled_out, rescaled_out)

                # generate 32P regressor txt file (WITHOUT threshold)
                regress_32p_txt_out = os.path.join(cur_working_root,
                                                   tsv_base.replace('.tsv', '_32p.txt'))
                den_sl.make_noise_matrix(tsv, regress_32p_txt_out)

                # generated regressed ts data (no despike)
                regressed_32p_out = rescaled_out.replace('.32k_fs_LR', '_regressed32P.32k_fs_LR')
                den_sl.regress_out_nuissance(regress_32p_txt_out, rescaled_out, regressed_32p_out)

                # estimate customized fd threshold
                tmp_fd_th = den_sl.get_custom_thres(tsv, regressed_32p_out)

                # get the run_x for session column
                match = re.search(r'run-(\d)', tmp_cifti)
                if match:
                    one_char = match.group(1)
                else:
                    one_char = "1"

                # organize custom fd threshold dataframe
                custom_fd_th_info = custom_fd_th_info.append({'session': f'{cur_scan}_run-{one_char}',
                                                              'fd_thres': tmp_fd_th}, ignore_index=True)

                # generate 32P regressor txt file (WITH threshold)
                regress_32p_despike_txt_out = os.path.join(cur_working_root,
                                                   tsv_base.replace('.tsv', '_32p_despike.txt'))
                print(f'{cur_scan}_run-{one_char} custom thres is {tmp_fd_th}')
                den_sl.make_noise_matrix(tsv, regress_32p_despike_txt_out, 50, tmp_fd_th)

                # generated regressed ts data (despiked)
                regressed_32p_despike_out = rescaled_out.replace('.32k_fs_LR',
                                                                 '_regressed32P_despiked.32k_fs_LR')
                den_sl.regress_out_nuissance(regress_32p_despike_txt_out,
                                             rescaled_out, regressed_32p_despike_out)

                # get sample mask
                cur_sample_mask = den_sl.get_sample_mask(tsv, tmp_fd_th)

                # # interpolate
                # regressed_32p_despike_intered_out = resample_rescale_out.replace('.32k_fs_LR',
                #                                                                  '_interpolated.32k_fs_LR')
                # den_sl.signal_interpolate(regressed_32p_despike_out, regressed_32p_despike_intered_out,
                #                           cur_sample_mask, t_r)


                # filtering
                filtered_out = regressed_32p_despike_out.replace('.32k_fs_LR', '_filtered.32k_fs_LR')
                den_sl.bandpass_ts(regressed_32p_despike_out, filtered_out,
                                   lowpass, highpass, t_r)

                # drop
                # note to take the inverse of sample mask
                dropped_out = filtered_out.replace('.32k_fs_LR', '_dropped.32k_fs_LR')
                den_sl.drop_high_motion_vols(filtered_out, dropped_out, ~cur_sample_mask)

custom_fd_th_info.to_csv('custom_fd_th_info.csv', index=False)
print('Done')
