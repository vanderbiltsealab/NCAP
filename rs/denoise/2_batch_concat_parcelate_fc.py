from glob import glob
import os
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

for root, dirs, files in os.walk(target_dir):
    for directory in dirs:
        if (directory == 'post_fmriprep'):
            print(root)
            print(directory)

            cur_func_root = os.path.join(root, directory)

            # glob fully denosied data
            cur_denoised = glob(f'{cur_func_root}/*_dropped.32k_fs_LR.dtseries.nii')

            if len(cur_denoised) == 0:
                print("=========ERROR========= if len(cur_denoised) == 0:")
                continue

            # concatenate
            cur_base_name = os.path.basename(cur_denoised[0])
            cur_out_name = cur_base_name[:cur_base_name.find("-restPA") + len("-restPA")]
            concat_out = os.path.join(cur_func_root, f'{cur_out_name}_concatenated.dtseries.nii')
            den_sl.concatenate_dtseries(cur_denoised, concat_out)

            # parcelate
            parcelate_out = concat_out.replace('.dtseries', '_gordon.32k_fs_LR.ptseries')
            den_sl.apply_parcels(concat_out, parcelate_out, atlas_dlabel)

            # FC
            fc_out = parcelate_out.replace('.ptseries.nii', '.pconn.nii')
            den_sl.fc_estimate(parcelate_out, fc_out)

print('Done')
