#!/bin/bash
# test xfms for human subject
subjuct=$1
# basic path for atlas and subject
MNI152="/data3/nh_data/20250204_DirectionalConnectome/dataset/atlas/human/mni_icbm152_t1_tal_nlin_sym_09a_converted.nii.gz"
subject_path="/data3/nh_data/20250204_DirectionalConnectome/dataset/data/human/${subjuct}"
xfms_path="/data3/nh_data/20250204_DirectionalConnectome/dataset/data/human/${subjuct}/xfms"
script_path="/data3/nh_data/20250204_DirectionalConnectome/scripts/Main/1_Brian_mapping"

xtract -native -bpx ${subject_path}/Diffusion.bedpostX -out ${subject_path}/xtract -stdwarp ${xfms_path}/standard2acpc_dc.nii.gz ${xfms_path}/acpc_dc2standard.nii.gz -species HUMAN -gpu
xtract_blueprint -bpx ${subject_path}/Diffusion.bedpostX -xtract ${subject_path}/xtract -out ${subject_path}/xtract/blueprint -seeds ${subject_path}/surf/${subjuct}.L.white_MSMAll.32k_fs_LR.surf.gii,${subject_path}/surf/${subjuct}.R.white_MSMAll.32k_fs_LR.surf.gii -target ${subject_path}/Diffusion/nodif_brain_mask.nii.gz -native -gpu -nsamples 1000 
bash $script_path/S2_xtract_blueprint_subcortical.sh ${subjuct} 
bash $script_path/S3_Preprocessing_mrtrix3_Human.sh ${subjuct} 