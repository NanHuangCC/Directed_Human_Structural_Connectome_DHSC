#!/bin/bash
subj="$1"
subjdir="/data3/nh_data/20250204_DirectionalConnectome/dataset/data/macaque/${subj}"
atlasdir="/data3/nh_data/20250204_DirectionalConnectome/dataset/atlas/macaque/D99_v2.0_dist/D99_v2.0_dist"
subcortical_regions="${subjdir}/anat/D99_LR_ATLAS_regions/subcortex_right"
script_path="/data3/nh_data/20250204_DirectionalConnectome/scripts/Main/1_Brian_mapping"

cd ${subjdir}
# cp ./dwi/${subj}_dwi.bval ./diffusion/bvals
# cp ./dwi/${subj}_dwi.bvec ./diffusion/bvecs
# get bedpostx
# bedpostx_gpu ./diffusion -n 3 -model 2 -NJOBS 8 
# perform xtract
xtract -native -bpx ./diffusion.bedpostX -out ./xtract -stdwarp ./xfms/F99toB0_fsl_warp.nii.gz ./xfms/B0toF99_fsl_warp.nii.gz -species MACAQUE -gpu 
xtract_blueprint -bpx ${subjdir}/diffusion.bedpostX -xtract ${subjdir}/xtract -out ${subjdir}/xtract/blueprint -seeds ${subjdir}/surf/${subj}/surf/lh.white.surf.gii,${subjdir}/surf/${subj}/surf/rh.white.surf.gii -target ${subjdir}/anat/B0_mask.nii.gz -native -nsamples 1000 -gpu
bash $script_path/S2_xtract_blueprint_subcortical.sh ${subj} 