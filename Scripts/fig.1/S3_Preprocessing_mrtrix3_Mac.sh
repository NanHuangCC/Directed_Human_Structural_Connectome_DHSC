#!/bin/bash

# conbvert niftl to native format MIF
subject="$1"
mkdir -p /data3/nh_data/20250204_DirectionalConnectome/dataset/data/macaque/${subject}/connectivity/mrtrix
workdir="/data3/nh_data/20250204_DirectionalConnectome/dataset/data/macaque/${subject}/connectivity/mrtrix"
anatdir="/data3/nh_data/20250204_DirectionalConnectome/dataset/data/macaque/${subject}/anat"
dwidir="/data3/nh_data/20250204_DirectionalConnectome/dataset/data/macaque/${subject}/diffusion"

cd  ${workdir}
bvecs=${dwidir}/bvecs
bvals=${dwidir}/bvals
dwi=${dwidir}/data.nii.gz

# generate .mif file
mrconvert -fslgrad ${bvecs} ${bvals} ${dwi} ${subject}.mif 

# creat mask for preproce
dwi2mask ${subject}.mif - | maskfilter - dilate preproce_mask.mif -npass 2

# denoising of diffusion MRI data
dwidenoise ${subject}.mif denoise.mif -noise noiselevel.mif -mask preproce_mask.mif -force

# Gibbs Ring artifeact correction
mrdegibbs denoise.mif degibbs.mif 

# motion and distortion correction
# creat b0 AP-PA pair 
# dwiextract raw1000AP.mif b0_AP.mif -bzero 
# dwiextract raw1000.mif - bzero | mrconvert - -coord 3 0 b0_PA.mif
# mrcat b0_AP.mif b0_PA.mif B0_pair.mif
# dwifslpreproc degibbs.mif geomorr.mif -pe_dir PA -rpe_pair -se_epi b0_pair.mif

# bias filed correction
# dwibiascorrect ants degibbs.mif biascorr.mif # if do geometry correction use geomorr.mif

# use recon-all segment all brain regions
# get get stop masks from reslut of freesufer
5ttgen freesurfer ${anatdir}/aseg.nii.gz 5tt.mif
5tt2gmwmi 5tt.mif 5tt_gmwm.mif

dwi2mask degibbs.mif - | maskfilter - dilate dwi_mask.mif
dwi2tensor -mask dwi_mask.mif degibbs.mif dt.mif
# caculate FA
tensor2metric -fa fa.mif -mask dwi_mask.mif dt.mif
# creat wm ROIs
mrthreshold -abs 0.4 fa.mif - | mrcalc - dwi_mask.mif -mult wm_mask.mif -force
# Estimate response function for CSD
dwi2response dhollander degibbs.mif response_wm.txt response_gm.txt response_csf.txt
### constrained spherical deconvolution (CSD) #####
# Perform CSD
dwi2fod msmt_csd degibbs.mif response_wm.txt wm.mif response_gm.txt gm.mif response_csf.txt csf.mif -mask dwi_mask.mif
# Perform tractography
tckgen -algorithm iFOD2 -act 5tt.mif -backtrack -crop_at_gmwmi -select 10M -angle 45 -cutoff 0.05 -minlength 2 -maxlength 300 -step 0.45 wm.mif tracks10M.tck -seed_image dwi_mask.mif -force
# gnerate connectome
# connectome ref D99 parcellation
tck2connectome -symmetric -zero_diagonal tracks10M.tck ${anatdir}/D99_atlas_in_T2.nii.gz connectome.csv -force 
tck2connectome tracks10M.tck ${anatdir}/D99_atlas_in_T2.nii.gz distances.csv -scale_length -stat_edge mean -force