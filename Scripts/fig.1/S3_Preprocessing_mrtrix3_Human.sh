#!/bin/bash
# apply surface based cortex prarcellation to niftii file
sub=$1
Human_dir="/data3/nh_data/20250204_DirectionalConnectome/dataset/data/human"
HCPex_folder="${Human_dir}/HCPex_subject/${sub}"
HCP_results_folder="${Human_dir}/${sub}"
Diffusion_folder="${HCP_results_folder}/Diffusion"
GlasserFolder="/data3/nh_data/20250204_DirectionalConnectome/dataset/atlas/human/Glasser_et_al_2016_HCP_MMP1.0_RVVG_2/HCP_PhaseTwo/Glasser_et_al_2016_HCP_MMP1.0_StudyData"


# map surface labels from Glasser parcellation to volume
for Hemisphere in R L; do   
    Glasser_label_gii="${GlasserFolder}/${Hemisphere}.label.gii" # used "wb_command -cifti-separate" to seperate the the CIFTI label into hemisphere wise gii label files
    white_fsLR32="${HCP_results_folder}/surf/${sub}.${Hemisphere}.white_MSMAll.32k_fs_LR.surf.gii"
    pial_fsLR32="${HCP_results_folder}/surf/${sub}.${Hemisphere}.pial_MSMAll.32k_fs_LR.surf.gii"
    T1_image="${HCP_results_folder}/anat/T1w_acpc_dc_restore_brain.nii.gz"
    output_volume_label="${HCP_results_folder}/anat/${Hemisphere}.HCP_MMP1.cortical_volume_labels.nii.gz"

    echo "Label to volume mapping hemisphere:  $Hemisphere"
    wb_command -label-to-volume-mapping $Glasser_label_gii $white_fsLR32 $T1_image $output_volume_label -ribbon-constrained $white_fsLR32 $pial_fsLR32
done

# get surfer based 
fslmaths ${HCPex_folder}/HCPex_in_T1w.nii.gz -thr 361 ${HCP_results_folder}/anat/subcortical_regions.nii.gz

fslmaths ${HCP_results_folder}/anat/subcortical_regions.nii.gz -binv ${HCP_results_folder}/anat/mask_zero.nii.gz
fslmaths ${HCP_results_folder}/anat/L.HCP_MMP1.cortical_volume_labels.nii.gz -mul ${HCP_results_folder}/anat/mask_zero.nii.gz ${HCP_results_folder}/anat/cortical_to_add.nii.gz
fslmaths ${HCP_results_folder}/anat/subcortical_regions.nii.gz -add ${HCP_results_folder}/anat/cortical_to_add.nii.gz ${HCP_results_folder}/anat/HCPex_surf_sub.nii.gz

fslmaths ${HCP_results_folder}/anat/HCPex_surf_sub.nii.gz -binv ${HCP_results_folder}/anat/mask_zero.nii.gz
fslmaths ${HCP_results_folder}/anat/R.HCP_MMP1.cortical_volume_labels.nii.gz -mul ${HCP_results_folder}/anat/mask_zero.nii.gz ${HCP_results_folder}/anat/cortical_to_add.nii.gz
fslmaths ${HCP_results_folder}/anat/HCPex_surf_sub.nii.gz -add ${HCP_results_folder}/anat/cortical_to_add.nii.gz ${HCP_results_folder}/anat/HCPex_surf_sub.nii.gz
rm ${HCP_results_folder}/anat/mask_zero.nii.gz ${HCP_results_folder}/anat/cortical_to_add.nii.gz
# conbvert niftl to native format MIF

mkdir -p ${HCP_results_folder}/connectivity/mrtrix
workdir="${HCP_results_folder}/connectivity/mrtrix"
anatdir="${HCP_results_folder}/anat"
dwidir="${HCP_results_folder}/Diffusion"

bvecs=${dwidir}/bvecs
bvals=${dwidir}/bvals
dwi=${dwidir}/data.nii.gz

# generate .mif file
mrconvert -fslgrad ${bvecs} ${bvals} ${dwi} ${workdir}/${sub}.mif -force
# creat mask for preproce
dwi2mask ${workdir}/${sub}.mif - | maskfilter - dilate ${workdir}/preproce_mask.mif -npass 2 -force
# denoising of diffusion MRI data
dwidenoise ${workdir}/${sub}.mif ${workdir}/denoise.mif -noise ${workdir}/noiselevel.mif -mask ${workdir}/preproce_mask.mif -force
# Gibbs Ring artifeact correction
mrdegibbs ${workdir}/denoise.mif ${workdir}/degibbs.mif -force
rm ${workdir}/denoise.mif ${workdir}/noiselevel.mif ${workdir}/preproce_mask.mif 
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
5ttgen freesurfer ${anatdir}/aparc+aseg.nii.gz  ${workdir}/5tt.mif -force
5tt2gmwmi ${workdir}/5tt.mif ${workdir}/5tt_gmwm.mif -force

dwi2mask ${workdir}/degibbs.mif - | maskfilter - dilate ${workdir}/dwi_mask.mif -force
dwi2tensor -mask ${workdir}/dwi_mask.mif ${workdir}/degibbs.mif ${workdir}/dt.mif -force
# caculate FA
tensor2metric -fa ${workdir}/fa.mif -mask ${workdir}/dwi_mask.mif ${workdir}/dt.mif -force
# creat wm ROIs
mrthreshold -abs 0.4 ${workdir}/fa.mif - | mrcalc - ${workdir}/dwi_mask.mif -mult ${workdir}/wm_mask.mif -force
# Estimate response function for CSD
dwi2response dhollander ${workdir}/degibbs.mif ${workdir}/response_wm.txt ${workdir}/response_gm.txt ${workdir}/response_csf.txt -force
### constrained spherical deconvolution (CSD) #####
# Perform CSD
dwi2fod msmt_csd ${workdir}/degibbs.mif ${workdir}/response_wm.txt ${workdir}/wm.mif ${workdir}/response_gm.txt ${workdir}/gm.mif ${workdir}/response_csf.txt ${workdir}/csf.mif -mask ${workdir}/dwi_mask.mif -force
# Perform tractography
tckgen -algorithm iFOD2 -act ${workdir}/5tt.mif -backtrack -crop_at_gmwmi -select 10M -angle 45 -cutoff 0.05 -minlength 5 -maxlength 300 -step 1.0 ${workdir}/wm.mif ${workdir}/tracks10M.tck -seed_image ${workdir}/dwi_mask.mif -force 
# gnerate connectome
# connectome ref D99 parcellation
tck2connectome -symmetric -zero_diagonal ${workdir}/tracks10M.tck ${anatdir}/HCPex_surf_sub.nii.gz ${workdir}/connectome.csv -force 
tck2connectome ${workdir}/tracks10M.tck ${anatdir}/HCPex_surf_sub.nii.gz ${workdir}/distances.csv -scale_length -stat_edge mean -force 

