#!/bin/bash
subj="$1"
data_path="/data3/nh_data/20250204_DirectionalConnectome/dataset/data/human/${subj}"
subcortical_regions="/data3/nh_data/20250204_DirectionalConnectome/dataset/data/human/HCPex_subject/${subj}/HCPex_ROIs_sym/subcortex"
blueprint_path=${data_path}/xtract/blueprint_subcortex
mkdir ${blueprint_path}

find $subcortical_regions -maxdepth 1 -type f -printf "%f\n" | while read filename; do
    name=$(echo "$filename" | awk -F'.' '{print $1}')
    mkdir -p ${blueprint_path}/${name}
    xtract_blueprint -bpx ${data_path}/Diffusion.bedpostX -xtract ${data_path}/xtract -out ${blueprint_path}/${name} -seeds ${subcortical_regions}/${filename} -target ${data_path}/Diffusion/nodif_brain_mask.nii.gz -native -nsamples 1000 -gpu
done
