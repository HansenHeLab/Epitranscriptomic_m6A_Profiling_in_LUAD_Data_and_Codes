#!/bin/bash

################
# for SLURM
#SBATCH -p long # only for repeat submission
#SABTCH -t 1:0:0
#SBATCH -e ./log/%x-%j.errout
#SBATCH -o ./log/%x-%j.sdout
#SBATCH --mem=30720M
################

##########################
## submission options (eg)
##########################

:<<!
#### Note the difference between the PBS and slurm for passing exteranl arguments #####
## rmdup=yes/no
## sample_pair_info.txt:   for sample with paired IP and Input,  sometime different IP could share the same Input

sbatch -J sub --export=sample_info=/path/of/sample_pair_info.txt,bam_dir=/path/of/bam/files/folder,out_dir=/output/directory,rmdup=yes_or_no,run_script=/cluster/projects/hansengroup/Yong/MeRIP-SAP/4_SLURM/step3_run_PeakCalling_Motif_h4h.sh  /cluster/projects/hansengroup/Yong/MeRIP-SAP/4_SLURM/sbatch_step3.sh


!

####################################
##  Read in Arguments for run_script
####################################

cd ${out_dir}

L=`wc -l ${sample_info} | awk '{print $1}'`	# sample size

for i in  $(seq ${L})	## {1..${L}}  doesn't work
do
	## assign sample name, IP, and Input
	Sample=`awk -v nr=${i} '{if(NR==nr){print $1}}' ${sample_info}`
	IP=`awk -v nr=${i} '{if(NR==nr){print $2}}' ${sample_info}`
	Input=`awk -v nr=${i} '{if(NR==nr){print $3}}' ${sample_info}`
	
	sbatch -J S3_${Sample} --export=bam_dir=${bam_dir},out_dir=${out_dir},rmdup=${rmdup},ip_name=${IP},input_name=${Input},name=${Sample}  ${run_script}
	echo "sbatch -J S3_${Sample} --export=bam_dir=${bam_dir},out_dir=${out_dir},rmdup=${rmdup},ip_name=${IP},input_name=${Input},name=${Sample}  ${run_script}"
	

done


