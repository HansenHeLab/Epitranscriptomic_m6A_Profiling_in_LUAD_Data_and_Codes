#!/bin/bash
#SBATCH -p long # only for repeat submission
#SABTCH -t 1:0:0
#SBATCH -e ./log/%x-%j.errout
#SBATCH -o ./log/%x-%j.sdout
#SBATCH --mem=30720M

##########################
## submission options (eg)
##########################

:<<!
#### Note the difference between the PBS and SGE for passing exteranl arguments #####
## strand=FR/RF

sbatch -J sub  --export=sample_info=/full/path/sample_seq_info.txt,out_dir=/out/bam/folder,strand_dir=RF,run_script=/full/path/to/step2_run_QC_Mapping_Pileup_Quant.sh  /full/path/to/sbatch_step2.sh 
!

####################################
##  Read in Arguments for run_script
####################################

cd ${out_dir}

L=`wc -l ${sample_info} | awk '{print $1}'`	# sample size

for i in  $(seq ${L})	## {1..${L}}  doesn't work
do
	## assign name, R1 and R2 
	seq_id=`awk -v nr=${i} '{if(NR==nr){print $2}}' ${sample_info}`
	R1=`awk -v nr=${i} '{if(NR==nr){print $3}}' ${sample_info}`
	R2=`awk -v nr=${i} '{if(NR==nr){print $4}}' ${sample_info}`

	sbatch -J S2_${seq_id}  --export=out_dir=${out_dir},strand=${strand_dir},mate1=$R1,mate2=$R2,name=${seq_id}  ${run_script}  
	echo "sbatch -J S2_${seq_id}  --export=out_dir=${out_dir},strand=${strand_dir},mate1=$R1,mate2=$R2,name=${seq_id}  ${run_script}" 

done
