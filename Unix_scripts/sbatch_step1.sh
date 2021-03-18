#!/bin/baqsh

################
# for SLURM
#SBATCH -p long	# only for repeat submission
#SABTCH -t 1:0:0
#SBATCH -e ./log/%x-%j.errout
#SBATCH -o ./log/%x-%j.sdout
#SBATCH --mem=30720M 
################

##########################
## submission options (eg)
##########################

:<<!
#########################
## sbatch submition e.g.:

sbatch -J sub  --export=sample_info=/full/path/sample_seq_info.txt,out_dir=/Directory/for/output,run_script=/cluster/projects/hansengroup/Yong/MeRIP-SAP/4_SLURM/step1_QC_test_h4h.sh  /cluster/projects/hansengroup/Yong/MeRIP-SAP/4_SLURM/sbatch_step1_QC_test.sh 
!

####################################
##  Read in Arguments for run_script
####################################

cd ${out_dir}

L=`wc -l ${sample_info} | awk '{print $1}'`	# sample size

for i in 1 ${L}    ##  $(seq ${L})	## test the first and last sample
do
	## assign name, R1 and R2 
	seq_id=`awk -v nr=${i} '{if(NR==nr){print $2}}' ${sample_info}`_test  ##addsuffix 
	R1=`awk -v nr=${i} '{if(NR==nr){print $3}}' ${sample_info}`
	R2=`awk -v nr=${i} '{if(NR==nr){print $4}}' ${sample_info}`

	## qsub the run_script: out_dir, mate1, mate2, name are predefiend in run_script

	sbatch -J ${seq_id} --export=out_dir=${out_dir},mate1=$R1,mate2=$R2,name=${seq_id} ${run_script}
	echo "sbatch -J ${seq_id} --export=out_dir=${out_dir},mate1=$R1,mate2=$R2,name=${seq_id}  ${run_script}"

done


