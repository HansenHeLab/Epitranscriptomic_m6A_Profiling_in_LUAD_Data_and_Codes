
#!/bin/bash

# for PBS
#PBS -V
#PBS -q himem
#PBS -e ./log/
#PBS -o ./log/
#PBS -l vmem=60g

#########################

## qsub -N meta -v out_dir=/cluster/projects/hansengroup/Yong/PCa_m6a/peak/rmdup,sample_list=/cluster/projects/hansengroup/Yong/PCa_m6a/bam/sample_pair_info_1.txt,name_out=Batch1_metagene_analysis  /cluster/projects/hansengroup/Yong/PCa_m6a/script/MeRIP-SAP/step4_run_metaGene_PeakDis_Pooled_h4h


cd ${out_dir}

#################
## reference data
#################

script=/cluster/projects/hansengroup/Yong/PCa_m6a/script/MeRIP-SAP                # for R scripts
db=/cluster/projects/hansengroup/Yong/DB/hg38_ek12
ref_txdb=${db}/Guitar_TxDb_hg38.Rdata   # for guitar plot
CDS_12=${db}/Gencode_v24_CDS_12.bed
bw_dir=/cluster/projects/hansengroup/Yong/PCa_m6a/bam


###########################################
## run Guitar plot for the m6A distribution
###########################################


module load R/3.3.0

Rscript --vanilla ${script}/guitar_plot_pool.R ${sample_list} ${name_out} ${ref_txdb}


:<<!  
## collapsed as over time
#######################################################
## meta gene analysis for merged rmdup reads
#######################################################

module load deeptools/2.3.4

## 5UTR: mean/3rd Qu : 200nt      3UTR: mean/3rd Qu : ~750
## color: IP : black   , Input: gray

L=`wc -l ${sample_list} | awk '{print $1}'`

samples=""
colors=""

for i in $(seq ${L})
do
 
      samples="$samples ${bw_dir}/`awk -v nr=${i} '{if(NR==nr){print $2}}' ${sample_list}`_minus_pileup_normed.bw  ${bw_dir}/`awk -v nr=${i} '{if(NR==nr){print $2}}' ${sample_list}`_plus_pileup_normed.bw ${bw_dir}/`awk -v nr=${i} '{if(NR==nr){print $3}}' ${sample_list}`_minus_pileup_normed.bw  ${bw_dir}/`awk -v nr=${i} '{if(NR==nr){print $3}}' ${sample_list}`_plus_pileup_normed.bw"

      colors="$colors black black gray gray"
done

computeMatrix scale-regions  --metagene  --skipZeros -b 1000 -m 4000  -a 2000  -R ${CDS_12}  -S ${samples}  -o ${name_out}_pooled_CM.gz
plotProfile -m ${name_out}_pooled_CM.gz  --perGroup  --startLabel "Start codon"  --endLabel  "Stop codon"  --legendLocation none --colors ${colors}  -out  ${name_out}_pooled_metagene.pdf

module unload deeptools/2.3.4    
!



