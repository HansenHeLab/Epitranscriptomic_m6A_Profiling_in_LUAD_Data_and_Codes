#!/bin/bash
#SBATCH -p himem
#SABTCH -t 3-0:0:0
#SBATCH -e ./log/%x-%j.errout
#SBATCH -o ./log/%x-%j.sdout
#SBATCH --mem=30000M

#### array job submision by SLURM
## eg 
## sbatch -J peak_m6a --array=2-126:2  --export=peak_bed=/path/to/normal_tumor_all_merge_labeled.bed,out_dir=/dir/of/output  /path/to/step4_merged_peak_m6a_level.sh

cd ${out_dir}

idx_1=${SLURM_ARRAY_TASK_ID}-2; 
idx_2=${SLURM_ARRAY_TASK_ID}-1; 

file_name=(normal01_input normal01_ip normal02_input normal02_ip normal03_input normal03_ip normal04_input normal04_ip normal05_input normal05_ip normal06_input normal06_ip normal07_input normal07_ip normal08_input normal08_ip normal09_input normal09_ip normal10_input normal10_ip tumor10_input tumor10_ip tumor11_input tumor11_ip tumor12_input tumor12_ip tumor13_input tumor13_ip tumor14_input tumor14_ip tumor15_input tumor15_ip tumor16_input tumor16_ip tumor17_input tumor17_ip tumor18_input tumor18_ip tumor19_input tumor19_ip tumor2_input tumor2_ip tumor20_input tumor20_ip tumor21_input tumor21_ip tumor22_input tumor22_ip tumor23_input tumor23_ip tumor24_input tumor24_ip tumor25_input tumor25_ip tumor26_input tumor26_ip tumor27_input tumor27_ip tumor28_input tumor28_ip tumor29_input tumor29_ip tumor3_input tumor3_ip tumor30_input tumor30_ip tumor31_input tumor31_ip tumor32_input tumor32_ip tumor33_input tumor33_ip tumor34_input tumor34_ip tumor35_input tumor35_ip tumor36_input tumor36_ip tumor37_input tumor37_ip tumor38_input tumor38_ip tumor39_input tumor39_ip tumor4_input tumor4_ip tumor40_input tumor40_ip tumor41_input tumor41_ip tumor42_input tumor42_ip tumor43_input tumor43_ip tumor44_input tumor44_ip tumor45_input tumor45_ip tumor47_input tumor47_ip tumor48_input tumor48_ip tumor49_input tumor49_ip tumor5_input tumor5_ip tumor50_input tumor50_ip tumor51_input tumor51_ip tumor52_input tumor52_ip tumor53_input tumor53_ip tumor54_input tumor54_ip tumor55_input tumor55_ip tumor6_input tumor6_ip tumor7_input tumor7_ip tumor8_input tumor8_ip tumor9_input tumor9_ip)

base_dir=/cluster/projects/hansengroup/Yong/LUAD_m6a/single_end_results/bw/

name_ip=${file_name[${idx_2}]};     
name_input=${file_name[${idx_1}]};  
                                    
ip_p=${base_dir}${name_ip}_plus_pileup_normed.bw                                                                                             
ip_m=${base_dir}${name_ip}_minus_pileup_normed.bw                                                                                            
input_p=${base_dir}${name_input}_plus_pileup_normed.bw                                                                                       
input_m=${base_dir}${name_input}_minus_pileup_normed.bw                                                                                      
                                                                                                                                             
                                                                                                                                             
#########################                                                                                                                                                                                                                                                                
module load R/3.3.0                                                                                                                          
module load ucsctools/378      #for slurm                                                                                                                    
module load bedtools/2.27.1

# file preparation                                                                                                                           
bigWigToBedGraph $ip_p ${name_ip}_plus_pileup.bdg                                                                                            
bigWigToBedGraph $ip_m ${name_ip}_minus_pileup.bdg                                                                                           
bigWigToBedGraph $input_p ${name_input}_plus_pileup.bdg                                                                                      
bigWigToBedGraph $input_m ${name_input}_minus_pileup.bdg                                                                                     
                                                                                                                                             
awk '{print $0 "\t" "*" "\t" "+" }'  ${name_ip}_plus_pileup.bdg > ${name_ip}_plus_pileup.bed                                                 

awk '{print $0 "\t" "*" "\t" "-" }'  ${name_ip}_minus_pileup.bdg > ${name_ip}_minus_pileup.bed                                               

cat ${name_ip}_plus_pileup.bed  ${name_ip}_minus_pileup.bed > ${name_ip}_pileup_normed.bed                                                   

rm ${name_ip}_*_pileup.bed ${name_ip}_*_pileup.bdg                                                                                           
                                                                                                                                             
awk '{print $0 "\t" "*" "\t" "+" }'  ${name_input}_plus_pileup.bdg > ${name_input}_plus_pileup.bed 

awk '{print $0 "\t" "*" "\t" "-" }'  ${name_input}_minus_pileup.bdg > ${name_input}_minus_pileup.bed                                 

cat ${name_input}_plus_pileup.bed ${name_input}_minus_pileup.bed > ${name_input}_pileup_normed.bed       

rm ${name_input}_*_pileup.bed ${name_input}_*_pileup.bdg                                                                                     


## for selected peak regions 
                                                                                                                                      
intersectBed -a ${peak_bed}  -b ${name_ip}_pileup_normed.bed -s  -wao > ${name_ip}_peak_intersect.cov                  

intersectBed -a ${peak_bed}  -b ${name_input}_pileup_normed.bed -s  -wao > ${name_input}_peak_intersect.cov            

rm ${name_ip}_pileup_normed.bed  ${name_input}_pileup_normed.bed



############################################
# for peak m6A measured by max, mean and sum 
###########################################

Rscript --vanilla /cluster/projects/hansengroup/Yong/MeRIP-SAP/4_SLURM/R_scripts/all_peak_m6a_level.R  ${peak_bed}   ${name_input}_peak_intersect.cov ${name_ip}_peak_intersect.cov all_peak_m6a_level_in_${name_ip}.txt

rm ${name_ip}_peak_intersect.cov
rm ${name_input}_peak_intersect.cov


