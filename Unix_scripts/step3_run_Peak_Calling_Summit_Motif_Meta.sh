#!/bin/bash
#SBATCH -p himem
#SABTCH -t 3-0:0:0
#SBATCH -e ./log/%x-%j.errout
#SBATCH -o ./log/%x-%j.sdout
#SBATCH --mem=61440M

#########################
cd ${out_dir}

#################
## reference data
#################
script=/cluster/projects/hansengroup/Yong/MeRIP-SAP/4_SLURM/R_scripts                # for R scripts
db=/cluster/projects/hansengroup/Yong/DB/hg38_ek12                      # for the DB
gtf=${db}/hg38_v25_ek12.gtf
CDS_12=${db}/Gencode_v24_CDS_12.bed 						# for CDS region only 
fa=${db}/hg38_ek12.fa
protein_list=${db}/protein_coding_gene_list.txt
lincRNA_list=${db}/lincRNA_gene_list.txt
ref_txdb=${db}/Guitar_TxDb_hg38.Rdata   # for guitar plot


#########################################
## initialization
## call peak with or without duplications
#########################################

dup_set=yes            #  remove duplication by default   

if [ ${rmdup} = ${dup_set} ]              ## white space is mandatory
then
	ip=${bam_dir}/${ip_name}_rmdup.bam
	input=${bam_dir}/${input_name}_rmdup.bam;
	ip_bw=${bam_dir}/${ip_name}_rmdup_pileup_normed.bw
    input_bw=${bam_dir}/${input_name}_rmdup_pileup_normed.bw
else 
	ip=${bam_dir}/${ip_name}_sorted_transcriptome.bam;
 	input=${bam_dir}/${input_name}_sorted_transcriptome.bam;
	ip_bw=${bam_dir}/${ip_name}_sorted_transcriptome_pileup_normed.bw
    input_bw=${bam_dir}/${input_name}_sorted_transcriptome_pileup_normed.bw
fi


#############################################
##  metagene plot for IP and input comparison
##  pereach sample
#############################################

module load deeptools/3.2.1

## 5UTR: mean/3rd Qu : 200nt      3UTR: mean/3rd Qu : ~750
computeMatrix scale-regions  --metagene  --skipZeros -b 1000 -m 4000  -a 2000  -R ${CDS_12}  -S ${ip_bw} ${input_bw}  -o ${name}_CM.gz
plotProfile -m ${name}_CM.gz  --perGroup  --startLabel "Start codon"  --endLabel  "Stop codon"  -out  ${name}_signal_metagene.pdf  

module unload deeptools/2.3.4      ## Inconsistent Phtyon version might cause mess


############################# 
## Metpeak: withdup or nodup
#############################

module load R/3.3.0

## call peak
## install MeTPeak under private R library __>>Following  MeTPeak GitHub instructions (R/3.3.0) 
Rscript --vanilla ${script}/metpeak.R $gtf $ip  $input  $name

## sort peaks
cp ./${name}/peak.bed  ./${name}_all_peaks.bed                
cp ./${name}/peak.xls  ./${name}_all_peaks.xls                
cp ./${name}/metpeak.Rdata  ./${name}_all_peaks.Rdata                
rm -r ${name}

## sort and peaks in human chromosome                                                                                                   
cp ${name}_all_peaks.xls ${name}_peaks_c.xls                                         
sed '/chr/!d' ${name}_peaks_c.xls  | sed '/chrM/d' | sed '/chrEK12/d' > ${name}_peaks_f.xls                            
sort -k 5 -g -r ${name}_peaks_f.xls  >  ${name}_peaks_sorted.xls                      
rm ${name}_peaks_c.xls  ${name}_peaks_f.xls                    
cp ${name}_all_peaks.bed ${name}_peaks_c.bed                                         
sed '/chr/!d' ${name}_peaks_c.bed  | sed '/chrM/d' | sed '/chrEK12/d' > ${name}_peaks_f.bed 
sort -k 5 -g  ${name}_peaks_f.bed  >  ${name}_peaks_sorted.bed    ## Ascending score,  as the score are positively correlated with p-value               
rm ${name}_peaks_c.bed   ${name}_peaks_f.bed 

awk -v OFS="\t" ' {i=i+1; $9=i; print $0}' ${name}_peaks_sorted.bed  > ${name}_peaks_sorted_1.bed

# transfer to gene_name
gene_list=${db}/all_gene_list.txt
Rscript --vanilla  ${script}/gene_id2name.R   ${name}_peaks_sorted_1.bed $gene_list  ${name}_peaks_sorted_named.bed  ${name}
rm ${name}_peaks_sorted_1.bed


############################
#i#  peaks in protein/lncRNA
############################i

# for protein coding gene
Rscript --vanilla ${script}/peak_summary.R $protein_list ${name}_peaks_sorted.bed ${name}_peaks_pc_gene_summary.txt

# for lincRNA
Rscript --vanilla ${script}/peak_summary.R $lincRNA_list ${name}_peaks_sorted.bed ${name}_peaks_lincRNA_gene_summary.txt


##################################
##  peaksmmit based ip/input ratio  
##################################

module load ucsctools/378
module load bedtools/2.27.1

ip_p=${bam_dir}/${ip_name}_plus_pileup_normed.bw
ip_m=${bam_dir}/${ip_name}_minus_pileup_normed.bw
input_p=${bam_dir}/${input_name}_plus_pileup_normed.bw
input_m=${bam_dir}/${input_name}_minus_pileup_normed.bw

## bw to bdg
bigWigToBedGraph $ip_p ${ip_name}_plus_pileup.bdg
bigWigToBedGraph $ip_m ${ip_name}_minus_pileup.bdg

bigWigToBedGraph $input_p ${input_name}_${ip_name}_plus_pileup.bdg        ## Avoiding mess in cases that sharing the same input
bigWigToBedGraph $input_m ${input_name}_${ip_name}_minus_pileup.bdg

## add strand information: for ++ --
## for IP
awk '{print $0 "\t" "*" "\t" "+" }'  ${ip_name}_plus_pileup.bdg > ${ip_name}_plus_pileup.bed
awk '{print $0 "\t" "*" "\t" "-" }'  ${ip_name}_minus_pileup.bdg > ${ip_name}_minus_pileup.bed
cat ${ip_name}_plus_pileup.bed  ${ip_name}_minus_pileup.bed > ${ip_name}_pileup_normed.bed
rm ${ip_name}_*_pileup.bed ${ip_name}_*_pileup.bdg

## for Input 
awk '{print $0 "\t" "*" "\t" "+" }'  ${input_name}_${ip_name}_plus_pileup.bdg > ${input_name}_${ip_name}_plus_pileup.bed
awk '{print $0 "\t" "*" "\t" "-" }'  ${input_name}_${ip_name}_minus_pileup.bdg > ${input_name}_${ip_name}_minus_pileup.bed

cat ${input_name}_${ip_name}_plus_pileup.bed  ${input_name}_${ip_name}_minus_pileup.bed > ${input_name}_${ip_name}_pileup_normed.bed

rm ${input_name}_${ip_name}*_pileup.bed ${input_name}_${ip_name}_*_pileup.bdg
bedGraphToBigWig  ${name}_minus_pileup_normed.bdg  $chrlen  ${name}_minus_pileup_normed.bw


## coverages 
intersectBed -a ${name}_peaks_sorted_named.bed  -b ${ip_name}_pileup_normed.bed -s  -wao > ${ip_name}_peak_intersect.cov

intersectBed -a ${name}_peaks_sorted_named.bed  -b ${input_name}_${ip_name}_pileup_normed.bed -s  -wao > ${input_name}_${ip_name}_peak_intersect.cov

rm ${ip_name}_pileup_normed.bed  ${input_name}_${ip_name}_pileup_normed.bed


#####################################
##  for ratio summit and distribution 
#####################################

Rscript --vanilla ${script}/ip_summit_ratio.R ${input_name}_${ip_name}_peak_intersect.cov ${ip_name}_peak_intersect.cov ${name}_peaks_summit.bed

rm ${input_name}_${ip_name}_peak_intersect.cov
rm ${ip_name}_peak_intersect.cov

## summit to genomefeatures
summit=${name}_peaks_summit.bed


############  distribution in 5 non-overlapping genomic features ############
intersectBed -wa -wb  -a  $summit  -b  ${db}/TSS_200.bed  ${db}/5UTR_trimmed.bed ${db}/CDS.bed  ${db}/stop_codon_up_200.bed ${db}/stop_codon.bed ${db}/stop_codon_down_200.bed  ${db}/3UTR.bed -names TSS_200 5UTR CDS stop_codon_up_200 stop_codon stop_codon_down_200 3UTR > ${name}_GenomicFeatures.int

echo Genomic_Feature    Peaks_cnt    Region_len    Density      ${name} > ${name}_GenomicFeatures.stat

awk '{if($7=="TSS_200"){print $9 "\t" $10 "\t" $5 }}' ${name}_GenomicFeatures.int |  sort -k3 -u | awk '{len=$2-$1; sum+=len} END {print "TSS_200" "\t" NR "\t" sum "\t " NR/sum}' >> ${name}_GenomicFeatures.stat

awk '{if($7=="5UTR"){print $9 "\t" $10 "\t" $5 }}' ${name}_GenomicFeatures.int | sort -k3 -u | awk '{len=$2-$1; sum+=len} END {print "5UTR" "\t" NR "\t" sum "\t" NR/sum}' >> ${name}_GenomicFeatures.stat

awk '{if($7=="CDS"){print $9 "\t" $10 "\t" $5}}' ${name}_GenomicFeatures.int | sort  -k3 -u| awk '{len=$2-$1; sum+=len} END {print "CDS" "\t"  NR "\t" sum "\t" NR/sum} ' >> ${name}_GenomicFeatures.stat

awk '{if($7=="stop_codon_up_200"){print $9 "\t" $10 "\t" $5 }}' ${name}_GenomicFeatures.int | sort -k3 -u | awk '{len=$2-$1; sum+=len} END {print "stop_codon_up_200" "\t"  NR "\t" sum "\t" NR/sum}' >> ${name}_GenomicFeatures.stat

awk '{if($7=="stop_codon"){print $9 "\t" $10 "\t" $5 }}' ${name}_GenomicFeatures.int |sort -k3 -u| awk '{len=$2-$1; sum+=len} END {print "stop_codon" "\t"  NR "\t" sum "\t" NR/sum}' >> ${name}_GenomicFeatures.stat

awk '{if($7=="stop_codon_down_200"){print $9 "\t" $10 "\t" $5 }}' ${name}_GenomicFeatures.int | sort -k3 -u | awk '{len=$2-$1; sum+=len} END {print "stop_codon_down_200" "\t"  NR "\t" sum "\t" NR/sum}' >> ${name}_GenomicFeatures.stat

awk '{if($7=="3UTR"){print $9 "\t" $10 "\t" $5 }}' ${name}_GenomicFeatures.int | sort -k3 -u | awk '{len=$2-$1; sum+=len} END {print "3UTR" "\t"  NR "\t" sum "\t" NR/sum}' >> ${name}_GenomicFeatures.stat


##################
## DREME for motif 
##  100 nt
##################

module load meme/4.10.2
module load ucsctools/378
module load bedtools/2.27.1

for peak_reg in 100 150 200      ## flanking region len
do

for peak_cnt in 1000 2000 5000
do

head -n ${peak_cnt} ${name}_peaks_summit.bed | awk -v len="${peak_reg}"  -v OFS="\t" '{$2=$2+1-len/2; $3=$3+len/2; print $0 }'  > ${name}_${peak_cnt}_summit_${peak_reg}.bed
fastaFromBed -s -fi $fa  -bed ${name}_${peak_cnt}_summit_${peak_reg}.bed  -fo ${name}_${peak_cnt}_summit_${peak_reg}.fa
cp ${name}_${peak_cnt}_summit_${peak_reg}.fa  ${name}_${peak_cnt}_summit_${peak_reg}_rna.fa 
sed -ie 's/T\|t/U/g' ${name}_${peak_cnt}_summit_${peak_reg}_rna.fa

## no argument of -rna for MEME 4.10.2 ..
dreme -p ${name}_${peak_cnt}_summit_${peak_reg}_rna.fa  -eps   -m 10 -norc -oc  ${name}_${peak_cnt}_dreme_${peak_reg}

rm ${name}_${peak_cnt}_summit_${peak_reg}.fa

centrimo --norc -oc ${name}_${peak_cnt}_centrimo_${peak_reg} ${name}_${peak_cnt}_summit_${peak_reg}_rna.fa  ./${name}_${peak_cnt}_dreme_${peak_reg}/dreme.txt

grep BEST ./${name}_${peak_cnt}_dreme_${peak_reg}/dreme.txt | awk -v cnt="$peak_cnt" '{print $0 "\t" $4/cnt}' > ${name}_${peak_cnt}_dreme_${peak_reg}_summary.txt

rm ${name}_${peak_cnt}_summit_${peak_reg}*

done

done


#############################################################################
############ m6A peak summits distribution : using GUITAR Plot ##############

module load R/3.3.0
summit=${name}_peaks_summit.bed
Rscript --vanilla ${script}/guitar_plot.R $summit  ${name}  ${ref_txdb}

!
