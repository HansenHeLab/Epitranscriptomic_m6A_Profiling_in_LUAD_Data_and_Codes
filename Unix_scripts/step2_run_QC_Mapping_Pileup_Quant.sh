#!/bin/bash
#SBATCH -p himem
#SABTCH -t 3-0:0:0
#SBATCH -e ./log/%x-%j.errout
#SBATCH -o ./log/%x-%j.sdout
#SBATCH --mem=61440M

################# 
cd ${out_dir}

#################
## initial set up   
#################

## Load tools ##################################
module load cutadapt/2.5	# for h4h
module load fastqc/0.11.5
module load STAR/2.4.2a
module load R/3.5.0
module load samtools/1.3.1


## Assign Reference and INDEX ##############################################
db=/cluster/projects/hansengroup/Yong/DB/hg38_ek12	# DB on h4h
fa=${db}/hg38_ek12.fa
gtf=${db}/hg38_v25_ek12.gtf
gtf_bed=${db}/hg38_v25_ek12_sorted.bed
gtf_exon=${db}/hg38_v25_ek12_sorted_exons.bed
gtf_p=${db}/protein_coding_sorted.bed      ## for protein coding gene only
genecode_bed=${db}/hg38_gencodev24.bed
index=${db}/star_index
chrlen=${db}/star_index/chrNameLength.txt

## CPC MeRip-seq: TruSeq CD Indexes:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[i7]ATCTCGTATGCCGTCTTCTGCTTG, length = 65
adapter=AGATCGGAAGAGC         # for common sequence for mate1 and mate2
min_len=50	              # minimum length for the left reads
ext_size=200		      ## pileup extend size


###########################
# remove adapter and FASTQC
###########################

cutadapt -a $adapter -m ${min_len} -o  ${name}_R1_rmadp.fastq.gz  ${mate1}
cutadapt -a $adapter -m ${min_len} -o  ${name}_R2_rmadp.fastq.gz  ${mate2}

mate1_f=${name}_R1_rmadp.fastq.gz
mate2_f=${name}_R2_rmadp.fastq.gz

fastqc  ${mate1_f}
fastqc  ${mate2_f}

#########################################################
# run star2 for alignment: uniquely mapped reads outputed
#########################################################

## for Mate 1
STAR --runMode alignReads --runThreadN 10  \
     --genomeDir $index  --sjdbGTFfile $gtf  --readFilesCommand gunzip -c  \
     --readFilesIn ${mate1_f}      --clip3pNbases 0       --clip5pNbases 0 \
     --outFilterMultimapNmax 1    --outFilterMismatchNoverReadLmax 0.04  --outSAMmapqUnique 255 \
     --outFileNamePrefix ${name}_R1_   \
     --outSAMtype BAM SortedByCoordinate  --outSAMattributes All 

STAR --runMode alignReads --runThreadN 10  \
     --genomeDir $index  --sjdbGTFfile $gtf   --readFilesCommand gunzip -c\
     --readFilesIn ${mate2_f}      --clip3pNbases 0       --clip5pNbases 0 \
     --outFilterMultimapNmax 1  --outFilterMismatchNoverReadLmax 0.04    --outSAMmapqUnique 255 \
     --outFileNamePrefix ${name}_R2_    \
     --outSAMtype BAM SortedByCoordinate  --outSAMattributes All 

rm ${mate1_f}  ${mate2_f}			## rm reads after removing adapter
rm -r ${name}_*_STARgenome ${name}_*_STARtmp;	
rm  ${name}_*Log.out ${name}_*.progress.out  ${name}_*SJ.out.tab

##################
## mapping summary
echo "############  R1 mapping summary  ############# " > ${name}_mapping_summary.txt
cat  ${name}_R1_Log.final.out  >>  ${name}_mapping_summary.txt

echo "" >> ${name}_mapping_summary.txt
echo "" >> ${name}_mapping_summary.txt
echo "############  R2 mapping summary  ############# " >>  ${name}_mapping_summary.txt
cat  ${name}_R2_Log.final.out  >>  ${name}_mapping_summary.txt


#####################
# RSeQC after mapping
#####################

module load RSeQC/3.0.1

samtools index ${name}_R1_Aligned.sortedByCoord.out.bam
samtools index ${name}_R2_Aligned.sortedByCoord.out.bam

## strand-specific?
echo "" >> ${name}_mapping_summary.txt
echo "" >> ${name}_mapping_summary.txt

## strand-specific?
echo "############ strand-test for R1 reads #############"    >> ${name}_mapping_summary.txt
infer_experiment.py -r ${gtf_bed} -i ${name}_R1_Aligned.sortedByCoord.out.bam >> ${name}_mapping_summary.txt

echo "" >> ${name}_mapping_summary.txt
echo "############ strand-test for R2 reads #############"    >> ${name}_mapping_summary.txt
infer_experiment.py -r ${gtf_bed} -i ${name}_R2_Aligned.sortedByCoord.out.bam >> ${name}_mapping_summary.txt


##########################################
## pool pair-end read: for strand-specific 
##########################################

################################
## stand specific infomation ###
strand_dir=FR                  # R1R2: FR 

if [ ${strand} == ${strand_dir} ]            ## white space is mandatory
then 
	mate_yes=${name}_R1_Aligned.sortedByCoord.out.bam
	mate_rev=${name}_R2_Aligned.sortedByCoord.out.bam
else
	mate_yes=${name}_R2_Aligned.sortedByCoord.out.bam
	mate_rev=${name}_R1_Aligned.sortedByCoord.out.bam
fi

##################################################################
## switch the flag for the revers mapped mate of the pair-end reads
## switching flag for reverse mate
samtools view -h ${mate_rev} > ${name}_rev.sam
sed 's/0\tchr/1\tchr/' ${name}_rev.sam > ${name}_rev_tmp1.sam
sed 's/16\tchr/0\tchr/' ${name}_rev_tmp1.sam > ${name}_rev_tmp2.sam
rm ${name}_rev_tmp1.sam  ${name}_rev.sam
sed 's/1\tchr/16\tchr/' ${name}_rev_tmp2.sam > ${name}_rev_switched.sam
rm ${name}_rev_tmp2.sam
samtools view -bS ${name}_rev_switched.sam > ${name}_rev_switched.bam
rm ${name}_rev_switched.sam

##  merging bam
samtools merge ${name}_merged.bam  ${mate_yes}  ${name}_rev_switched.bam
rm ${mate_yes}  ${mate_rev}  ${name}_rev_switched.bam  ${mate_yes}.bai  ${mate_rev}.bai
samtools sort ${name}_merged.bam -o ${name}_sorted.bam
samtools index ${name}_sorted.bam
rm ${name}_merged.bam


################
## extrat EK12
################

samtools view -H ${name}_sorted.bam > ${name}_ek12.bam
samtools view  ${name}_sorted.bam chrEK12 >> ${name}_ek12.bam


############################################
## retrive transcriptome reads only 
## remove duplications and split strand bams
############################################

module load bedtools/2.27.1

intersectBed -abam ${name}_sorted.bam -b ${gtf_exon} > ${name}_sorted_transcriptome.bam
samtools index ${name}_sorted_transcriptome.bam

## remove duplication
samtools rmdup -s ${name}_sorted_transcriptome.bam  ${name}_rmdup.bam
samtools index ${name}_rmdup.bam

## split into corresponding strand
## minus strand
samtools view -H ${name}_rmdup.bam > ${name}_rmdup_minus.bam
samtools view  -f 16  ${name}_rmdup.bam >> ${name}_rmdup_minus.bam
samtools sort ${name}_rmdup_minus.bam -o ${name}_rmdup_minus_sorted.bam
samtools index ${name}_rmdup_minus_sorted.bam

## plus strand
samtools view -H ${name}_rmdup.bam > ${name}_rmdup_plus.bam
samtools view -F 16   ${name}_rmdup.bam >> ${name}_rmdup_plus.bam
samtools sort ${name}_rmdup_plus.bam -o ${name}_rmdup_plus_sorted.bam
samtools index ${name}_rmdup_plus_sorted.bam

rm ${name}_rmdup_plus.bam ${name}_rmdup_minus.bam


##################
## mapping summary
##################

echo "" >> ${name}_mapping_summary.txt
echo "" >> ${name}_mapping_summary.txt

# echo "###### all_input_reads_rmadp\all_uniqe_mapped_reads\rmdup_reads\rmdup_plus_reads\rmdup_minus_reads\ek12_reads#########"   >> ${name}_mapping_summary.txt

echo all_input_reads_rmadp: `grep "Number of input reads"  ${name}_R*_Log.final.out | awk '{sum += $7} END {print sum}'`   >>  ${name}_mapping_summary.txt

echo all_uniquely_mapped_read: `samtools view -c ${name}_sorted.bam` >> ${name}_mapping_summary.txt

echo rmdup_reads: `samtools view -c ${name}_rmdup.bam`  >> ${name}_mapping_summary.txt

echo rmdup_plus_strand_reads: `samtools view -c ${name}_rmdup_plus_sorted.bam` >> ${name}_mapping_summary.txt

echo rmdup_minus_strand_reads: `samtools view -c ${name}_rmdup_minus_sorted.bam` >> ${name}_mapping_summary.txt

echo ek12_reads: `samtools view  -c ${name}_ek12.bam` >> ${name}_mapping_summary.txt


################################
## pileup in normalized bw files
################################

module load MACS/2.2.5
module load ucsctools/378
module load RSeQC/3.0.1

## bdg pileups 
macs2 pileup -i ${name}_sorted_transcriptome.bam -o ${name}_sorted_pileup1.bdg --outdir ./ -f BAM --extsize ${ext_size}
macs2 pileup -i ${name}_rmdup.bam -o ${name}_rmdup_pileup1.bdg --outdir ./ -f BAM --extsize ${ext_size}
macs2 pileup -i ${name}_rmdup_minus_sorted.bam -o ${name}_minus_pileup1.bdg --outdir ./ -f BAM --extsize ${ext_size}
macs2 pileup -i ${name}_rmdup_plus_sorted.bam -o ${name}_plus_pileup1.bdg --outdir ./ -f BAM --extsize ${ext_size}


## fillting ambiguous reads 
grep  "chr" ${name}_sorted_pileup1.bdg > ${name}_sorted_pileup2.bdg
grep  "chr" ${name}_rmdup_pileup1.bdg > ${name}_rmdup_pileup2.bdg
grep  "chr" ${name}_minus_pileup1.bdg > ${name}_minus_pileup2.bdg
grep  "chr" ${name}_plus_pileup1.bdg > ${name}_plus_pileup2.bdg

grep -v "_\|chrEBV" ${name}_sorted_pileup2.bdg > ${name}_sorted_pileup3.bdg
grep -v "_\|chrEBV" ${name}_rmdup_pileup2.bdg > ${name}_rmdup_pileup3.bdg
grep -v "_\|chrEBV" ${name}_minus_pileup2.bdg > ${name}_minus_pileup3.bdg
grep -v "_\|chrEBV" ${name}_plus_pileup2.bdg > ${name}_plus_pileup3.bdg

sort -k1,1  -k2n,2  ${name}_sorted_pileup3.bdg > ${name}_sorted_pileup.bdg
sort -k1,1  -k2n,2  ${name}_rmdup_pileup3.bdg > ${name}_rmdup_pileup.bdg
sort -k1,1  -k2n,2  ${name}_minus_pileup3.bdg > ${name}_minus_pileup.bdg
sort -k1,1  -k2n,2  ${name}_plus_pileup3.bdg > ${name}_plus_pileup.bdg

rm ${name}_*_pileup1.bdg ${name}_*_pileup2.bdg ${name}_*_pileup3.bdg

## bdg to bw 
bedGraphToBigWig  ${name}_sorted_pileup.bdg  $chrlen  ${name}_sorted_transcriptome_pileup.bw
bedGraphToBigWig  ${name}_rmdup_pileup.bdg  $chrlen  ${name}_rmdup_pileup.bw
bedGraphToBigWig  ${name}_minus_pileup.bdg  $chrlen  ${name}_minus_pileup.bw
bedGraphToBigWig  ${name}_plus_pileup.bdg  $chrlen  ${name}_plus_pileup.bw

rm ${name}_*_pileup.bdg

## normalization to reads/1M
normalize_bigwig.py  -i ${name}_sorted_transcriptome_pileup.bw -o ${name}_sorted_transcriptome_pileup_normed.bdg
normalize_bigwig.py  -i ${name}_rmdup_pileup.bw -o ${name}_rmdup_pileup_normed.bdg
normalize_bigwig.py  -i ${name}_plus_pileup.bw -o ${name}_plus_pileup_normed.bdg
normalize_bigwig.py  -i ${name}_minus_pileup.bw -o ${name}_minus_pileup_normed.bdg

## bdg to bw
bedGraphToBigWig  ${name}_sorted_transcriptome_pileup_normed.bdg  $chrlen  ${name}_sorted_transcriptome_pileup_normed.bw
bedGraphToBigWig  ${name}_rmdup_pileup_normed.bdg  $chrlen  ${name}_rmdup_pileup_normed.bw
bedGraphToBigWig  ${name}_plus_pileup_normed.bdg  $chrlen  ${name}_plus_pileup_normed.bw
bedGraphToBigWig  ${name}_minus_pileup_normed.bdg  $chrlen  ${name}_minus_pileup_normed.bw

rm ${name}_*_pileup_normed.bdg;
rm ${name}_*_pileup.bw


#############################
## bw based metagene analysis
#############################

# module load deeptools/3.2.1
computeMatrix scale-regions  -R ${gtf_p} -S ${name}_sorted_pileup_normed.bw ${name}_rmdup_pileup_normed.bw -o ${name}


#############################
## HTseq for gene reads count
############################
module purge                  ## unload previously loaded python  
module load HTSeq/0.11.0

## based on merged alignment bam files
## have alreadly switched for reverse mate,  -s yes

htseq-count -t gene -i gene_name -f bam -r pos -s yes  -q ${name}_sorted_transcriptome.bam  $gtf  > ${name}_htseq_count


#############################################
## bam2cram: saving nearly half storage space
############################################
samtools view -T $fa -C -o ${name}_sorted.cram  ${name}_sorted.bam
rm ${name}_sorted.bam
