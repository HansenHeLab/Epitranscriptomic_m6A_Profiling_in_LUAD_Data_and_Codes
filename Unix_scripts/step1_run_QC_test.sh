#!/bin/bash

################################################
## QC test by subsmapling 1M reads from fastq.gz 
## exreact strand-specific information


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
module load cutadapt/2.5 
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
genecode_bed=${db}/hg38_gencodev24.bed
index=${db}/star_index

## CPC MeRip-seq: TruSeq CD Indexes:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[i7]ATCTCGTATGCCGTCTTCTGCTTG, length = 65

adapter=AGATCGGAAGAGC         # for common sequence for mate1 and mate2
min_len=50	              # minimum length for the left reads


###########################
# remove adapter and FASTQC
###########################

echo Running cutadpat ...

## 1M reads for testing 
head -n 1000000 ${mate1}  > ${name}_R1.fastq.gz
head -n 1000000 ${mate2}  > ${name}_R2.fastq.gz

cutadapt -a $adapter -m ${min_len} -o  ${name}_R1_rmadp.fastq.gz  ${name}_R1.fastq.gz
cutadapt -a $adapter -m ${min_len} -o  ${name}_R2_rmadp.fastq.gz  ${name}_R2.fastq.gz

rm ${name}_R1.fastq.gz ${name}_R2.fastq.gz

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

rm ${name}*_Log.final.out


#####################
# RSeQC after mapping
#####################
module load RSeQC/3.0.1

samtools index ${name}_R1_Aligned.sortedByCoord.out.bam
samtools index ${name}_R2_Aligned.sortedByCoord.out.bam

## whole genome bam to transcriptome bam files
module load bedtools 

intersectBed -abam ${name}_R1_Aligned.sortedByCoord.out.bam -b ${gtf_exon} > ${name}_R1_Aligned.sortedByCoord.out.transcriptome.bam
intersectBed -abam ${name}_R2_Aligned.sortedByCoord.out.bam -b ${gtf_exon} > ${name}_R2_Aligned.sortedByCoord.out.transcriptome.bam

samtools index ${name}_R1_Aligned.sortedByCoord.out.transcriptome.bam
samtools index ${name}_R2_Aligned.sortedByCoord.out.transcriptome.bam

## strand-specific?
echo "" >> ${name}_mapping_summary.txt
echo "" >> ${name}_mapping_summary.txt

## strand-specific?
echo "############ strand-test for R1 reads #############"    >> ${name}_mapping_summary.txt
infer_experiment.py -r ${gtf_bed} -i ${name}_R1_Aligned.sortedByCoord.out.bam >> ${name}_mapping_summary.txt

echo "" >> ${name}_mapping_summary.txt
echo "############ strand-test for R2 reads #############"    >> ${name}_mapping_summary.txt
infer_experiment.py -r ${gtf_bed} -i ${name}_R2_Aligned.sortedByCoord.out.bam >> ${name}_mapping_summary.txt


