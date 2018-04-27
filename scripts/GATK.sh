#! /bin/bash


#########################INPUT FILES#########################################
#
# sample : single or paired end sample(.fq) files
# ref : Reference genome (.fa)
# Knownsites : dbsnp version - hg38
# snpeff_db : hg38
#############################################################################

sample="liver.fastq"
ref="hg38.fa"
picard="java -jar /path/to/picard.jar"
GATK="java -jar /path/to/GenomeAnalysisTK.jar"
knownsites="/path/to/dbsnp_138.hg38.vcf.gz"
snpeff_db="hg38"
gene=""

#############################################################################

## STEP-1: Alignment -- BWA
## -M : Need to provide the -M flag to BWA, this tells it to consider split reads as secondary, need this for GATK variant calling/Picard support.
## -R : Readgroup info is provided with the -R flag. This information is key for downstream GATK functionality. GATK will not work without a read group tag.
echo "BWA Indexing ..."

bwa index $ref

bwa mem -M -R '@RG\tID:sample_1\tLB:sample_1\tPL:ILLUMINA\tPM:HISEQ\tSM:sample_1' $ref $sample > aligned_reads.sam

## STEP-2: Picard Tools
## Sort SAM file by coordinate, convert to BAM

$picard SortSam INPUT=aligned_reads.sam OUTPUT=sorted_reads.bam SORT_ORDER=coordinate

## STEP-3: Collect Alignment & Insert Size Metrics (optional)
## output files --> alignment_metrics.txt, insert_metrics.txt, insert_size_histogram.pdf, depth_out.txt

$picard CollectAlignmentSummaryMetrics R=$ref I=sorted_reads.bam O=alignment_metrics.txt

$picard CollectInsertSizeMetrics INPUT=sorted_reads.bam OUTPUT=insert_metrics.txt HISTOGRAM_FILE=insert_size_histogram.pdf


## STEP-4: Mark Duplicates

$picard MarkDuplicates INPUT=sorted_reads.bam OUTPUT=dedup_reads.bam METRICS_FILE=metrics.txt

## STEP-5: Build BAM Index

$picard BuildBamIndex INPUT=dedup_reads.bam

$picard CreateSequenceDictionary \
        R=$ref \
	O=$ref.dict

samtools faidx $ref

tabix -f -p vcf $knownsites

$gatk -T BaseRecalibrator \
        -I dedup_reads.bam \
        -R $ref \
        -knownSites $knownsites \
	-o recal_data.table

$gatk -T PrintReads \
        -R $ref \
        -I dedup_reads.bam \
        -BQSR recal_data.table \
	-o recal_reads.bam


gatk -T HaplotypeCaller \
        -R $ref \
        -I recal_reads.bam \
	-o raw_variants.vcf

gatk -T SelectVariants \
        -R $ref \
        -V raw_variants.vcf \
        -selectType SNP \
        -o raw_snps.vcf

gatk -T SelectVariants \
        -R $ref \
        -V raw_variants.vcf \
        -selectType INDEL \
	-o raw_indels.vcf

gatk -T VariantFiltration \
        -R $ref \
        -V raw_snps.vcf \
        -filterName "QD_filter" \
        -filter "QD'<'2.0" \
        -filterName "FS_filter" \
        -filter "FS'>'60.0" \
        -filterName "MQ_filter" \
        -filter "MQ'<'40.0" \
        -filterName "SOR_filter" \
        -filter "SOR'>'4.0" \
	-o filtered_snps.vcf

gatk -T VariantFiltration \
        -R $ref \
        -V raw_indels.vcf \
        -filterName "QD_filter" \
        -filter "QD'<'2.0" \
        -filterName "FS_filter" \
        -filter "FS'>'200.0" \
        -filterName "SOR_filter" \
        -filter "SOR'>'10.0" \
	-o filtered_indels.vcf


java -jar snpEff.jar -v $snpeff_db filtered_snps.vcf > filtered_snps_ann.vcf

java -jar snpEff.jar -v $snpeff_db filtered_snps.vcf > filtered_snps_ann.vcf
