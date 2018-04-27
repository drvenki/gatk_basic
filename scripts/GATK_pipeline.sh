#! /bin/bash

#########################INPUT FILES#########################################
# sample : single or paired end sample(.fq) files
# ref : Reference genome (.fa)
# Knownsites : dbsnp version - hg38
# snpeff_db : hg38
#############################################################################
# input fastq files. SE: sample="read.fastq", PE: sample="read_1.fastq read_2.fastq", collapsed PE in single file: sample="reads.fastq"
sample="read_1.fq read_2.fq"
ref="hg38.fa"
picard="java -jar /path/to/picard.jar"
gatk="java -jar /path/to/GenomeAnalysisTK.jar"
knownsites="/path/to/dbsnp_138.hg38.vcf.gz"
snpEff="java -jar /path/to/snpEff.jar"
snpSift="java -jar /path/to/SnpSift.jar"
snpeff_db="GRCh38.86"

#############################################################################
# start of shell script

#############################################################################

## Fastq raw data Quality Check
#  FastQC reads a set of sequence files and produces from each one a quality control report consisting of a number of different modules, each one of which will help to identify a different potential type of problem in your data.
#
date
mkdir qc_report
fastqc -o qc_report $sample

## STEP-1: Alignment -- BWA
## Create algorithm indexes for reference fa file(.fai) to align reads
date
bwa index $ref

## -M : Need to provide the -M flag to BWA, this tells it to consider split reads as secondary, need this for GATK variant calling/Picard support.
## -R : Readgroup info is provided with the -R flag. This information is key for downstream GATK functionality. GATK will not work without a read group tag. For platforms other than illumina, please provide appropriate info.
## for example for pacbio \tPL:PACBIO\tPM:SMRT\ or for oxford nanopore \tPL:OXNANOPORE\tPM:GRIDION\

date
bwa mem -M -R '@RG\tID:sample_1\tLB:sample_1\tPL:ILLUMINA\tPM:HISEQ\tSM:sample_1' $ref $sample > aligned_reads.sam

## STEP-2: Picard Tools
## Sort SAM file by coordinate, convert to BAM
date
$picard SortSam INPUT=aligned_reads.sam OUTPUT=sorted_reads.bam SORT_ORDER=coordinate

## STEP-3: Collect Alignment & Insert Size Metrics (optional)
## output files --> alignment_metrics.txt, insert_metrics.txt, insert_size_histogram.pdf, depth_out.txt
date
$picard CollectAlignmentSummaryMetrics R=$ref I=sorted_reads.bam O=alignment_metrics.txt

$picard CollectInsertSizeMetrics INPUT=sorted_reads.bam OUTPUT=insert_metrics.txt HISTOGRAM_FILE=insert_size_histogram.pdf


## STEP-4: Mark Duplicates
## This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA.
## Duplicates can arise during sample preparation e.g. library construction using PCR.
date
$picard MarkDuplicates INPUT=sorted_reads.bam OUTPUT=dedup_reads.bam METRICS_FILE=metrics.txt

## STEP-5: Build BAM Index
date
$picard BuildBamIndex INPUT=dedup_reads.bam

## STEP-6: Creates a sequence dictionary for a reference sequence.
date
$picard CreateSequenceDictionary \
        R=$ref \
	O=$ref.dict

## CREATING INDEXES FOR PATTERN MATCHING TO HELP VARIANT CALLING
## STEP-7:  Create indexes for reference fa file(.fai)
date
samtools faidx $ref

## STEP-8: Create indexes for knownsites vcf file.
date
tabix -f -p vcf $knownsites

## STEP-9: Generates recalibration table for Base Quality Score Recalibration (BQSR)
date
$gatk -T BaseRecalibrator \
        -I dedup_reads.bam \
        -R $ref \
        -knownSites $knownsites \
	-o recal_data.table

## STEP-10: Print reads in the SAM/BAM/CRAM file
date
$gatk -T PrintReads \
        -R $ref \
        -I dedup_reads.bam \
        -BQSR recal_data.table \
	-o recal_reads.bam

## STEP-11: Call germline SNPs and indels via local re-assembly of haplotypes
date
$gatk -T HaplotypeCaller \
        -R $ref \
        -I recal_reads.bam \
	-o raw_variants.vcf

## STEP-12: Select SNPs from a VCF file
date
$gatk -T SelectVariants \
        -R $ref \
        -V raw_variants.vcf \
        -selectType SNP \
        -o raw_snps.vcf

## STEP-13: Select INDELs from a VCF file
date
$gatk -T SelectVariants \
        -R $ref \
        -V raw_variants.vcf \
        -selectType INDEL \
	-o raw_indels.vcf

## STEP-14: Filter variant calls(SNPs) based on INFO and/or FORMAT annotations.
## filter values for QD, FS, MQ, SOR are based on GATK best practices and could be changed to optimal values according to data
date
$gatk -T VariantFiltration \
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

## STEP-15: Filter variant calls(INDELs) based on INFO and/or FORMAT annotations
date
$gatk -T VariantFiltration \
        -R $ref \
        -V raw_indels.vcf \
        -filterName "QD_filter" \
        -filter "QD'<'2.0" \
        -filterName "FS_filter" \
        -filter "FS'>'200.0" \
        -filterName "SOR_filter" \
        -filter "SOR'>'10.0" \
	-o filtered_indels.vcf

### STEP-16: Building Database from RefSeq table from UCSC
date
$snpEff build -refSeq -v $snpeff_db

### STEP-17: variant annotation with snpEff, writing out snps and indels into separate files
date
$snpEff -v $snpeff_db filtered_snps.vcf > filtered_snps_ann.vcf

date
$snpEff -v $snpeff_db filtered_indels.vcf > filtered_indels_ann.vcf

### STEP-18: Extract missense variants
date
$snpSift filter "(ANN[*].EFFECT = 'missense_variant')" filtered_snp_refseq.vcf > missense_snp.vcf
