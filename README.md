# gatk_basic
A very simple and basic shell script to run GATK best practice pipeline

This is a basic variant calling and annotation pipeline for beginners using haplotype caller and snpeff. 
It is written as a single shell script. 

There are two scripts included in this repository. The scripts can be found in the gatk_basic/scripts folder.

GATK.sh script runs the pipeline for variant calling and variant annotation step by step.

geneReport.sh gets the missense variants in a particular gene. 

Tools Required:
  1. BWA 
  2. Samtools  
  3. Picard
  4. Tabix
  5. GATK
  6. SnpEff

To Run:
  1. Install tools to be required
  2. Create a directory and place the shell script(GATK.sh) inside the directory
  3. Place samples inside the directory or you can specify the path of the file.
  4. Provide samples name manually inside the script,
      if single end -> sample="read.fq"
      if paired end -> sample="read_1.fq read_2.fq"
  5. Download and Provide the corresponding reference genome for your data coordinates
  6. Specify the path of jar file for Picard and GATK tools.
  7. Select corresponding SnpEff database for your data (eg: hg38, hg19)

  command to run the script:
  -------------------------

  chmod +x GATK.sh

  ./GATK.sh

  Gene Report:
  -----------
  --> For searching your gene of interest in annotated vcf file and showing results.
  --> while running this script, you have to give vcf file path and gene of interest (one or many as comma separated value).
  
  Input --> Annotated vcf file, Gene of Interest

  Usage:
  -----
  chmod +x geneReport.sh

  ./geneReport.sh
