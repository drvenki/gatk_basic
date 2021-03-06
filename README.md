# gatk_basic
A very simple and basic shell script to run GATK best practice pipeline for Ubuntu environment.

This is a basic variant calling and annotation pipeline for beginners using the GATK haplotype caller for variant calling and snpeff for variant annotation. 
It is written as a single shell script <b>to help people new to variant calling and bioinformatics</b>. This script will get the job done for basic variant calling and annotation purpose. This pipeline is not to be used for clinical samples as this is written only for educational purposes.

There are two scripts included in this repository. The scripts can be found in the gatk_basic/scripts folder.

<b>GATK_pipeline.sh</b> script runs the pipeline for variant calling and variant annotation step by step.

<b>geneReport.sh</b> gets the missense variants in a particular gene. 

<b>Tools Required:</b> <br>
BWA, Samtools, Picard, Tabix, GATK, SnpEff

 
 
 To Run:
------
  1. Install required tools.
  2. Create a directory and place the shell scripts (GATK_pipeline.sh & geneReport.sh) and the "header.txt" file inside the same directory (> mkdir scripts)
  3. Place samples inside the directory or you can specify the path of the file (> mkdir samples)
  4. Provide samples name manually inside the script,
      if single end or a single file with collapsed paired-end reads -> sample="read.fq"
      if paired end -> sample="read_1.fq read_2.fq"
  5. Download and provide the corresponding reference genome for your data coordinates
  6. Specify the path of jar file for Picard and GATK tools.
  7. Select corresponding SnpEff database for your data (eg: hg38, hg19)
 
  
  Tools Required and Installation guide (Linux-Ubuntu):
----------------------------------------------------

Open the <b>Terminal</b> program on your Ubuntu and copy paste these commands in same sequence

  1. BWA  (Alignment tool)
    
    > sudo apt-get install bwa

  2. Samtools  (variant calling and sequence manipulation)
    
    > sudo apt-get install samtools

  3. Picard (required JAVA v1.9) (sequence data manipulation and cleaning)

    > wget https://github.com/broadinstitute/picard/releases/download/2.18.1/picard.jar
    (place this jar file in a folder where you want so that you can reuse it whenever you want to run Picard tool. PS: Picard is the tool used for getting quality metrics from the sequence data and ALSO to mark & PCR duplicates).

  4. Tabix (required JAVA v1.9) (data transformation)

    > sudo apt-get install tabix 

  5. GATK (required JAVA v1.9) (variant calling)

    Download the jar file from (https://software.broadinstitute.org/gatk/download/archive)

  6. SnpEff & SnpSift (required JAVA v1.9) (variant annotation)

    1. Download the jar file from (http://snpeff.sourceforge.net/download_donate.html)
    
    2. Building Database for reference genome from UCSC(refseq ids: NM_*, NP_*)
      * Please follow the instructions in the given link
      (http://snpeff.sourceforge.net/SnpEff_manual.html#databases) - Option 3: Building a database from RefSeq table from UCSC 

  7. FastQC (fast raw data quality control)

    > sudo apt-get install fastqc

  8. XML parser for Refseq Protein id annotation (xml parsing for fetching and extracting NP ids from NCBI)

    > sudo apt-get install libxml2-utils

Links to download annotation source files.
------------------------------

Reference genome --> http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

Knownsites for GATK-Recalibration --> ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_138.hg38.vcf.gz


  Giving permissions to the scripts to be run by users:
  -------------------------

  > chmod +x GATK_pipeline.sh
  
  > chmod +x geneReport.sh
  
  Scripts notes:
  --------------
  --> GATK_pipeline.sh is for doing QC, aligning the samples to reference, variant calling with GATK haplotype caller and annotating the variants for functional effects with snpeff.
  
  --> geneReport.sh is for searching your gene of interest in annotated vcf file which is the result of GATK_pipeline.sh script.
  
  --> you have to first run GATK_pipeline.sh script. This will produce a set of results files and out of which give annotated vcf file (gene.ann.vcf) and gene of interest (one or many as comma separated value) as Input.

  Variant Calling:
  -------------------------
  > ./GATK_pipeline.sh 
  
  OR 
  
  > bash GATK_pipeline.sh

  
  Gene wise Report:
  -----------------
  NOTE: for this step, make sure the file <b>header.txt is in the same folder as the script</b>
  > ./geneReport.sh
  
  OR
  
  > bash geneReport.sh

  The results from this script will be gene.txt, gene.report.txt, gene.missense.txt. The file gene.missense.txt will contain only missense variants in the given gene.


