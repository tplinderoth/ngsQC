readCleaner
===========

readCleaner is for quality control of next-generation sequencing reads in fastq format.

## Installation

readCleaner requires the Boost C++ libraries to compile. To install Boost on Ubuntu you can simply use apt-get:

	sudo apt-get install libboost-all-dev

Then once you have Boost installed, compile readCleaner:

	cd ./ngsQC/readCleaner
	make

## External program dependencies

readCleaner relies on some external program to be installed and in the user's PATH. The following programs as well as the versions that readCleaner has been tested with follow. The name enclosed in [ ] is what the name of the binary that readCleaner looks for.

* super_deduper v 2.0 [super_deduper]
* cutadapt v 1.18 [cutadapt]
* PEAR v 0.9.10 [pear]
* bowtie 2 v 2.3.4.1 [bowtie2]
* trimmomatic v 0.36 [full path to trimmomatic jar file, e.g. /usr/local/src/trimmomatic-0.36/dist/jar/trimmomatic-0.36.jar]
* FastQC v 0.11.4 [fastqc] (this is optional, as it required only for post-filtering quality evaluation)

## Usage

Run readCleaner without arguments to get usage information.

	./readCleaner

	readCleaner version 0.0.2

	Usage: readCleaner [arguments]

	EXECUTABLE DEPENDENCIES:
	super_deduper, cutadapt, pear, bowtie2, fastqc (optional)
	Note: Program executable names should be exactly as listed above

	JAR DEPENDENCIES:
	trimmomatic

	INPUT:
	-fastqdir     <string>      Directory containing fastq files []
	-outdir       <string>      Directory to place cleaned fastq files in []
	-contam       <string>      Fasta format file of potential contaminant sequences []
	-adapter      <string>      Adapter type: truseq, nextera, smallrna [truseq]
	-trimjar      <string>      Full path to Trimmomatic jar file []
	-minqual      <int>         Minimum base quality for trimming ends of reads using BWA algorithm [20]
	-minlength    <int>         Minimum read length [36]
	-dust         <double>      Maximum read DUST score for removing low complexity reads [4]
	-missing      <double>      Maximum proportion of 'N' bases in reads [0.5]
	-quality_base <int>         Phred quality score offset (int = 33 or 64) [33]
	-nthreads     <int>         Number of threads for multithreaded programs to use [1]
	-pe           <int>         Data is paired end (int = 1) or unpaired (int = 0) [1]
	-eval         <int>         Perform fastQC evaluation of cleaned reads (int = 1), or do not (int = 0) [1]

The input are fastq format files. Use `--fastqdir` to direct readCleaner to the directory where all fastq files that you want filtered exist.

The output will be quality-controlled fastq files in a directory called 'reads' named by their original prefix with '_u_final.txt' appended to unpaired read files, '_1_final.txt' appended to forward read files, and '_2_final.txt' appended to reverse reads files. Separate forward and reverse reads files will only be generated with `pe 1`. A directory labeled 'evaluation' will contain the FastQC quality evaluation files for the filtered fastq files. Lastly, a directory called 'filtered' will contain fastq files that failed quality control to due being identified as sourced from contaminants, having low-complexity, or representing duplicates. 
