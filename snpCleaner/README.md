snpCleaner
==========

snpCleaner.pl takes a VCF format file as input and filters it based on various quality metrics. Indels are ignored. It outputs a VCF and associated bed-format file of sites that pass QC as well as a VCF-like file of sites that do not pass with flags indicating which filters each site failed. For a description of the  different filters run snpCleaner.pl without arguments.

	./snpCleaner.pl
	
	#####################
	# SNPcleaner v2.4.3 #
	#####################
	This scripts works with bcftools vcf file format to filter SNPs. This script is only for 
	SNP filtering and will ignore INDELs. 
	
	#######################################
	Usage: snpCleaner.pl [options] <in.vcf>
	or
	cat <in.vcf> | snpCleaner.pl [options]
	
	Options:
	-?	 help
	-2	 keep non-biallelic sites
	-q INT   ploidy [2]
	-c INT   minimum high-quality site read depth (VCF DP4 sum) [0] 
	-d INT   minimum raw site read depth (VCF DP flag) [1]
	-D INT   maximum raw site read depth (VCF DP flag)[1000000]
	-k INT   minimum number of 'covered' individuals (requires -u) [1] 
	-u INT   minimum read depth for an individual to be considered 'covered' (requires -k) [0]
	-a INT   minimum number of high-quality alternate alleles for site [0]
	-Q INT   minimum RMS mapping quality for SNPS [10]
	-S FLOAT min p-value for strand bias [0.0001]
	-b FLOAT min p-value for base quality bias [1e-100]
	-f FLOAT min p-value for map quality bias [0]
	-e FLOAT min p-value for end distance bias [0.0001]
	-h FLOAT min p-value for exact test of HWE [0.0001]
	-F FLOAT inbreeding coefficient value [0]
	-H FLOAT min p-value for exact test of excess of heterozygous [0]
	-L FLOAT min p-value for exact test of defect of heterozygous [0]
	-A FILE  ANCESTRAL fasta file (with FAI on same folder)
	-M CHAR  mutation type(s) to remove (ex. 'CT_GA') (requires -A)
	-o FILE  where to write VCF format output of sites that pass QC, 'none' for no output [STDOUT]
	-B FILE  file of sites that passed filters
	-p FILE  name of file to dump sites that failed filters (bziped)
	-r FILE  list of contigs to exclude
	-X FILE  BED file of exonic regions (sorted from lowest to highest contig)
	-t       filter non-exonic sites (requires -X)
	-g       filter exons with SNPs out of HWE (requires -X)
	-v	 process nonvariants
	-w       print only SNP positions to VCF output
	
	Generating input VCF with BCFtools:
	-run 'bcftools mpileup' with '-a SP,DP' to provide depth and strand bias info
	-run 'bcftools call' using '-f GC -c'
	
	An example for how to generate the input VCF file:
	bcftools mpileup -f <reference fasta> -b <bam list> -a SP,DP | bcftools call -f GQ -c - > unfiltered_sites.vcf 
	
	Consider using 'bcftools call --skip-variants indels' since these sites will
	be ignored by snpCleaner anyways.
	
	Output Notes:
	-bed format file is zero-based.
	-FILTER field flags in the failed sites VCF (dumped with option -p) indicate filters that the sites
	 failed to pass and correspond to the option flags (e.g. 'S' indicates strand bias).
	
## Input VCF format

Some of the filters rely on BCFtools annotations and so when generating input VCFs
* run `bcftools mpileup` with `-a SP,DP` to provide depth and strand bias info
* run `bcftools call` using `-f GC -c`

Example for how to generate the input VCF file:

	bcftools mpileup -f <reference fasta> -b <bam list> -a SP,DP | bcftools call -f GQ -c - > unfiltered_sites.vcf

Example for how to pipe directly into snpCleaner without dumping an unfiltered VCF file:

	bcftools mpileup -b <bam list file> -f <reference file> -a SP,DP | bcftools call -f GQ -c | snpCleaner.pl -v -B <qc-passed sites bed file name> -p <bzipped failed sites file name> > <qc-passed sites VCF file name>

When generating input using BCFtools you may want to increase the max per-file depth with `bcftools mpileup -d` to ensure all reads are used. Also consider using 'bcftools call --skip-variants indels' since these sites will be ignored by snpCleaner anyways.

## Output Notes
* Indels are ignored.
* The file of sites that passed filters (-B) lists regions as CHR:START-END and single sites as CHR:POS
* FILTER field flags in the failed sites VCF (-p) indicate filters that the sites failed to pass and correspond to the option flags when running the program (e.g. a site that fails minimum site read depth and strand bias will be flaged 'd,S').
* The failed sites VCF file is bzipped.
