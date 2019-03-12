snpCleaner
==========

snpCleaner.pl takes a VCF format file as input and filters it based on various quality metrics. Indels are ignored. It outputs a VCF and associated bed-format file of sites that pass QC as well as a VCF-like file of sites that do not pass with flags indicating which filters each site failed. For a description of the  different filters run snpCleaner.pl without arguments.

	./snpCleaner.pl

	Usage: snpCleaner.pl [options] <in.vcf>
	or
	cat <in.vcf> | snpCleaner.pl [options]

	Options:
	-?	 help
	-2	 keep non-biallelic sites
	-q INT   ploidy [2]
	-d INT   minimum site read depth [2]
	-D INT   maximum site read depth [1000000]
	-k INT   minimum number of individuals covered by at least -u reads (requires -u) [1] 
	-u INT   minimum read depth for an individual (required by -k) [0]
	-a INT   minimum number of alternate alleles for site [2]
	-Q INT   minimum RMS mapping quality for SNPS [10]
	-S FLOAT min p-value for strand bias [0.0001]
	-b FLOAT min p-value for base quality bias [1e-100]
	-f FLOAT min p-value for map quality bias [0]
	-e FLOAT min p-value for end distance bias [0.0001]
	-h FLOAT min p-value for exact test of HWE [0.0001]
	-F FLOAT inbreeding coefficient value [0]
	-H FLOAT min p-value for exact test of heterozygote excess [0]
	-L FLOAT min p-value for exact test of heterozygote deficit [0]
	-A FILE  ANCESTRAL fasta file (with *.fai index in same directory)
	-M CHAR  mutation type(s) to remove (e.g. 'CT_GA' removes C->T SNPs) (requires -A)
	-B CHAR  name of bed-format file of sites that pass filters
	-p CHAR  name of file to dump sites that failed filters (bziped)
	-r FILE  list of contigs to exclude
	-X FILE  bed file of exonic regions (sorted from lowest to highest contig)
	-t       filter non-exonic sites (requires -X)
	-g	 filter exons with SNPs out of Hardy-Weinberg Equilibrium (requires -X)
	-v	 process nonvariants

	Generating input VCF with BCFtools:
	-run 'bcftools mpileup' with '-a SP,DP' to provide depth and strand bias info
	-run 'bcftools call' using '-f GC -c'

	An example for how to generate the input VCF file:
	bcftools mpileup -f <reference fasta> -b <bam list> -a SP,DP | bcftools call -f GQ -c - > unfiltered_sites.vcf 

	Consider using 'bcftools call --skip-variants indels' since these sites will
	be ignored by snpCleaner anyways.

	Output Notes:
	-bed format file is zero-based.
	-Characters preceeding filtered sites (dumped with option -p) indicate filters that the sites
	 failed to pass and correspond to the option flags (e.g. 'dS' indicates min site depth and strand bias).

## Input VCF format

Some of the filters rely on BCFtools annotations and so when generating input VCFs
* run `bcftools mpileup` with `-a SP,DP` to provide depth and strand bias info
* run `bcftools call` using `-f GC -c`

Example for how to generate the input VCF file:
	bcftools mpileup -f <reference fasta> -b <bam list> -a SP,DP | bcftools call -f GQ -c - > unfiltered_sites.vcf

Example for how to pipe directly into snpCleaner without dumping an unfiltered VCF file
	bcftools mpileup -b <bam list file> -f <reference file> -a SP,DP | bcftools call -f GQ -c | snpCleaner.pl -v -B <qc-passed sites bed file name> -p <bzipped failed sites file name> > <qc-passed sites VCF file name>

When generating input using BCFtools you may want to increase the max per-file depth with `bcftools mpileup -d` to ensure all reads are used. Also consider using 'bcftools call --skip-variants indels' since these sites will be ignored by snpCleaner anyways.

## Output Notes
* Indels are ignored
* bed format file is zero-based.
* Characters preceeding filtered sites (dumped with option -p) indicate filters that the sites failed to pass and correspond to the option flags when running the program (e.g. a site that fails minimum site read depth and strand bias will be flaged 'dS')
* The failed sites file is bzipped
