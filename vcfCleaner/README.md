vcfCleaner
==========

The vcfCleaner directory contains programs for VCF quality filtering and file format conversion useful for generating quality-controlled input files used by other programs. 

## Installation

	cd ./ngsQC/vcfCleaner
	make

## vcfCleaner

vcfCleaner takes a VCF file as input and quality filters it based on annotations from different VCF production pipelines. For an overview of running vcfCleaner on VCFs produced by different programs run vcfCleaner without arguments:

	./vcfCleaner

	vcfCleaner [VCF type] [arguments]

	VCF types:
	gatk        vcf produced by GATK HaplotypeCaller/GenotypeGVCFs workflow

### Filtering VCFs produced by the GATK HaplotypeCaller/GenotypeGVCFs workflow

`vcfCleaner gatk` filters VCFs with GATK annotations.
* Input VCF should be generated by GATK [GenotypeGVCFs](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php).
* If filtering monomorphic sites in addition to SNPs, run `GenotypeGVCFs` with the `-allSites` option. vcfCleaner does not work on gVCF blocks.
* Indels are ignored.

Run `vcfCleaner gatk` without arguments to display usage information:

	./vcfCleaner gatk

	This program is for filtering VCFs produced by GATK's HaplotypeCaller/GenotypeGVCFs workflow.
	Takes uncompressed and gzipped/bgzipped VCFs, however bgzipped reading can be buggy so use at your own risk.
	Indels are ignored and not included in files containing filtered sites.

	vcfCleaner gatk [arguments]

	-vcf                FILE    Input VCF file to filter ('-' for STDIN)
	-out                STRING  Output file name (prefix)
	-biallelic          0|1     0: keep multiallelic sites, 1: keep only biallelic sites, A [1]
	-allsites           0|1     0: process only SNPs, 1: process all sites (except indels) [0]
	-allsites_vcf       0|1     Output VCF contains QC-passed 0: SNPs only, 1: all sites [0]
	-maxcov             INT     Max total site depth, D [1000000]
	-mincov             INT     Min total site depth, C [2]
	-minind_cov         INT     Min individual coverage [1]
	-minind             INT     Min number of individuals covered by at least -minind_cov reads, U [1]
	-mingeno            INT     Min number of individuals with called genotype, G [1]
	-rms_mapq           FLOAT   Min RMS mapping quality, R [40]
	-mapqRankSum        FLOAT   Max absolute Wilcoxon rank sum test Z-score of alt vs. ref read map quality, M [12.5]
	-posbias            FLOAT   Max absolute Wilcoxon rank sum test Z-score of alt vs. ref read position bias, P [8]
	-strandbias         FLOAT   Max Phred-scaled Fisher's exact test p-value of strand bias, S [60]
	-baseqbias          FLOAT   Max absolute Wilcoxon rank sum test Z-score of alt vs ref base qualities, B [21]
	-qual               FLOAT   Min Phred-scaled quality score of ALT assertion, Q [30]
	-varq_depth         FLOAT   Min variant Phred-scaled confidence/quality by depth, V [2]
	-hetexcess          FLOAT   Max Phred-scaled p-value for exact test of excess heterozygosity, H [40]
	-varonly            0|1     The INFO AF for SNP-only output VCF must be in range 0=[0,1], or 1=(0,1) [1]
	-maf                FLOAT   Minor allele frequency lower bound for SNP-only VCF [0]
	-verbose            0|1|2   Level of warnings to issue: 0 = suppress all, 1 = site-level, 2 = individual-level [2]

	Other site QC fail flags:
	F, Unkown reference allele
	N, No data for any individuals

### Output

* VCF containing sites that pass quality controls
* A *.pos position file of sites that pass quality control with syntax [chr]:[start-stop] (1-based)
* 1-based position file of sites that fail quality control flagged for each failed filter


Fail flags:
* A: Biallelic
* C: Min total site depth
* D: Max total site depth
* U: Min number of individuals covered by at least `-minind_cov` reads
* G: Min number of individuals with called genotype
* R: RNS mapping quality
* M: Alternate vs. reference allele map quality
* P: Alternate vs reference allele position bias
* S: Strand bias
* B: Alternate vs. reference base quality bias
* Q: Min QUAL
* V: Min qualiity by depth
* H: Excess heterozygosity
* F: Unknown reference allele
* N: No data for any individuals


## expandrf

expandrf converts files in the vcfCleaner QC-passed region (*.pos) format to either a 1-based position file with syntax [chr position] or bed format. Run expandrf without arguments for usage info.

	./expandrf

	./expandrf [format: bed|pos] [input region file]

	formats:
	bed - chr region_start region_end (0-indexed)
	pos - chr position (1-indexed)

## vcf2bimbam

vcf2bimbam converts a VCF file with genotype likelihood information (`PL` FORMAT subfield) to bimbam format with posterior expected genotypes. Missing data is encoded as 'NA'. Run vcf2bimbam without arguments for usage information.

	./vcf2bimbam

	vcf2bimbam [input vcf] [genotype prior]

	priors:
	hwe   Hardy-Weinberg Equilibrium probability of the genotypes given the ALT allele frequency
	unif  Uniform probability for genotypes

	BIMBAM output:
	(1) SNP, (2) ALT, (3) REF, (4..N) Posterior expected genotype
	Missing genotypes encoded as 'NA'
