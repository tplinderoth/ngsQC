#!/usr/bin/perl -w
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

    SNPcleaner.pl

=head1 SYNOPSIS

Usage: 
SNPcleaner.pl [OPTIONS] <infile.vcf>
or
cat <infile.vcf> | SNPcleaner.pl [OPTIONS]

OPTIONS:

--help|-?	this help screen

input files:

--pop_info|-P	FILE	input file with population information (each line should list \'sample_name\tpop_ID\') for HWE filtering
--anc|-A	FILE	ancestral-state fasta file (with FAI in same directory)
--pileup|-G	FILE	pileup format file (all sites and individuals in input vcf must be in pileup)
--exons|-X	FILE	BED format file of exonic regions (sorted from lowest to highest numbered contig)

coverage filters:

--minDepth|-d	INT	minimum site read depth [$opts{d}]
--maxDepth|-D	INT	maximum site read depth [$opts{D}]
--minIndiv|-k	INT	minimum number of individuals with at least [-u INT]X coverage (requires SNPcleaner -u and mpileup -D) [$opts{k}] 
--minIndiv_cov|-u	INT	minimum individual coverage threshold used for -k (requires SNPcleaner -k and mpileup -D) [$opts{u}]
--minalt|-a	INT	minimum number of alternate alleles per site [$opts{a}]

bias and other quality-aspect filters:

--minRMSmap|-Q	INT	minimum RMS mapping quality for SNPs [$opts{Q}]
--mapqual|-f	FLOAT	min p-value for map quality bias [$opts{f}]
--strand_ind|-S	FLOAT	min p-value for strand bias from combining p-values across individuals [$opts{S}]
--strand_site|-s	FLOAT	min p-value for strand bias determined from read counts summed over all individuals at the site [$opts{s}]
--allele_bias|-T	FLOAT	min p-value for allele bias in potential heterozygotes [$opts{T}]
--hetero_llh|-R	FLOAT	skip called homozygotes with heterozygote likelihood less than -R FLOAT for allele bias filter [$opts{R}] 	
--basequal|-b	FLOAT	min p-value for base quality bias [$opts{b}]
--endbias|-e	FLOAT	min p-value for biased distance of alternate bases from ends of reads (indication of misalignment) [$opts{e}]

Hardy-Weinberg equilibrium filters:

--hwe|-h	FLOAT	min p-value for exact test of HWE (two-tailed) [$opts{h}]
--hetero_excess|-H	FLOAT	min p-value for exact test of excess heterozygotes [$opts{H}]
--hetero_deficit|-L	FLOAT	min p-value for exact test of deficient heterozygotes [$opts{L}]
--inbreed_coef|-F	FLOAT	inbreeding coefficient value [$opts{F}]
--rmv_nonHWEexons|-g	filter-out exons containing at least one SNP out of HWE (requires -X)

mutation type filters:

--mutation_rmv|-M	STR	mutation type(s) to remove (ex: '-M CT_GA' means remove C<=>T and G<=>A)
--one_dir|-w	remove mutations defined by -M in one direction (ex: '-M CT_GA' means remove C=>T and G=>A) (requires -A)
--alt_excess|-J	FLOAT	min p-value for excess substitutions defined by -M in sites called as nonvariable (requires -A if -w) [$opts{J}] 
--error|-E	FLOAT	sequencing error rate (required by -J for filtering mutation types) [$opts{E}]

general filters:

--nonvar|-v	process nonvariant sites (in addition to varients)
--keep_nonbinary|-2	keep non-biallelic sites
--exclude_contigs|-r	FILE	list of contigs\/chromosomes to exclude (each line should list a contig name exactly as it appears in the input vcf file)
--rmv_nonexonic|-t	filter-out non-exonic sites (requires -X)

output:

--bed|-B	FILE	name of dumped BED format file for sites that pass all filters
--failed_sites|-p	FILE	name of dumped file containing sites that failed at least one filter (bziped)
--vcfout|-o	FILE	name of dumped vcf file containing sites that passed filters
--ind_depth|-I	FILE	dumped file with individual mean depth

=head1 DESCRIPTION

    This script will read a VCF file and filter SNPs based on a set of rules. It is similar to SAMTOOLS 
    "varfilter.pl" script but extends some of the filters.
    NOTE: This script is only for SNP filtering and will ignore INDELs!

=head1 AUTHOR

    Tyler Linderoth - tylerp.linderoth _at_ gmail _dot_ com
    Filipe G. Vieira - fgarrettvieira _at_ gmail _dot_ com
    Matteo Fumagalli - 

=head1 CONTRIBUTORS

    Additional contributors names and emails here

=cut


# Let the code begin...

use strict;
use warnings;
use Getopt::Long;
use IO::Compress::Bzip2;
use Statistics::Distributions; # from Michael Kospach http://search.cpan.org/~mikek/Statistics-Distributions-1.02/Distributions.pm

my $version = "2.4.2";

Getopt::Long::Configure(qw{no_auto_abbrev no_ignore_case_always});

my %opts = ('?'=>0, 
'2'=>undef, 
'd'=>2, 
'D'=>1000000, 
'k'=>1, 
'u'=>0, 
'a'=>0, 
'Q'=>10, 
'S'=>1e-4, 
's'=>1e-4, 
'b'=>1e-100, 
'f'=>0, 
'e'=>1e-4, 
'h'=>0, 
'E'=>0.01, 
'J'=>1e-6, 
'F'=>0, 
'H'=>0, 
'L'=>0, 
'T'=>0.001, 
'R'=>0.05, 
'A'=>undef, 
'M'=>undef, 
'B'=>undef, 
'p'=>undef, 
'r'=>undef, 
'X'=>undef, 
't'=>undef, 
'o'=>undef, 
'g'=>undef, 
'G'=>undef, 
'v'=>undef, 
'w'=>undef, 
'I'=>undef, 
'P'=>"output");

#getopts('?2d:D:k:u:a:Q:S:s:b:f:e:h:E:F:G:H:J:L:A:M:B:p:r:X:I:Z:P:T:R:tgvw', \%opts);

GetOptions('help|?!'=> \$opts{'?'}, 
'keep_nonbinary|2!' => \$opts{2},
'rmv_nonexonic|t!' => \$opts{t},
'rmv_nonHWEexons|g!' => \$opts{g},
'nonvar|v!' => \$opts{v},
'one_dir|w!' => \$opts{w},
'minDepth|d=i' => \$opts{d},
'maxDepth|D=i' => \$opts{D},
'minIndiv|k=i' => \$opts{k},
'minIndiv_cov|u=i' => \$opts{u},
'minalt|a=i' => \$opts{a},
'minRMSmap|Q=f' => \$opts{Q},
'strand_ind|S=f' => \$opts{S},
'strand_site|s=f' => \$opts{s},
'basequal|b=f' => \$opts{b},
'mapqual|f=f' => \$opts{f},
'endbias|e=f' => \$opts{e},
'hwe|h=f' => \$opts{h},
'hetero_excess|H=f' => \$opts{H},
'hetero_deficit|L=f' => \$opts{L},
'allele_bias|T=f' => \$opts{T},
'inbreed_coef|F=f' => \$opts{F},
'anc|A=s' => \$opts{A},
'mutation_rmv|M=s' => \$opts{M},
'alt_excess|J=f' => \$opts{J},
'error|E=f' => \$opts{E},
'bed|B=s' => \$opts{B},
'failed_sites|p=s' => \$opts{p},
'exclude_contigs|r=s' => \$opts{r},
'exons|X=s' => \$opts{X},
'hetero_llh|R=f' => \$opts{R},
'pileup|G=s' => \$opts{G},
'ind_depth|I=s' => \$opts{I},
'pop_info|P=s' => \$opts{P},
'vcfout|o=s' => \$opts{o}
);

die (qq/
#####################
# SNPcleaner v$version #
#####################

Usage: 
SNPcleaner.pl [OPTIONS] <infile.vcf>
or
cat <infile.vcf> | SNPcleaner.pl [OPTIONS]

############ OPTIONS ############

--help|-?	this help screen\n
input files:\n
--pop_info|-P	FILE	input file with population information (each line should list \'sample_name\tpop_ID\') for HWE filtering
--anc|-A	FILE	ancestral-state fasta file (with FAI in same directory)
--pileup|-G	FILE	pileup format file (all sites and individuals in input vcf must be in pileup)
--exons|-X	FILE 	BED format file of exonic regions (sorted from lowest to highest numbered contig)\n
coverage filters:\n
--minDepth|-d	INT	minimum site read depth [$opts{d}]
--maxDepth|-D	INT	maximum site read depth [$opts{D}]
--minIndiv|-k	INT	minimum number of individuals with at least [-u INT]X coverage (requires SNPcleaner -u and mpileup -D) [$opts{k}] 
--minIndiv_cov|-u	INT	minimum individual coverage threshold used for -k (requires SNPcleaner -k and mpileup -D) [$opts{u}]
--minalt|-a	INT	minimum number of alternate alleles per site [$opts{a}]\n
bias and other quality-aspect filters:\n
--minRMSmap|-Q	INT	minimum RMS mapping quality for SNPs [$opts{Q}]
--mapqual|-f	FLOAT	min p-value for map quality bias [$opts{f}]
--strand_ind|-S	FLOAT	min p-value for strand bias from combining p-values across individuals [$opts{S}]
--strand_site|-s	FLOAT	min p-value for strand bias determined from read counts summed over all individuals at the site [$opts{s}]
--allele_bias|-T	FLOAT	min p-value for allele bias in potential heterozygotes [$opts{T}]
--hetero_llh|-R	FLOAT	skip called homozygotes with heterozygote likelihood less than -R FLOAT for allele bias filter [$opts{R}] 	
--basequal|-b	FLOAT	min p-value for base quality bias [$opts{b}]
--endbias|-e	FLOAT	min p-value for biased distance of alternate bases from ends of reads (indication of misalignment) [$opts{e}]\n
Hardy-Weinberg equilibrium filters:\n
--hwe|-h		FLOAT	min p-value for exact test of HWE (two-tailed) [$opts{h}]
--hetero_excess|-H	FLOAT	min p-value for exact test of excess heterozygotes [$opts{H}]
--hetero_deficit|-L	FLOAT	min p-value for exact test of deficient heterozygotes [$opts{L}]
--inbreed_coef|-F	FLOAT	inbreeding coefficient value [$opts{F}]
--rmv_nonHWEexons|-g	filter-out exons containing at least one SNP out of HWE (requires -X)\n
mutation type filters:\n
--mutation_rmv|-M	STR	mutation type(s) to remove (ex: '-M CT_GA' means remove C<=>T and G<=>A)
--one_dir|-w	remove mutations defined by -M in one direction (ex: '-M CT_GA' means remove C=>T and G=>A) (requires -A)
--alt_excess|-J	FLOAT	min p-value for excess substitutions defined by -M in sites called as nonvariable (requires -A if -w) [$opts{J}] 
--error|-E	FLOAT	sequencing error rate (required by -J for filtering mutation types) [$opts{E}]\n
general filters:\n
--nonvar|-v	process nonvariant sites (in addition to varients)
--keep_nonbinary|-2	keep non-biallelic sites
--exclude_contigs|-r	FILE	list of contigs\/chromosomes to exclude (each line should list a contig name exactly as it appears in the input vcf file)
--rmv_nonexonic|-t	filter-out non-exonic sites (requires -X)\n
output:\n
--bed|-B	FILE	name of dumped BED format file for sites that pass all filters
--failed_sites|-p	FILE	name of dumped file containing sites that failed at least one filter (bziped)
--vcfout|-o	FILE	name of dumped vcf file containing sites that passed filters
--ind_depth|-I	FILE	dumped file with individual mean depth

#################################

Notes:\n
Some of the filters rely on annotations generated by SAMtools\/BCFtools.
To use the eveness-of-coverage filters (options -k and -u), -D must be used with samtools mpileup.
If option -s or -S is set to 0, that particular strand bias filter is not performed.
It's recomended to use mpileup -I to ignore indels.
Characters in front of filtered sites (dumped with option -p) indicate filters that the site failed to pass.
\n/) if($opts{'?'} || (!$ARGV[0] && -t STDIN));

my $time_specs = localtime;

# Argument check
if($opts{t} || $opts{g}) {
    die(qq/option -X (exonic region BED file) is required for options -t and -g\n/) unless ($opts{X});
}
if($opts{M}) {
    die(qq/option -A (ANCESTRAL fasta file) is required for option -M if -w is activated\n/) if ($opts{w} && !$opts{A});
}
if($opts{k} !~ m/^\d+$/g || $opts{u} !~ m/^\d+$/g) {
    die(qq/option -k and -u are inter-dependent (min number of -k INT individuals with less than -u INT coverage)\n/);
}



my ($excont, @t, @seq, @snp_buffer, %exons, %ancestral, @ind_pop, @ind_depth, %flag_pos, %het_index);
my ($n_sites, $prev_contig, $cur_contig, $pos, $forward_reads, $alt_reads, $total_reads) = (0, 'start', 'start', 0, 0, 0, 0);



#### Open some necessary filehandles

#### Open output vcf format file
open(OUTVCF, '>', $opts{o}) or die("ERROR: Could not open OUTVCF file: $!");

#### Open and set up filter log file
open(LOG, '>', "$opts{o}.log");

{	my $ofh = select LOG;
	$| = 1;
	print LOG "SNPcleaner version $version\n\ninput commands\n";
	my $commands = '';
	foreach (keys(%opts)) {
		$commands .= " -$_ $opts{$_}" if $opts{$_};
	}
	$commands =~ s/^\s//;
	print LOG "$commands\n\nSTART TIME: $time_specs\n";
	select $ofh;
}	

#### Open output bed file

open(BED, '>', $opts{B}) or die("ERROR: could not create BED file: $!") if($opts{B});
my $bz2 = new IO::Compress::Bzip2($opts{p}, 'Append' => 0) if($opts{p});



#### Read file with list of excluded contigs
if($opts{r}) {
    open(RMV, '<', $opts{r}) or die("ERROR: could not open excluded CONTIGS file: $!");
    $excont = join("\t", <RMV>);
    close(RMV);
}



#### Read FASTA index file
if($opts{A}) {
    open(FASTA, '<', $opts{A}) or die("ERROR: could not open FASTA file: $!");
    binmode(FASTA);
    
    open(FAI, '<', $opts{A}.".fai") or die("ERROR: could not open FAI file: $!");
    while(<FAI>){
	chomp;
	my @fai = split(/\t/);
	$ancestral{$fai[0]}{'length'} = $fai[1];
	$ancestral{$fai[0]}{'start'} = $fai[2];
	$ancestral{$fai[0]}{'n_chars'} = $fai[3];
	$ancestral{$fai[0]}{'n_bytes'} = $fai[4];
    }
    close(FAI);
}



#### Read exon/intron information file
if($opts{X}) {
    open(EXON, '<', $opts{X}) or die("ERROR: could not open EXON file: $!");
    while (<EXON>) {
	my @exon = split(/\s+/);
	push( @{$exons{$exon[0]}}, {'start' => $exon[1], 'end' => $exon[2], 'HWE' => 0} );
    }
    close(EXON);
}


#### Open pileup file and create index -- current indexing method is memory intensive
my %pileup_index;
open(PILEUP, '<', $opts{G}) or die("ERROR: Could not open PILEUP file: $!");
while (<PILEUP>) {
	my $position = $1 if $_ =~ /^(\S+\s+\d+)\s/;
	$position =~ s/\s+/\t/; # make sure hash key is standardized 
	$pileup_index{$position} = tell(PILEUP) - length($_);
}


#### Open output file with individual mean depth
if($opts{I}) {
    open (ICOV, '>', $opts{I}) or die("ERROR: could not open ICOV file: $!");
}


my $npop=1;
my $output_pop_HWE_distrib = (($opts{h} == 0 && $opts{H} == 0 && $opts{L} == 0) ? (1) : (0));
my @pops = (0);
my %ind_pop;
if($opts{P} ne "output") {
    ### Open input file with population definition
    open (POPIN, '<', $opts{P}) or die("ERROR: could not open population file: $!");
    %ind_pop = map { chomp; split(/\t/) } <POPIN>;
    close(POPIN);
    
    ### Count number of populations
    #@pops = sort keys %{{ map {$_ => 1} values(%ind_pop) }};
    my %found;
    @pops = grep { ! $found{$_}++ } values(%ind_pop);
    $npop = $#pops+1;
    print "Found ", $npop, " subpopulations.", "\n";
}

if($output_pop_HWE_distrib) {
    ### Open output file with per site HWE
    open (HWEOUT, '>', $opts{P}.".hwe_distrib") or die("ERROR: could not open HWEOUT file: $!");
    ### Prepare header
    my @buffer;
    foreach my $pop (@pops) {
	next if($pop eq "N/A");
	push(@buffer, "pop".$pop."_pHWE", "pop".$pop."_pHWE_H", "pop".$pop."_pHWE_L");
    }
    print (HWEOUT join("\t", @buffer, "min_pHWE", "min_pHWE_H", "min_pHWE_L")."\n");
}

$" = "\t"; #for formatting printed output

# test for VCF data from STDIN or @ARGV

if (-t STDIN) {
	die ("ERROR: No input VCF\n") unless @ARGV;
}

# The core loop
while (<>) {

    my $violate = ''; # for flagging filter violations
    
    @t = split;
    
    ##### Get indiv information
    @ind_pop = (0)x($#t-9+1);
    if(m/^#CHROM/i) {
	foreach my $ind (@t[9 .. $#t]) {
	    push(@ind_pop, $ind_pop{$ind});
	}
    }

    ##### Print (and skip) header lines
    if ($_ =~ m/^#/) {
    	print OUTVCF $_;
    	next;
    }
    
	#### get pileup line
    
    seek PILEUP, $pileup_index{"$t[0]\t$t[1]"}, 0;
    my $pileup_line = <PILEUP>;

    ##### Update position
    $pos = $t[1];

    ##### If contig changed
    if($t[0] ne $cur_contig){
	$prev_contig = $cur_contig;
	$cur_contig = $t[0];
    }

    ##### Skip sites with unknown ref
    $violate .= 'N' if ($t[3] eq 'N');

    ##### Skip non-variable sites
    $violate .= 'v' if ($t[4] eq '.' && !$opts{v});

    ##### Skip Indels
    my $indel_size = 1;
    $indel_size *= map { length($_) } split(/,/,$t[3]);
    if (length($t[3]) > 1 || $indel_size > 1) {
	$violate .= 'I';
	next;
    }
    
    ##### Skip non-biallelic Sites
    
    my $alternate_num = $t[4] eq '.' ? 0 : scalar(split(/,/,$t[4]));
    
    if( $alternate_num > 1 && !$opts{2} ) {
	$violate .= '2';
    }

    ##### Skip sites from excluded contigs
    $violate .= 'r' if ($opts{r} && $excont =~ /\b$cur_contig\b/);
    
    # get read type counts
    
    if ($t[7] =~ m/DP4=(\d+),(\d+),(\d+),(\d+)/i) {
		$forward_reads = $1 + $3;
		$total_reads = $1 + $2 + $3 +$4;
		$alt_reads = $3 + $4;
    }

    ##### Read ANC base from FASTA
    my $anc_base = '';
    if($opts{A}) {
	my $n_lines = int($pos / $ancestral{$cur_contig}{'n_chars'} - 1e-6);
	my $extra_bytes_per_line = $ancestral{$cur_contig}{'n_bytes'} - $ancestral{$cur_contig}{'n_chars'};
	seek(FASTA, $ancestral{$cur_contig}{'start'} + $pos - 1 + $n_lines*$extra_bytes_per_line, 0);
	read(FASTA, $anc_base, 1);
	$anc_base = 'N' if($anc_base =~ m/[RYSWKMBDHV]/i);
	warn("WARNING: invalid ancestral base at ",$cur_contig,", pos ",$pos,": ",$anc_base,".\n") if($anc_base !~ m/[ACTGN]/i);

	# Skip non-biallelic sites (major and minor differ from ANCESTRAL)
	$violate .= "m($anc_base)" if($anc_base !~ m/\Q$t[3]\E|\Q$t[4]\E|N/i && $t[4] ne '.' && !$opts{2});
    }

    ##### Skip sites with specified mutation types
    if ($opts{M} && $t[3] =~ m/[ATCG]/i && $t[4] =~ m/[ATCG]/i) { # if site is variable
	my $refalt=$t[3].$t[4];
	my $altref=$t[4].$t[3];

	if ($opts{M} =~ /($refalt|$altref)/i) {
	    my @mutype = split("", $1);
	    if ($opts{w}) {
	    	$violate .= "M" if ($anc_base !~ m/$mutype[1]/i);
	    } else {
	    	$violate .= "M";
	    }
	}	
    }
    
    # secondary mutation type (DNA damage) filter
    
    if ($alt_reads > 0 && $alternate_num == 0) {
		if ($opts{M} =~ /$t[3]/i) {
			my @m = $opts{M} =~ /([agct]$t[3]|$t[3]\s{0}[acgt])/gi;
			my ($alt_allele, $mutation_pval) = DNAdamage($pileup_line, $opts{E}, $t[3], \@m);
			if ($alt_allele ne 'NA' ) {
				if ($opts{w}) {
					my $mutation_rmv = $1 if $opts{M} =~ /($alt_allele$t[3]|$t[3]$alt_allele)/;
					if ($mutation_rmv =~ /\w(\w)/) {
						$violate .= 'J' if ($1 !~ /$anc_base/i && $mutation_pval < $opts{J});
					}		
				} else {
					$violate .= 'J' if ($mutation_pval < $opts{J});
				}
			}
		}
	}
        
	# get individual info flag indexes
    
    if (!exists($flag_pos{$t[8]})) {
    	my @format = split(':', $t[8]);
    	map { $flag_pos{$t[8]}{$format[$_]} = $_ } 0 .. $#format; 
    }
    
    # get index of heterozygote likelihoods in individual information array
    if ($alternate_num > 0) {
    	if (!exists($het_index{$alternate_num})) {
    		my @homo_index;
    		push @homo_index, 0;
    		my ($offset, $homo_pos) = (2, 0);
    		for (my $i = 2; $i <= $alternate_num + 1; $i++) {
    			$homo_pos += $offset;
    			$offset++;
    			push @homo_index, $homo_pos;
    		}
    		my $h = shift @homo_index;
    		my $last_index = $homo_index[-1];
    		for (my $j = 0; $j <= $last_index; $j++) {
    			if ($j != $h) {
    				push @{$het_index{$alternate_num}}, $j;
    			} else {
    				$h = shift @homo_index;
   	 			}
   	 		}
   	 	}
    }
    
    # collect individual information
    
    my @genoinfo = @t[9 .. $#t];
    my $ind_genotypes = {};
    
    ##### count how many individuals have coverage >= $opts{u} and collect genotypic information
    
    my $covcount = 0;
    if( exists($flag_pos{$t[8]}{DP}) ) { # if DP is missing in vcf skip even coverage filter
		for(my $i=0; $i <= $#genoinfo; $i++) {
		    my @ind_info = split(":", $genoinfo[$i]);
		  	if ($alternate_num > 0) {
		  		if ($ind_info[$flag_pos{$t[8]}{GT}] eq './.') {
		    		push @{$$ind_genotypes{$i}}, 'M'; # missing genotype
		    		next;
		    	} else {
		    		push @{$$ind_genotypes{$i}}, $ind_info[$flag_pos{$t[8]}{GT}];
		    		my @geno_llh = split(/,/, $ind_info[$flag_pos{$t[8]}{PL}]);
		    		# get greatest heterozygote likelihood - for base type bias filter
		    		my $max_llh = $geno_llh[${$het_index{$alternate_num}}[0]];
		    		map { $max_llh = $_ if $_ < $max_llh } @geno_llh[@{$het_index{$alternate_num}}];
		    		push @{$$ind_genotypes{$i}}, (10 ** (-$max_llh / 10));
				}
		  	}
		    $covcount++ if ($ind_info[$flag_pos{$t[8]}{DP}] >= $opts{u});
		    $ind_depth[$i] += $ind_info[$flag_pos{$t[8]}{DP}];
		}
	$violate .= 'k' if ($covcount < $opts{k});	
    } else {
		die("ERROR: no individual depth information at $cur_contig pos $pos. Check for \"samtools mpileup -D\" option.");
    }

    ##### Site coverage 

    $violate .= 'd' if ($total_reads < $opts{d});
    $violate .= 'D' if ($total_reads > $opts{D});
    $violate .= 'a' if ($total_reads > 0 && $alt_reads < $opts{a});

    ##### allele type bias  
    
    if ($alternate_num > 0) {
    	my $allele_state_bias = type_test($pileup_line, $ind_genotypes, $opts{R});
    	$violate .= 'T' if ($allele_state_bias >= 0 && $allele_state_bias < $opts{T});
    }
    
    ##### Root-mean-square mapping quality of covering reads
    my $mq = $1 if ($t[7] =~ m/MQ=(\d+)/i);
    $violate .= 'Q' if ($mq && $mq < $opts{Q});

    ##### baseQ, mapQ, and tail distance bias
    my ($strand, $baseqb, $mapqb, $tail_dist);
    if ($t[7] =~ m/PV4=([^,]+),([^,]+),([^,]+),([^,;\t]+)/) {
	$strand = $1;
	$baseqb = $2;
	$mapqb = $3;
	$tail_dist = $4;
    }
    # $violate .= 'S' if ($strand && $strand < $opts{S}); # deprecated, preffer custom strand bias test
    $violate .= 'b' if ($baseqb && $baseqb < $opts{b});
    $violate .= 'f' if ($mapqb && $mapqb < $opts{f});
    $violate .= 'e' if ($tail_dist && $tail_dist < $opts{e});
	
	##### strand bias 
	
	# site-based strand bias by combining p-values across individuals
	
	if ($opts{S} > 0) {
		my $individual_strand_bias = individualStrandBias($pileup_line);
		$violate .= 'S' if ($individual_strand_bias >= 0 && $individual_strand_bias < $opts{S});
	}
	
	# site-based strand bias by considering reads counts for the site as a whole
	
	if ($opts{s} > 0) {
		my $site_strand_bias = binom_test($forward_reads, $total_reads, 0.5, 'two.sided');
		$violate .= 's' if ($site_strand_bias < $opts{s});
	}
	
    ##### Identify non-exonic regions
    my $exon_id = -1;
    if ($opts{X})  {
	my $n_exons = scalar(@{$exons{$cur_contig}});
	for($exon_id = 0; $exon_id < $n_exons; $exon_id++) {
	    last if ($pos >= $exons{$cur_contig}[$exon_id]{'start'} && 
		     $pos <= $exons{$cur_contig}[$exon_id]{'end'});
	}
	$exon_id = -1 if($exon_id >= $n_exons);
	$violate .= 't' if ( $opts{t} && $exon_id < 0 );
    }
    
    ##### HWE exact test per population
    if ($t[4] ne '.') {
	my %count_geno; # 3 values for each pop
	my ($min_pHWE, $min_pHI, $min_pLOW, @buffer) = (1,1,1);
	
        for(my $ind=0; $ind <= $#genoinfo; $ind++) {
	    if    ($genoinfo[$ind] =~ m/0\/0:/)       { $count_geno{$ind_pop[$ind]}[0]++; }
	    elsif ($genoinfo[$ind] =~ m/0\/1:|1\/0:/) { $count_geno{$ind_pop[$ind]}[1]++; }
	    elsif ($genoinfo[$ind] =~ m/1\/1:/)       { $count_geno{$ind_pop[$ind]}[2]++; }
        }
		
        foreach my $pop ( sort keys(%count_geno) ) {

	    next if ($pop eq "N/A");
	    my ($pHWE, $pHI, $pLOW) = hwe_exact($count_geno{$pop}[1], $count_geno{$pop}[0], $count_geno{$pop}[2], $opts{F});
	    die(qq/Genotype counts less than 0\n/) if $pHWE == -1;
	    
	    push(@buffer, $pHWE, $pHI, $pLOW);
	    $min_pHWE = (($min_pHWE < $pHWE) ? ($min_pHWE) : ($pHWE));
	    $min_pHI  = (($min_pHI  < $pHI)  ? ($min_pHI)  : ($pHI));
	    $min_pLOW = (($min_pLOW < $pLOW) ? ($min_pLOW) : ($pLOW));
	    
	    if (!$output_pop_HWE_distrib) {
		if ($pHWE < $opts{h}) {
		    $violate .= "h(p=$pHWE;$pop)";
		    $exons{$cur_contig}[$exon_id]{'HWE'} = 1 if($exon_id > -1);
		}
		if ($pHI < $opts{H}) {
		    $violate .= "H(p=$pHI;$pop)";
		    $exons{$cur_contig}[$exon_id]{'HWE'} = 1 if($exon_id > -1);
		}
		if ($pLOW < $opts{L}) {
		    $violate .= "L(p=$pLOW;$pop)";
		    $exons{$cur_contig}[$exon_id]{'HWE'} = 1 if($exon_id > -1);
		}
	    }
	}
	if ($output_pop_HWE_distrib && $violate eq '') {print(HWEOUT join("\t", @buffer, $min_pHWE, $min_pHI, $min_pLOW)."\n");  }
    }

    ##### Remove exons with SNPs out of HWE
    unless( $#snp_buffer < 0 || 
	    ($exon_id >= 0 &&
	     $cur_contig eq $snp_buffer[0]{'contig'} && 
	     $exon_id == $snp_buffer[0]{'exon_id'}) ) {

	print_buffer($opts{g}, \@snp_buffer, \%exons);
	undef @snp_buffer;
    }

    
    push( @snp_buffer, {'contig' => $cur_contig, 'exon_id' => $exon_id, 'pos' => $pos, 'violate' => $violate, 'vcf' => join("\t",@t)} );
    $n_sites++;
}
# Force printing of last entry
print_buffer($opts{g}, \@snp_buffer, \%exons);

print(STDERR $n_sites." sites processed!\n");
# Print per-individual depth
for(my $i=0; $i <= $#ind_depth; $i++) {
    print(STDERR "Ind ".($i+1)." depth:\t".($ind_depth[$i]/$n_sites)."\n");
    if($opts{I}) {
       print(ICOV "Ind ".($i+1)." depth:\t".($ind_depth[$i]/$n_sites)."\n");
    } 
}


close FASTA if($opts{A});
close BED if($opts{B});
$bz2->close if($opts{p});
close ICOV if($opts{I});
close HWEOUT if ($opts{P});
close PILEUP if ($opts{G});
close OUTVCF if ($opts{o});

$time_specs = localtime;
print LOG "END TIME: $time_specs\n";
close LOG;

exit(0);


#################
### Functions ###
#################

sub print_buffer {
    my ($opts_g, $buffer, $exons) = @_;

    foreach my $g (@$buffer) {
	$g->{'violate'} .= 'g' if($opts_g && $g->{'exon_id'} >= 0 && $exons->{$g->{'contig'}}[$g->{'exon_id'}]{'HWE'} == 1);
	
	if( !$g->{'violate'} ){
	    print(OUTVCF $g->{'vcf'}."\n") if $opts{o};
	    print(BED $g->{'contig'}."\t".$g->{'pos'}."\t".($g->{'pos'}+1)."\n") if $opts{B};
	}else{
	    $bz2->print($g->{'violate'}."\t".$g->{'vcf'}."\n") if($opts{p});
	}
    }
}



sub hwe_exact {

# Citation:
# Implements an exact SNP test of Hardy-Weinberg Equilibrium as described in Wigginton et al. 2005
# note that probabilities are calculated from the midpoint in order to take advantage of the recurrence
# relationships recognized in Guo and Thompson (1992) in the implementation of their MCMC sampler

# Funtion return 3 p-values for deviations from HWE: two-tail, one-tail (excess of heterozygotes) and one-tail (excess of homozygotes)

    my ($obs_hets, $obs_homa, $obs_homb, $F) = @_;

    $obs_hets = 0 if( !defined($obs_hets) );
    $obs_homa = 0 if( !defined($obs_homa) );
    $obs_homb = 0 if( !defined($obs_homb) );

    return(-1) if ($obs_hets < 0 || $obs_homa < 0 || $obs_homb < 0);

    my $obs_homr; #rare homozygote
    my $obs_homc; #commmon homozygote
    my $n = $obs_homa + $obs_homb + $obs_hets; # total number genotypes

    # define common and rare homozygotes
    if ($obs_homa > $obs_homb) {
	$obs_homc = $obs_homa;
	$obs_homr = $obs_homb;
    } elsif ($obs_homa < $obs_homb) {
	$obs_homc = $obs_homb;
	$obs_homr = $obs_homa;
    } elsif ($obs_homa == $obs_homb) { # need to check how matching number homos affects algorithm
	$obs_homc = $obs_homa;
	$obs_homr = $obs_homb;
    }

    my $rare = 2 * $obs_homr + $obs_hets; # number of minor alleles

    # theta for inbreeding HWE calculations
    my $pc = 1 - $rare/(2*$n);
    my $pr = 1 - $pc;
    my $pCC = $pc**2 + $pc*$pr*$F;
    my $pCR = 2*$pc*$pr - 2*$pc*$pr*$F;
    my $pRR = $pr**2 + $pc*$pr*$F;
    $pRR = 1e-6 if($pRR == 0);
    my $theta = ($pCR**2)/($pCC*$pRR);
    $theta = 1e-6 if($theta == 0);

    # initialize heterozygote probability array
    my @probs;
    for (my $i = 0; $i <= $rare; $i++) {
	$probs[$i] = 0.0;
    }

    # find midpoint of the minor allele count distribution
    my $mid = int($rare * (2 * $n - $rare) / (2 * $n));
    $mid = $mid + 1 if ( ($mid % 2) != ($rare % 2) ); # ensures number minor alleles and midpoint have parity

    my $curr_hets = $mid;
    my $curr_homr = ($rare - $mid) / 2;
    my $curr_homc = $n - $curr_hets - $curr_homr;

    $probs[$mid] = 1.0;
    my $sum = $probs[$mid];

    # calculate probabilities from midpoint down 
    for ($curr_hets = $mid; $curr_hets > 1; $curr_hets -= 2) {
	$probs[$curr_hets - 2] = $probs[$curr_hets] * $curr_hets * ($curr_hets - 1) / 
	    ($theta * ($curr_homr + 1) * ($curr_homc + 1));
	$sum += $probs[$curr_hets - 2];

	# 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
	$curr_homr++;
	$curr_homc++;
    }

    # calculate  probabilities from midpoint up
    $curr_hets = $mid;
    $curr_homr = ($rare - $mid) / 2;
    $curr_homc = $n - $curr_hets - $curr_homr;

    for ($curr_hets = $mid; $curr_hets <= $rare - 2; $curr_hets += 2) {
	$probs[$curr_hets + 2] = $probs[$curr_hets] * $theta * $curr_homr * $curr_homc /
	    (($curr_hets + 2) * ($curr_hets + 1));
	$sum += $probs[$curr_hets + 2];
	
	# add 2 heterozygotes for next interation -> subtract one rare, one common homozygote
	$curr_homr--;
	$curr_homc--;
    }
    for (my $i = 0; $i <= $rare; $i++) {
	$probs[$i] /= $sum;
    }

    # p-value calculation for hwe
    my $p_hwe = 0.0;
    for (my $i = 0; $i <= $rare; $i++) {
	next if ($probs[$i] > $probs[$obs_hets]);
	$p_hwe += $probs[$i];
    }
    $p_hwe = 1.0 if ($p_hwe > 1);

    # alternate p-value calculation for p_hi/p_low heterozygous
    my $p_hi = $probs[$obs_hets];
    for (my $i = $obs_hets + 1; $i <= $rare; $i++) {
	$p_hi += $probs[$i];
    }
    my $p_low = $probs[$obs_hets];
    for (my $i = $obs_hets - 1; $i >= 0; $i--) {
	$p_low += $probs[$i];
    }
#    my $p_hi_low;
#    if ($p_hi < $p_low) {
#	$p_hi_low = 2 * $p_hi;
#    } else {
#	$p_hi_low = 2 * $p_low
#    }

    return($p_hwe, $p_hi, $p_low);
} 

sub choose {
    my ($n, $k) = @_;
    my ($num, $j) = (1, 1);
   
    return 0 if $k > $n || $k < 0; # shouldn't need this
    $k = ($n - $k) if ($n - $k) < $k;
   
    while ($j <= $k) {
        $num *= $n--;
        $num /= $j++;
    }
    return $num;
}

sub binom_test {
    my ($k, $n, $p, $tail) = @_;
    my ($prob_low, $prob_hi, $i, $j) = (0, 0, 0, $k);
    
    # calculate probability in left tail
    
    map { $prob_low += choose($n, $_) * ($p ** $_) * ((1 - $p) ** ($n - $_)) } 0 .. $k;
   
   return ($prob_low) if $tail eq 'lower';
   
   # calculate probability in right tail 
   
    map { $prob_hi += choose($n, $_) * ($p ** $_) * ((1 - $p) ** ($n - $_)) } $k .. $n;
   
   return ($prob_hi) if $tail eq 'upper';
   
    # return two-tailed p-value
   
   if ($tail eq 'two.sided') {
    	return $prob_low < $prob_hi ? 2 * $prob_low : 2 * $prob_hi;
   }
  }


sub type_test {
	
	my $p = 0.5;
	chomp $_[0];
	my @pileup = split(/\s+/, $_[0]);
	my ($geno, $het_thresh) = @_[1..2];
	
	# check that reference base is A,C,T,G
	
	return (-1) if $pileup[2] !~ /[ACGT]/i;
	
	my $pval = 0;
	my $ind_field = 1;
	my $k = 0;
	
	# get allelic type counts for each individual
	
	my $individual = -1; # individual key in %ind_genotypes
	for (my $i = 3; $i <= $#pileup; $i++) {
		
		if ($ind_field == 1) {
			$individual++;
			
			# check if individual is called a homozygote and the likelihood of it being a heterozygote
			
			if (${$$geno{$individual}}[0] =~ /(\d+)[\/\|](\d+)/) {
				if ($1 == $2 && ${$$geno{$individual}}[1] <= $het_thresh) {
					$i += 2 unless $pileup[$i+1] =~ /^\d+/;
					next;
				}
			} elsif (${$$geno{$individual}}[0] eq 'M') {
				$i += 2 unless $pileup[$i+1] =~ /^\d+/;
				next;
			}
			
			# skip to next individual if read depth is zero
			
			if ($pileup[$i] == 0) {
				$i += 2 if $pileup[$i+1] eq '*';
				next;
				}
			
			$k++; # degrees of freedom for chi-square test
			
		} elsif ($ind_field == 2) {
			
			# first remove any symbol that isn't a ref or alternate base
			
			while ( ($pileup[$i] =~ /[+|-](\d+)[ACGTNacgtn]/) == 1) {
				$pileup[$i] =~ s/[+|-]$1[ACGTNacgtn]{$1}//;
			}
			$pileup[$i] =~ s/\^\S|[^ACGT\.,]//gi;
			
			# count alternative bases and pick the most frequent as alternate
			
			my $alt = '';
			my $alt_count = 0;
			foreach ('A', 'C', 'G', 'T') {
				 my $c = () = $pileup[$i] =~ /$_/ig;
				 if ($c > $alt_count) {
				 	$alt = $_;
				 	$alt_count = $c;
				 }	 
			}
			
			# count the number of reference base reads & get binomial p-value w/ respect to most common alternate allele
			$pileup[$i] =~ s/[^$alt\.,]//ig;	 
			my $depth = length($pileup[$i]);
			my $ref_counts = () = $pileup[$i] =~ /[\.,]/ig;
			$pval += log(binom_test($ref_counts, $depth, $p, 'two.sided')) if ($ref_counts != 0);
		}
		
		$ind_field++;
		$ind_field = 1 if $ind_field > 3;
	
	}
	
	# finish calculation for combining the p-values using Fisher's method & return smallest pval
	
	my $chisq = -2 * $pval;
	my $df = 2 * $k;
	return(-1) if $df < 1;
	my $chi_pvalue = Statistics::Distributions::chisqrprob($df, $chisq);
	
	return($chi_pvalue);
}

sub individualStrandBias {
	
	chomp $_[0];
	my @pileup = split(/\s+/, $_[0]);
	my $p = 0.5;
	my ($ind_field, $k, $pval) = (1, 0, 0);
	
		for (my $i = 3; $i <= $#pileup; $i++) {
		
		if ($ind_field == 1) {
			
			# skip to next individual if read depth is zero
			
			if ($pileup[$i] == 0) {
				$i += 2 if $pileup[$i+1] eq '*';
				next;
				}
			
			$k++; # degrees of freedom for chi-square test
			
		} elsif ($ind_field == 2) {
			
			# first remove any symbol that isn't a ref or alternate base
			
			while ( ($pileup[$i] =~ /[+|-](\d+)[ACGTNacgtn]/) == 1) {
				$pileup[$i] =~ s/[+|-]$1[ACGTNacgtn]{$1}//;
			}
			$pileup[$i] =~ s/\^\S|[^ACGT\.,]//gi;
			
			# count the number of forward reads and perform a binomial test on the individual
				 
			my $depth = length($pileup[$i]);
			my $forward_counts = () = $pileup[$i] =~ /[\.ACGT]/g;
			
			$pval += log(binom_test($forward_counts, $depth, $p, 'two.sided'));
		}
		
		$ind_field++;
		$ind_field = 1 if $ind_field > 3;
	
	}
	
	my $chisq = -2 * $pval;
	my $df = 2 * $k;
	return(-1) if $df < 1;
	my $chi_pvalue = Statistics::Distributions::chisqrprob($df, $chisq);
	return ($chi_pvalue);

}

sub DNAdamage {
	
	chomp $_[0];
	my @pileup = split(/\s+/, $_[0]);
	my ($error, $ref_allele, $mutations) = @_[1 .. 3];
	my (%pval, %df);
	foreach (@$mutations) {
		my $alt = $1 if $_ =~ /($ref_allele\w|\w$ref_allele)/;
		$alt =~ s/$ref_allele//;
		$pval{$alt} = 0 unless exists($pval{$alt});
		$df{$alt} = 0 unless exists($df{$alt});
	}
	my $ind_field = 1;
	
	for (my $i = 3; $i <= $#pileup; $i++) {
		
		if ($ind_field == 1) {
			
			# skip to next individual if read depth is zero
			
			if ($pileup[$i] == 0) {
				$i += 2 if $pileup[$i+1] eq '*';
				next;
				}
						
		} elsif ($ind_field == 2) {
			
			# first remove any symbol that isn't a ref or alternate base
			
			while ( ($pileup[$i] =~ /[+|-](\d+)[ACGTNacgtn]/) == 1) {
				$pileup[$i] =~ s/[+|-]$1[ACGTNacgtn]{$1}//;
			}
			$pileup[$i] =~ s/\^\S|[^ACGT\.,]//gi;
			
			# count the number of alternate reads and perform a binomial test on the individual
				 
			my $depth = length($pileup[$i]);
			foreach (keys(%pval)) {
				my $alt_counts = () = $pileup[$i] =~ /$_/ig;
				if ($alt_counts > 0) {
					$pval{$_} += log(binom_test($alt_counts, $depth, $error, 'upper'));
					$df{$_}++;
				}
			}
		}
		$ind_field++;
		$ind_field = 1 if $ind_field > 3;
	
	}
	
	my $minp = 1;
	my $nonref = 'NA';
	
	foreach (keys(%pval)) {
		my $chisq = -2 * $pval{$_};
		my $k = $df{$_} * 2;
		next if $k < 1;
		my $chi_pvalue = Statistics::Distributions::chisqrprob($k, $chisq);
		if ($chi_pvalue < $minp) {
			$minp = $chi_pvalue;
			$nonref = $_;
		}
	}
	
	
	return ($nonref, $minp);
	
}

