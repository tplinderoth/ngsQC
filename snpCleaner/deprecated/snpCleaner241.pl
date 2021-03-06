#!/usr/bin/perl

# SNPcleaner.pl
# Author: Tyler Linderoth, tylerp.linderoth@gmail.com
my $version = "2.4.1";

# TODO:
#    Add argument check
#    Output smaller BED file (join contiguous regions)

use strict;
use warnings;
use Getopt::Std;
use IO::Compress::Bzip2;

my %opts = ('?'=> 0, 'd'=>2, 'D'=>1000000, 'k'=>1, 'u'=>0, 'a'=>2, 'Q'=>10, 'S'=>1e-4, 'b'=>1e-100, 'f'=>0, 'e'=>1e-4, 'h'=>1e-4, 'F'=>0, 'H'=>0, 'L'=>0, 'A'=>undef,
	    'M'=>undef, 'B'=>undef, 'p'=>undef, 'r'=>undef, 'X'=>undef, 't'=>undef, 'g'=>undef, 'v'=>undef, '2'=>undef, 'q'=>2);

getopts('?2d:D:k:u:a:Q:S:b:f:e:h:F:H:L:A:M:B:p:r:X:q:tgv', \%opts);

die (qq/
#####################
# SNPcleaner v$version #
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
-q INT   ploidy [$opts{q}]
-d INT   minimum site read depth [$opts{d}]
-D INT   maximum site read depth [$opts{D}]
-k INT   minimum number of 'covered' individuals (requires -u) [$opts{k}] 
-u INT   minimum read depth for an individual to be considered 'covered' (requires -k) [$opts{u}]
-a INT   minimum number of alternate alleles for site [$opts{a}]
-Q INT   minimum RMS mapping quality for SNPS [$opts{Q}]
-S FLOAT min p-value for strand bias [$opts{S}]
-b FLOAT min p-value for base quality bias [$opts{b}]
-f FLOAT min p-value for map quality bias [$opts{f}]
-e FLOAT min p-value for end distance bias [$opts{e}]
-h FLOAT min p-value for exact test of HWE [$opts{h}]
-F FLOAT inbreeding coefficient value [$opts{F}]
-H FLOAT min p-value for exact test of excess of heterozygous [$opts{H}]
-L FLOAT min p-value for exact test of defect of heterozygous [$opts{L}]
-A FILE  ANCESTRAL fasta file (with FAI on same folder)
-M CHAR  mutation type(s) to remove (ex. 'CT_GA') (requires -A)
-B CHAR  name of BED file for kept SNP positions
-p CHAR  name of file to dump sites that failed filters (bziped)
-r FILE  list of contigs to exclude
-X FILE  BED file of exonic regions (sorted from lowest to highest contig)
-t       filter non-exonic sites (requires -X)
-g	 filter exons with SNPs out of HWE (requires -X)
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
 failed to pass and correspond to the option flags (e.g. 'S' indicates strand bias).
\n/) if ($opts{'?'} || (!$ARGV[0] && -t STDIN));

# Argument check
if($opts{t} || $opts{g}) {
    die(qq/option -X (exonic region BED file) is required for options -t and -g\n/) unless ($opts{X});
}
if($opts{M}) {
    die(qq/option -A (ANCESTRAL fasta file) is required for option -M\n/) unless ($opts{A});
}
if($opts{k} !~ m/^\d+$/g || $opts{u} !~ m/^\d+$/g) {
    die(qq/option -k and -u are inter-dependent (min number of -k INT individuals with less than -u INT coverage)\n/);
}

if ($opts{q} < 1) {
	die(qq/Invalid ploidy option -q: $opts{q} ... ploidy should be >= 1\n/);
}
else
{
	if ($opts{q} != 2 && ($opts{h} + $opts{H} + $opts{L}) > 0) {
		print STDERR "Warning: Hardy-Weinberg filtering only applicable for diploid data (-q 2)\n";
	}
}

my ($excont, @t, @seq, @buffer, %exons, %ancestral);
my ($n_sites, $ind_dp, $prev_format, @format, @ind_depth) = (0,0,"");

#open some necessary filehandles
open(BED, '>', $opts{B}) or die("ERROR: could not create BED file: $!") if($opts{B});
my $bz2 = new IO::Compress::Bzip2($opts{p}, 'Append' => 0) if($opts{p});


# Read file with list of excluded contigs
if($opts{r}) {
    open(RMV, '<', $opts{r}) or die("ERROR: could not open excluded CONTIGS file: $!");
    $excont = join("\t", <RMV>);
    close(RMV);
}

# Read FASTA index file
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

# Read exon/intron information file
if($opts{X}) {
    open(EXON, '<', $opts{X}) or die("ERROR: could not open EXON file: $!");
    while (<EXON>) {
	my @exon = split(/\s+/);
	push( @{$exons{$exon[0]}}, {'start' => $exon[1], 'end' => $exon[2], 'HWE' => 0} );
    }
    close(EXON);
}


my ($prev_contig, $cur_contig, $pos) = ('start','start',0);
$" = "\t"; #for formatting printed output

# check if input VCF is from @ARGV or STDIN
if (-t STDIN)
{
	die ("No input VCF\n") unless @ARGV;
}

# The core loop
while (<>) {
    my $violate = ''; # for flagging filter violations

    @t = split;

    # Skip (and print) header lines
    (print, next) if(m/^#/);

    # Update position
    $pos = $t[1];

    # If contig changed
    if($t[0] ne $cur_contig){
	$prev_contig = $cur_contig;
	$cur_contig = $t[0];
    }

    # Skip sites with unknown ref
    $violate .= 'N' if ($t[3] eq 'N');

    # Skip non-variable sites
    $violate .= 'v' if ($t[4] eq '.' && !$opts{v});

    # Skip Indels
    if (length($t[3]) > 1 || length($t[4]) > 1) {
	$violate .= 'I';
	next;
    }

    # Skip sites from excluded contigs
    $violate .= 'r' if ($opts{r} && $excont =~ /\b$cur_contig\b/);

    # Skip non-biallelic Sites

    if( scalar(split(/,/,$t[4])) > 1 && !$opts{2} ) {
	$violate .= '2';
    }

    # Read ANC base from FASTA
    my $anc_base;
    if($opts{A}) {
	my $n_lines = int($pos / $ancestral{$cur_contig}{'n_chars'} - 1e-6);
	my $extra_bytes_per_line = $ancestral{$cur_contig}{'n_bytes'} - $ancestral{$cur_contig}{'n_chars'};
	seek(FASTA, $ancestral{$cur_contig}{'start'} + $pos - 1 + $n_lines*$extra_bytes_per_line, 0);
	read(FASTA, $anc_base, 1);
	$anc_base = 'N' if($anc_base =~ m/[RYSWKMBDHV]/i);
	warn("WARNING: invalid ancestral base at ",$cur_contig,", pos ",$pos,": ",$anc_base,".\n") if($anc_base !~ m/[ACTGN]/i);

	# Skip non-biallelic sites (major and minor differ from ANCESTRAL)
	$violate .= "m($anc_base)" if($anc_base !~ m/\Q$t[3]\E|\Q$t[4]\E|N/i && $t[4] ne '.');
    }

    # Skip sites with specified mutation types
    if ($opts{M} && $t[3] =~ m/[ATCG]/i && $t[4] =~ m/[ATCG]/i) { # if site is variable
	my $refalt=$t[3].$t[4];
	my $altref=$t[4].$t[3];

	if ($opts{M} =~ /($refalt|$altref)/i) {
	    my @mutype = split("", $1);
	    $violate .= "M" if ($anc_base !~ m/$mutype[1]/i);
	}	
    }
    
    # Eveness across individuals coverage filter
    # check where in vcf FORMAT the DP ID is	
    if( $prev_format ne $t[8] ){
	$ind_dp = 0;
	$prev_format = $t[8];
	@format = split(":",$t[8]);
	foreach (@format) {
	    $ind_dp++ if($_ ne 'DP');
	    last if($_ eq 'DP');
	}
    }

    # count how many individuals have coverage >= $opts{u}
    my $covcount = 0;
    my @genoinfo = @t[9 .. $#t];
    if( $ind_dp <= $#format ) { # if DP is missing in vcf skip even coverage filter
	for(my $i=0; $i <= $#genoinfo; $i++) {
	    my @ind_info = split(":", $genoinfo[$i]);
	    $covcount++ if ($ind_info[$ind_dp] >= $opts{u});
	    $ind_depth[$i] += $ind_info[$ind_dp];
	}
	$violate .= 'k' if ($covcount < $opts{k});	
    }else {
	warn("WARNING: no individual depth information at ",$cur_contig,", pos ",$pos,". Check for \"samtools mpileup -D\" option.");
    }

    # Site coverage 
    my ($dp, $dp_alt) = (0,0);
    if ($t[7] =~ m/DP4=(\d+),(\d+),(\d+),(\d+)/i) {
	$dp = $1 + $2 + $3 + $4;
	$dp_alt = $3 + $4;
    }
    $violate .= 'd' if ($dp < $opts{d});
    $violate .= 'D' if ($dp > $opts{D});
    $violate .= 'a' if ($dp > 0 && $dp_alt < $opts{a});

    # Root-mean-square mapping quality of covering reads
    my $mq = $1 if ($t[7] =~ m/MQ=(\d+)/i);
    $violate .= 'Q' if ($mq && $mq < $opts{Q});

    # Strand, baseQ, mapQ, and tail distance bias
    my ($strand, $baseqb, $mapqb, $tail_dist);
    if ($t[7] =~ m/PV4=([^,]+),([^,]+),([^,]+),([^,;\t]+)/) {
	$strand = $1;
	$baseqb = $2;
	$mapqb = $3;
	$tail_dist = $4;
    }
    $violate .= 'S' if ($strand && $strand < $opts{S});
    $violate .= 'b' if ($baseqb && $baseqb < $opts{b});
    $violate .= 'f' if ($mapqb && $mapqb < $opts{f});
    $violate .= 'e' if ($tail_dist && $tail_dist < $opts{e});

    # Identify non-exonic regions
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
    
    # HWE exact test
    if ($opts{q} == 2 && $t[4] ne '.') {
	my %genocount = (homoa => 0, homob => 0, het => 0);
	foreach (@genoinfo) {
	    if (/0\/0:/) {
		$genocount{homoa}++;
	    } elsif (/1\/1:/) {
		$genocount{homob}++;
	    } elsif (/0\/1:|1\/0:/) {
		$genocount{het}++;
	    }
	}

	my ($pHWE, $pHI, $pLOW) = hwe_exact($genocount{het},$genocount{homoa},$genocount{homob},$opts{F});
	die(qq/Genotype counts less than 0\n/) if $pHWE == -1;
	if ($pHWE < $opts{h}) {
	    $violate .= "h(p=$pHWE)";
	    $exons{$cur_contig}[$exon_id]{'HWE'} = 1 if($exon_id > -1);
	}
	if ($pHI < $opts{H}) {
	    $violate .= "H(p=$pHI)";
	    $exons{$cur_contig}[$exon_id]{'HWE'} = 1 if($exon_id > -1);
	}
	if ($pLOW < $opts{L}) {
	    $violate .= "L(p=$pLOW)";
	    $exons{$cur_contig}[$exon_id]{'HWE'} = 1 if($exon_id > -1);
	}
    }

    # remove exons with SNPs out of HWE
    unless( $#buffer < 0 || 
	    ($exon_id >= 0 &&
	     $cur_contig eq $buffer[0]{'contig'} && 
	     $exon_id == $buffer[0]{'exon_id'}) ) {

	print_buffer($opts{g}, \@buffer, \%exons);
	undef @buffer;
    }
    
    push( @buffer, {'contig' => $cur_contig, 'exon_id' => $exon_id, 'pos' => $pos, 'violate' => $violate, 'vcf' => join("\t",@t)} );
    $n_sites++;
}
# Force printing of last entry
print_buffer($opts{g}, \@buffer, \%exons);

print(STDERR $n_sites." sites processed!\n");
# Print per-individual depth
for(my $i=0; $i <= $#ind_depth; $i++) {
    print(STDERR "Ind ".($i+1)." depth:\t".($ind_depth[$i]/$n_sites)."\n");
}


close FASTA if($opts{A});
close BED if($opts{B});
$bz2->close if($opts{p});


exit(0);


#################
### Functions ###
#################
sub print_buffer {
    my ($opts_g, $buffer, $exons) = @_;

    foreach my $g (@$buffer) {
	$g->{'violate'} .= 'g' if($opts_g && $g->{'exon_id'} >= 0 && $exons->{$g->{'contig'}}[$g->{'exon_id'}]{'HWE'} == 1);
	
	if( !$g->{'violate'} ){
	    print($g->{'vcf'}."\n");
	    print(BED $g->{'contig'}."\t".($g->{'pos'}-1)."\t".$g->{'pos'}."\n") if $opts{B};
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

    my ($obs_hets, $obs_homa, $obs_homb, $F) = @_;

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
