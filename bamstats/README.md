bamstats
========

## Description:

Extracts the following stats from a BAM file:

* Total depth
* Average base and mapping qualities
* Root mean square (RMS) base and mapping qualities
* Fraction of base quality zero and mapping quality zero reads
* Number of samples with data

## Compile:

Once in directory containing bamstats.cpp

`g++ -O3 -o bamstats bamstats.cpp`

## Usage:

```
bamstats [options] < BAM | -b BAM_LIST >

The BAM list supplied with -b should contain one BAM file per row.

Native options:
--minind INT                    Only print sites for which at least INT individuals have data. [0]
--qual_offset INT               Subtract INT from the raw Phred score values. [33]

SAMtools mpileup options (set to SAMtools defaults):
-A, --count-orphans             Keep anamolous read pairs.
-a                              Output all positions, including zero depth sites.
-aa                             Output all positions, including zero depth sites and unused reference sequences.
-B, --no-BAQ                    Disable base alignment quality.
-C, --adjust-MQ INT             Adjust mapping quality (0: disable).
-A, --count-orphans             Keep anamolous read pairs.
-d, --max-depth INT             Maximum per BAM file depth.
-E, --redo-BAQ                  Recalculate base alignment qualities.
-f, --fasta-ref FILE            Indexed reference sequence in fasta format.
-G, --exclude-RG FILE           Exclude read groups listed in FILE.
-r, --region STRING             Region to analyze.
-l, --positions FILE            Skip unlisted positions in BED region or "chr position" format.
-q, --min-MQ INT                Skip alignments with map quality less than INT.
-Q, --min-BQ INT                Skip alignments with base quality less than INT.
-A, --count-orphans             Keep anamolous read pairs.
-R, --ignore-RG                 Ignore RG tags.
--rf, --incl-flags STRING|INT   Only keep reads with any of the mask bits set (STRING is comma delimited).
--ff, --excl-flags STRING|INT   Skip reads with any of the mask bits set (STRING is comma delimited).
-x, --ignore-overlaps, --disable-overlap-removal  Disable reads pair overlap detection and removal.
-X, --customized-index FILE     Use custome index files.

Notes:
*Options -a and -aa are incomptable with with --minind values greater than zero.
*This program expects SAMtools executable to be in the users's PATH.

```
