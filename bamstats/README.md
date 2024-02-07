bamstats
========

## Description:

Extracts the following stats from a BAM file:

* Total depth
* Average map quality
* RMS map quality
* Fraction of reads with map quality of zero
* Average base quality
* RMS base quality
* Fraction of reads with base quality of zero

## Compile:

Once in directory containing bamstats.cpp

`g++ -O3 -std=c++11 -o bamstats bamstats.cpp`

## Usage:

`bamstats [options to samtools mpileup] <bam file>`

The only additional mpileup fields that this script works with is map quality
produced using samtools mpileup -s/--output-MQ.
