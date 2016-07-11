/*
 * readCleaner.h
 *
 *  Created on: Jun 12, 2016
 *      Author: tyler
 *
 * Tested with the following dependency versions:
 * cutadapt v1.10
 * super_deduper v1.4
 * pear v0.9.10
 * bowtie2 v2.2.9
 * fastqc v0.11.4
 */

#ifndef READCLEANER_H_
#define READCLEANER_H_

#include "parseArg.h"
#include <utility>

const char* version = "0.0.1";

typedef std::pair<std::string, std::string> adaptseq;

struct Adapter {
public:
	Adapter();
	std::string seq (std::string type, bool rc); // return particular adapter sequence or its reverse compliment
	adaptseq truseq;
	adaptseq nextera;
	adaptseq smallrna;
private:
	// Adapters
	const char* _truseq; // truseq adapter sequence adjacent to insert
	const char* _truseqrc; // reverse compliment of truseq adapter
	const char* _nextera; // nextera adapter sequence adjacent to insert
	const char* _nexterarc; // reverse compliment of nextera adapter sequence
	const char* _smallrna; // illumina small RNA adapter sequence adjacent to insert
	const char* _smallrnarc; // reverse compliment of illumina small RNA adapter sequence
};

void exitMessage ();
int depCheck (int doqual); // check if dependencies are executable
int cleanFastq (parseArg* userarg); // wrapper to clean fastq files
std::string libName (const bfs::directory_iterator* fqiter, bool pe); // extracts library name
bool isForwardRead (const bfs::directory_iterator* iter); // checks if read is the forward read (using R1/R2 naming)
int emptyFiles (const std::string fnames [3], const int nfiles); // check if paired files or unpaired file is empty
int setFastqNames1 (std::string names [3], std::string file, const bool pe); // set input forward/reverse or unpaired fastq file names
const std::string& readID (const std::string& read_header, std::string& read_id); // extracts unique part of read header
int filterComplexity (std::string infiles [3], const double dustcutoff, const std::string* outdir, const std::string* lib, const bool pe); // filter out low complexity reads
int findLowComplexity (std::fstream& infile, std::fstream& outfile, std::vector<std::string>* flagreads, const double cutoff,
		const bool calcdust, const std::string* in_name, std::vector<std::string>* r2flag = NULL); // flag reads with high DUST score in a fastq file
bool isLowComplexityRead (const std::string& read_id, const std::vector<std::string>* flagged_reads); // checks if read is stored in vector of reads flagged as low complexity
int addLowComplexityRead (const std::string& read_header, std::vector<std::string>* flagged_reads); // insert read identifier into vector of low complexity reads
void clearFileNames (std::string names [3]); // clears array of file names
int deleteComplexityFiles (const std::string* dir, const std::string* lib, const bool pe); // delete complexity filter output files
int deleteFailMessage(const char* file); // failure to delete file error message
int setFinalCleaned (std::string fqfiles [3], const std::string& dir, const std::string* lib, const bool pe); // generate totally cleaned files
int newFileName (std::string& oldname, const char* newname); // change name of file and set the oldname string to the new name

namespace superdeduper {
int rmvDuplicates (std::string infiles [3], const int nfiles, const std::string* outdir, const std::string* lib, const bool pe); // remove duplicate reads with super_deduper
void dupFastqOut (std::string names [3], const std::string* outdir, const std::string* lib, const bool pe); // set super_deduper output file names
int deleteDupFiles (const std::string* outdir, const std::string* lib, const bool pe); // delete super_deduper output files
}

namespace cutadapt {
int trimAdapters (std::string infiles [3], const int nfiles, const std::string adapter_type, const std::string* outdir,
		const std::string* lib, const int minlength, const int minqual, const double pmissing, bool pe); // trim for adapters and quality
void setTrimOutput (std::string output [3], const std::string* dir, const std::string* lib, const char round, const bool pe); // set trimmed output file names
void trimFastqOut (std::string a [3], std::string b [3]); // copies names from array b into array a
int deleteTrimFiles (const std::string* dir, const std::string* lib, const bool pe); // delete temporary cutadapt files
}

namespace pear {
int mergeReads (std::string infiles [3], std::string* outdir, std::string* lib, int minlength, double missing, int phredbase, int nthreads); // merge overlapping paired end reads
void mergeFastqOut (std::string files [3], const std::string* outdir, const std::string* lib); // set PEAR output fastq file names
int deleteMergeFiles (const std::string* dir, const std::string* lib); // delete temporary pear files
}

namespace bowtie2 {
int bowtie2Build (std::string file); // build bowtie2 index
int rmvContamination (std::string infiles [3], const int nfiles, const std::string& contamfile, const std::string* outdir,
		const std::string* lib, int phredbase, int nthreads, const bool pe); // remove contaminant reads
int doPairedUnmapped (const std::string& mapfile, const std::string fqfiles [3], const std::string* outdir, const std::string* lib); // Process results from mapping PE reads to contaminant sequences
int outputUnmappedPE (std::fstream& infile, std::fstream& outfile, std::vector<std::string>* badreads, const std::string& in_name); // Find noncontaminant PE reads and write them to output
int flagMapped (std::fstream& fs, std::vector<std::string>* reads, const std::string& filename); // flag read pairs that mapped concordantly to contaminant file
int doSingleUnmapped (const std::string& mapfile, const std::string* outdir, const std::string* lib); // Process results from mapping unpaired reads to contaminant sequences
void contamFastqOut (std::string files [3], const std::string* outdir, const std::string* lib, const bool pe); // set names of contamination-free fastq files
int deleteMappingFiles (const std::string * dir, const std::string* lib, const bool pe); // delete temporary bowtie2 files
}

namespace fastqc {
int qualityEvaluation (const std::string infiles [3], const std::string* outdir, const bool pe); // call fastQC to evaluate read quality
}

#endif /* READCLEANER_H_ */
