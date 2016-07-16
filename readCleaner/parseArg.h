/*
 * parseArg.h
 *
 *  Created on: Jun 12, 2016
 *      Author: tyler
 */

#ifndef PARSEARG_H_
#define PARSEARG_H_

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include <boost/filesystem.hpp>
#include <string>

namespace bfs = boost::filesystem;

class parseArg {
public:
	// FUNCTIONS //
	parseArg ();
	int getInput (const int c, char** v, const char* version);
	static bool bowtie2Index (const char* contamfile); // check if fasta is indexed by bowtie2
	int eval () const; // return _eval
	std::string fqdir () const; // return _fqdir
	std::string outdir ()const ; // return _outdir
	std::string contamfile () const; // return _contam
	std::string adapter () const; // return _adapter
	std::string trimjar () const; // return _trimjar
	double missingThresh () const; // return _missing
	double dust () const; // return _dust
	int minlength () const; // return _minlen
	int minquality () const; // return _minqual
	int pereads () const; // return _pe
	int nthreads () const; // return _nthreads
	int phred_offset() const; // return _phredbase
	bool contamIndex () const; // return _contamidx
	void setContamIndex (bool val); // set _contamidx
private:
	// FUNCTIONS //
	void helpinfo (const char* version);
	int lastChar (std::string str, const char c);
	// MEMBER VARIABLES //
	std::string _fqdir; // directory with fastq files
	std::string _outdir; // output directory
	std::string _contamf; // contaminants file
	std::string _adapter; // type of adapter
	std::string _trimjar; // trimmomatic executable
	int _minqual; // minimum quality for read trimming
	int _minlen; // minimum read length
	double _dust; // low complexity threshold
	double _missing; // missing data threshold
	int _pe; // is data paired end?
	int _eval; // perform evaluation with fastQC?
	int _phredbase; // base PHRED quality score
	int _nthreads; // number of threads to use
	bool _contamidx; // is contaminants file indexed?
};


#endif /* PARSEARG_H_ */
