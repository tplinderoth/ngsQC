/*
 * parseArg.cpp
 *
 *  Created on: Jun 12, 2016
 *      Author: tyler
 */

#include "parseArg.h"
#include <iomanip>
#include <sstream>

parseArg::parseArg()
	: _adapter("truseq"),
	  _minqual(20),
	  _minlen(36),
	  _dust(4.0),
	  _missing(0.5),
	  _pe(1),
	  _eval(1),
	  _phredbase(33),
	  _nthreads(1),
	  _contamidx(true)
{}

int parseArg::getInput (const int c, char** v, const char* version) {
	int argPos = 1;
	int inc = 0;
	std::string test;
	const std::string adapters [] = {"truseq", "nextera", "smallrna"};
	int nadapters = sizeof(adapters)/sizeof(*adapters);

	if (c < 2)
	{
		helpinfo(version);
		return 1;
	}

	while (argPos < c)
	{
		if (strcmp(v[argPos],"-fastqdir") == 0)
		{
			_fqdir = v[argPos+1];
			if (lastChar(_fqdir, '/') == 1)
				_fqdir = _fqdir.substr(0,_fqdir.length()-1);
			if (bfs::exists(_fqdir))
			{
				if (!bfs::is_directory(_fqdir))
				{
					std::cerr << _fqdir << " is not a directory\n";
					return -1;
				}
			}
			else
			{
				std::cerr << _fqdir << " does not exists\n";
				return -1;
			}
		}
		else if (strcmp(v[argPos],"-outdir") == 0)
		{
			_outdir = v[argPos+1];
			if (lastChar(_outdir, '/') == 1)
				_outdir = _outdir.substr(0,_outdir.length()-1);
			if (bfs::exists(_outdir))
			{
				if (!bfs::is_directory(_outdir))
				{
					std::cerr << _outdir << " is not a directory\n";
					return -1;
				}
			}
			else
			{
				if (!bfs::create_directory(_outdir))
				{
					std::cerr << "Couldn't create output directory " << _outdir << "\n";
					return -1;
				}
			}
		}
		else if (strcmp(v[argPos],"-contam") == 0)
		{
			_contamf = v[argPos+1];
			if (bfs::exists(_contamf))
			{
				if (!bfs::is_regular_file(_contamf))
				{
					std::cerr << _contamf << " is not a regular file\n";
					return -1;
				}
			}
			else
			{
				std::cerr << "Contaminants file " << _contamf << " does not exist\n";
				return -1;
			}
			setContamIndex(bowtie2Index(v[argPos+1]));
		}
		else if (strcmp(v[argPos],"-adapter") == 0)
		{
			_adapter = v[argPos+1];
			bool isadapter = false;
			for (int i=0; i<nadapters; ++i)
			{
				if (_adapter == adapters[i])
					isadapter = true;
			}
			if (!isadapter)
			{
				std::cerr << _adapter << " is an invalid type of adapter\n";
				return -1;
			}
		}
		else if (strcmp(v[argPos],"-minqual") == 0)
		{
			_minqual = atoi(v[argPos+1]);
			if (_minqual < 0)
			{
				std::cerr << "Minimum quality must be >= 0, check -minqual\n";
				return -1;
			}
		}
		else if (strcmp(v[argPos],"-minlength") == 0)
		{
			_minlen = atof(v[argPos+1]);
			if (_minlen < 1)
			{
				std::cerr << "Minimum read length must be >= 1, check -minlength\n";
				return -1;
			}
		}
		else if (strcmp(v[argPos],"-dust") == 0)
		{
			_dust = atof(v[argPos+1]);
			if (_dust < 0)
			{
				std::cerr << "DUST score must be >= 0, check -dust\n";
				return -1;
			}
		}
		else if (strcmp(v[argPos],"-missing") == 0)
		{
			_missing = atof(v[argPos+1]);
			if (_missing < 0.0 || _missing > 1.0)
			{
				std::cerr << "Proportion of missing data out of range, check -missing\n";
				return -1;
			}
		}
		else if (strcmp(v[argPos],"-quality_base") == 0)
		{
			_phredbase = atoi(v[argPos+1]);
			switch(_phredbase)
			{
				case 33 :
					break;
				case 64 :
					break;
				default :
				{
					std::cerr << "Only Phred quality scaling of 33 or 64 is valid, check -quality_base\n";
					return -1;
				}
			}
		}
		else if (strcmp(v[argPos],"-nthreads") == 0)
		{
			_nthreads = atoi(v[argPos+1]);
			if (_nthreads < 0)
			{
				std::cerr << "Number of threads cannot be negative, check -nthreads\n";
				return -1;
			}
		}
		else if (strcmp(v[argPos],"-pe"))
		{
			_pe = atoi(v[argPos+1]);
			switch(_pe)
			{
				case 0 :
					break;
				case 1 :
					break;
				default :
					std::cerr << "Is data paired end (-pe 1) or single end (-pe 0)?\n";
					return -1;
			}
		}
		else if (strcmp(v[argPos],"-eval"))
		{
			_eval = atoi(v[argPos+1]);
			switch (_eval)
			{
				case 0 :
					break;
				case 1 :
					break;
				default :
					std::cerr << "Run FastQC over cleaned data (-eval 1) or not (-eval 0)?\n";
					return -1;
			}
		}
		argPos += 2 + inc;
		inc = 0;
	}

	return 0;
}

bool parseArg::bowtie2Index (const char* contamfile) {
    std::string index(contamfile);
    std::string revindex(contamfile);
    int baselen = index.length();
    index += ".x.bt2";
    revindex += ".rev.x.bt2";
    int i = 1;

    for(i=1; i<=4; ++i)
    {
    		index[baselen+1] = static_cast<std::ostringstream*>( &(std::ostringstream() << i) )->str()[0];
    		if (!bfs::exists(index))
    			return false;
    }

    for (i=1; i<=2; ++i)
    {
    	revindex[baselen+5] = static_cast<std::ostringstream*>( &(std::ostringstream() << i) )->str()[0];
    	if (!bfs::exists(revindex))
    		return false;
    }

    return true;
}

void parseArg::helpinfo (const char* version) {
	int w = 14;
	std::cerr << "readCleaner version " << version << "\n";
	std::cerr << "\nUsage: readCleaner [arguments]\n"
	<< "\nDEPENDENCIES:\n"
	<< "super_deduper, cutadapt, pear, bowtie2, fastqc (optional)\n"
	<< "Note: Program executable names should be exactly as listed above\n"
	<< "\nINPUT:"
	<< "\n" << std::setw(w) << std::left << "-fastqdir" << std::setw(w) << "<string>" << "Directory containing fastq files [" << _fqdir << "]"
	<< "\n" << std::setw(w) << std::left << "-outdir" << std::setw(w) << "<string>" << "Directory to place cleaned fastq files in [" << _outdir << "]"
	<< "\n" << std::setw(w) << std::left << "-contam" << std::setw(w) << "<string>" << "Fasta format file of potential contaminant sequences [" << _contamf << "]"
	<< "\n" << std::setw(w) << std::left << "-adapter" << std::setw(w) << "<string>" << "Adapter type: truseq, nextera, smallrna [" << _adapter << "]"
	<< "\n" << std::setw(w) << std::left << "-minqual" << std::setw(w) << "<int>" << "Minimum base quality for trimming ends of reads using BWA algorithm [" << _minqual << "]"
	<< "\n" << std::setw(w) << std::left << "-minlength" << std::setw(w) << "<int>" << "Minimum read length [" << _minlen << "]"
	<< "\n" << std::setw(w) << std::left << "-dust" << std::setw(w) << "<double>" << "Maximum read DUST score for removing low complexity reads [" << _dust << "]"
	<< "\n" << std::setw(w) << std::left << "-missing" << std::setw(w) << "<double>" << "Maximum proportion of 'N' bases in reads [" << _missing << "]"
	<< "\n" << std::setw(w) << std::left << "-quality_base" << std::setw(w) << "<int>" << "Phred quality score offset (int = 33 or 64) [" << _phredbase << "]"
	<< "\n" << std::setw(w) << std::left << "-nthreads" << std::setw(w) << "<int>" << "Number of threads for multithreaded programs to use [" << _nthreads << "]"
	<< "\n" << std::setw(w) << std::left << "-pe" << std::setw(w) << "<int>" << "Data is paired end (int = 1) or unpaired (int = 0) [" << _pe << "]"
	<< "\n" << std::setw(w) << std::left << "-eval" << std::setw(w) << "<int>" << "Perform fastQC evaluation of cleaned reads (int = 1), or do not (int = 0) [" << _eval << "]"
	<< "\n\n";
}

int parseArg::lastChar (std::string str, const char c) {
	return (!str.empty() && str.length() > 0 && str[str.length()-1] == c) ? 1 : 0;
}

int parseArg::eval () const {return _eval;}

std::string parseArg::fqdir () const {return _fqdir;}

std::string parseArg::outdir () const {return _outdir;}

std::string parseArg::contamfile () const {return _contamf;}

std::string parseArg::adapter () const {return _adapter;}

double parseArg::missingThresh () const {return _missing;}

double parseArg::dust () const {return _dust;}

int parseArg::minlength () const {return _minlen;}

int parseArg::minquality () const {return _minqual;}

int parseArg::pereads () const {return _pe;}

bool parseArg::contamIndex () const {return _contamidx;}

void parseArg::setContamIndex (bool val) {_contamidx = val;}

int parseArg::nthreads () const {return _nthreads;}

int parseArg::phred_offset() const {return _phredbase;}
