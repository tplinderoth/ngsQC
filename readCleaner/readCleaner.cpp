/*
 * readCleaner.cpp
 *
 *  Created on: Jun 12, 2016
 *      Author: tyler
 */

#include "readCleaner.h"
#include "generalUtils.h"
#include "seqstat.h"
#include <sstream>

int main (int argc, char** argv) {
	//define variables
	int rc = 0;
	parseArg userin;

	// parse user input
	if ((rc=userin.getInput(argc, argv, version)) == 1)
		return 0;
	else if (rc < 0)
	{
		exitMessage();
		return 1;
	}

	// check dependencies
	std::cerr << "Checking for external dependencies...\n";
	if ((rc=depCheck(userin.eval())) != 0)
	{
		exitMessage();
		return rc;
	}

	// process fastq files
	if ((rc=cleanFastq(&userin)) != 0)
	{
		exitMessage();
		return rc;
	}

	return 0;
}

int depCheck (int doqual) {

	const char* exec = NULL;

	if (system("which super_deduper"))
		exec = "super_deduper";
	if (system("which cutadapt"))
		exec = "cutadapt";
	if (system("which pear"))
		exec = "pear";
	if (system("which bowtie2"))
		exec = "bowtie2";
	if (system("which bowtie2-build"))
		exec = "bowtie2-build";
	if (doqual)
	{
		if (system("which fastqc"))
			exec = "fastqc";
	}

	if (exec)
	{
		std::cerr << "Dependency '"<< exec << "' is not executable\n";
		return 1;
	}

	return 0;
}

void exitMessage ()
{
	std::cerr << "--> exiting\n";
}


int cleanFastq (parseArg* userarg) {
	bfs::path dir(userarg->fqdir());
	std::string scratchdir (userarg->outdir() + "/cleantmp");
	std::string evaldir (userarg->outdir() + "/evaluation");
	std::string lib;
	std::string fqfiles [3];
	bool ispair = userarg->pereads();
	int nfiles;

	if (bfs::exists(dir))
	{
		if (!bfs::is_directory(dir))
		{
			std::cerr << dir.string() << " is not a directory\n";
			return -1;
		}
	}
	else
	{
		std::cerr << "Directory with fastq files '" << dir.string() << "' doesn't exist\n";
		return -1;
	}

	// index contaminants file if necessary
	if (!userarg->contamIndex())
	{
		if (bowtie2::bowtie2Build(userarg->contamfile()))
			return -1;
	}

	// create working directory
	if (!bfs::create_directory(scratchdir))
	{
		std::cerr << "Couldn't create temporary directory " << scratchdir << "\n";
		return -1;
	}

	// create evaluation directory
	if (userarg->eval())
	{
		if (!bfs::create_directory(evaldir))
		{
			std::cerr << "Couldn't create directory for storing FastQC output: " << evaldir << "\n";
			return -1;
		}
	}

	// loop through fastq files
	for (bfs::directory_iterator diriter(dir); diriter != bfs::directory_iterator(); ++diriter)
	{
		if (ispair)
		{
			if (!isForwardRead(&diriter))
				continue;
		}

		if (!lib.empty())
			lib.clear();
		lib = libName(&diriter, ispair);
		nfiles = ispair ? 2 : 1;
		std::cerr << "Processing library " << lib << "\n";

		// remove duplicate reads
		if (setFastqNames1(fqfiles, diriter->path().string(), ispair))
			return -1;
		if (superdeduper::rmvDuplicates(fqfiles, nfiles, &scratchdir, &lib, ispair))
		{
			std::cerr << "Duplicate removal for library " << lib << " failed\n";
			return -1;
		}

		// trim for adapters and quality
		if (cutadapt::trimAdapters(fqfiles, nfiles, userarg->adapter(), &scratchdir, &lib, userarg->minlength(), userarg->minquality(), userarg->missingThresh(), ispair))
		{
			std::cerr << "Adapter and quality trimming for library " << lib << " failed\n";
			return -1;
		}
		if (superdeduper::deleteDupFiles(&scratchdir, &lib, ispair))
			return -1;

		// filter low complexity reads
		if (filterComplexity(fqfiles, userarg->dust(), &scratchdir, &lib, ispair))
		{
			std::cerr << "Removing low complexity reads from library " << lib << " failed\n";
			return -1;
		}
		if (cutadapt::deleteTrimFiles(&scratchdir, &lib, ispair))
			return -1;

		// merge overlapping reads
		if (ispair)
		{
			if (pear::mergeReads(fqfiles, &scratchdir, &lib, userarg->minlength(), userarg->missingThresh(), userarg->phred_offset(), userarg->nthreads()))
			{
				std::cerr << "Merging reads for library " << lib << " failed\n";
				return -1;
			}
			++nfiles;

			if (deleteComplexityFiles(&scratchdir, &lib, ispair))
				return -1;
		}

		// remove contaminant reads
		if (bowtie2::rmvContamination(fqfiles, nfiles, userarg->contamfile(), &scratchdir, &lib, userarg->phred_offset(), userarg->nthreads(), ispair))
		{
			std::cerr << "Removing contaminant reads for library " << lib << " failed\n";
			return -1;
		}

		if (ispair)
		{
			if (pear::deleteMergeFiles(&scratchdir, &lib))
				return -1;
		}
		else
		{
			if (deleteComplexityFiles(&scratchdir, &lib, ispair))
				return -1;
		}

		// name cleaned fastq files
		if (setFinalCleaned(fqfiles, userarg->outdir(), &lib, ispair))
		{
			std::cerr << "Generating final, cleaned fastq files for library " << lib << " failed\n";
			return -1;
		}

		if (bowtie2::deleteMappingFiles(&scratchdir, &lib, ispair))
			return -1;

		// perform evaluation of cleaned reads with fastqc
		if (userarg->eval())
		{
			if (fastqc::qualityEvaluation(fqfiles, &evaldir, ispair))
			{
				std::cerr << "Quality evaluation of library " << lib << " failed\n";
				return -1;
			}
		}
	}

	// remove working directory
	if (!bfs::remove_all(scratchdir))
		std::cerr << "Warning: Couldn't remove temporary directory " << scratchdir << "\n";

	std::cerr << "Finished cleaning reads!\n";
	return 0;
}

int setFinalCleaned (std::string fqfiles [3], const std::string& dir, const std::string* lib, const bool pe) {
	std::string prefix = dir + "/" + *lib;

	// check for good input and that files exist
	for (int i= 0; i<3; ++i)
	{
		if (!pe && i>0)
			break;
		if (!fqfiles[i].empty())
		{
			if (!bfs::exists(fqfiles[i]))
			{
				std::cerr << fqfiles[i] << " does not exist\n";
				return -1;
			}
		}
		else
		{
			std::cerr << "Missing fastq file name in call to setFinalCleaned\n";
			return -1;
		}
	}

	// Generate final, cleaned files
	if (newFileName(fqfiles[0], pe ? (prefix + "_clean_R1.fastq").c_str() : (prefix + "_clean.fastq").c_str()))
		return -1;

	if (pe)
	{
		if (newFileName(fqfiles[1], (prefix + "_clean_R2.fastq").c_str()))
			return -1;
		if (newFileName(fqfiles[2], (prefix + "_clean_u.fastq").c_str()))
			return -1;
	}

	return 0;
}

int newFileName (std::string& oldname, const char* newname) {
	if (rename(oldname.c_str(), newname))
	{
		std::cerr << "Unable to change " << oldname << " to " << newname << "\n";
		return -1;
	}
	oldname.clear();
	oldname = newname;
	return 0;
}

int superdeduper::rmvDuplicates (std::string infiles [3], const int nfiles, const std::string* outdir, const std::string* lib, const bool pe)
{
	static std::stringstream rmvdupCommand;
	static std::string outfile;
	int start = 5;
	int length = 10;

	// check that input files are not empty
	if (emptyFiles(infiles, nfiles))
		return -1;

	if (!rmvdupCommand.str().empty())
	{
		rmvdupCommand.clear();
		rmvdupCommand.str(std::string());
	}

	if (!outfile.empty())
		outfile.clear();
	outfile = *outdir + "/" + *lib;

	if (pe)
		rmvdupCommand << "super_deduper -1 " << infiles[0] << " -2 " << infiles[1] << " -p " << outfile << " -s " << start << " -l " << length;
	else
		rmvdupCommand << "super_deduper -U " << infiles[0] << " -p " << outfile << " -s " << start << " -l " << length;
	std::cerr << rmvdupCommand.str() << "\n";

	// run super_deduper
	if(system(rmvdupCommand.str().c_str()))
		return -1;

	// set outfile names
	dupFastqOut(infiles, outdir, lib, pe);

	return 0;
}

void superdeduper::dupFastqOut (std::string names [3], const std::string* outdir, const std::string* lib, const bool pe) {
	clearFileNames(names);
	names[0] = *outdir + "/" + *lib + "_nodup_PE1.fastq";
	if (pe)
		names[1] = *outdir + "/" + *lib + "_nodup_PE2.fastq";
}

int superdeduper::deleteDupFiles (const std::string* dir, const std::string* lib, const bool pe) {
	std::stringstream file;
	for(int i=0; i<2; ++i)
	{
		if (!pe && i>0)
			break;
		file << *dir << "/" << *lib << "_nodup_PE" << i+1 << ".fastq";
		if (remove(file.str().c_str()))
		{
			std::cerr << "Couldn't delete temporary file " << file.str() << "\n";
			return -1;
		}
		file.str(std::string());
	}
	return 0;
}

int cutadapt::deleteTrimFiles (const std::string* dir, const std::string* lib, const bool pe) {
	const char mode [] = {'a', 'g'};
	const int nmode = sizeof(mode)/sizeof(*mode);
	int i,j;
	std::stringstream file;
	for (i=0; i<2; ++i)
	{
		if (!pe && i>0)
			break;
		for (j=0; j<nmode; ++j)
		{
			file << *dir << "/" << *lib << "_trimmed" << mode[j] << "_R" << i+1 << ".fastq";
			if (remove(file.str().c_str()))
			{
				std::cerr << "Couldn't delete temporary file " << file.str() << "\n";
				return -1;
			}
			file.str(std::string());
		}
	}
	return 0;
}

int deleteComplexityFiles (const std::string* dir, const std::string* lib, const bool pe) {
	const char* file;

	file = (*dir + "/" + *lib + "_complex_R1.fastq").c_str();
	if (remove(file))
		return deleteFailMessage(file);

	if (pe)
	{
		file = (*dir + "/" + *lib + "_complex_R2.fastq").c_str();
		if (remove(file))
			return deleteFailMessage(file);
	}

	return 0;
}

int pear::deleteMergeFiles (const std::string* dir, const std::string* lib) {
	const char* file;
	std::string prefix(*dir + "/" + *lib);

	file = (prefix + ".unassembled.forward.fastq").c_str();
	if (remove(file))
		return deleteFailMessage(file);

	file = (prefix + ".unassembled.reverse.fastq").c_str();
	if (remove(file))
		return deleteFailMessage(file);

	file=(prefix + ".assembled.fastq").c_str();
	if (remove(file))
		return deleteFailMessage(file);

	file=(prefix + ".discarded.fastq").c_str();
	if (remove(file))
		return deleteFailMessage(file);

	return 0;
}

int bowtie2::deleteMappingFiles (const std::string * dir, const std::string* lib, const bool pe) {
	const char* file;
	std::string prefix(*dir + "/" + *lib);

	// delete concordantly mapped paired read fastq files
	if (pe)
	{
		file = (prefix + "_paired.1.fastq").c_str();
		if (remove(file))
			return deleteFailMessage(file);
		file = (prefix + "_paired.2.fastq").c_str();
		if (remove(file))
			return deleteFailMessage(file);
	}

	// delete SAM file
	file = (prefix + ".sam").c_str();
	if (remove(file))
		return deleteFailMessage(file);

	return 0;
}

int deleteFailMessage(const char* file) {
	std::cerr << "Could not delete temporary file " << file << "\n";
	return -1;
}

int cutadapt::trimAdapters (std::string infiles [3], const int nfiles, const std::string adapter_type, const std::string* outdir, const std::string* lib, const int minlength, const int minqual, const double pmissing, bool pe) {

	// adapter trimming parameters
	static Adapter adapter;
	const double errtol = 0.1; // Error tolerance, # allowable mismatches = length of match * errtol (-e)
	const int overlap = 1; // minimum number of bases to match between adapter and read to consider trimming (-O)
	const int ntimes = 2; // search for up to ntimes occurrences of the adapter in a read (-n)
	const char trimtype [] = {'a', 'g'}; // cutadapt read trimming mode
	const int nrounds = sizeof(trimtype)/sizeof(*trimtype);

	// read filtering parameters
	const std::string pairfilter ("any"); // which of the paired end reads have to match the filtering criterion to filter

	// trim
	static std::stringstream trimCommand;
	static std::string outnames [3];
	bool revcompliment;
	int rv;

	for (int i=0; i<nrounds; ++i)
	{
		// check for good input
		if (emptyFiles(infiles, nfiles))
			return -1;

		// set output names
		setTrimOutput(outnames, outdir, lib, trimtype[i], pe);

		// prepare command buffer for new command
		if (!trimCommand.str().empty())
		{
			trimCommand.clear();
			trimCommand.str(std::string());
		}

		revcompliment = trimtype[i] == 'a' ? false : true;
		if (pe)
		{
			trimCommand << "cutadapt -" << trimtype[i] << " " << adapter.seq(adapter_type, revcompliment) << " -" << static_cast<char>(toupper(trimtype[i]))
					<< " " << adapter.seq(adapter_type, revcompliment) << " -n " << ntimes << " -e " << errtol << " -O " << overlap
					<< " -m " << minlength << " -q " << minqual << "," << minqual << " --max-n " << pmissing << " --pair-filter " << pairfilter
					<< " -o " << outnames[0] << " -p " << outnames[1] << " " << infiles[0] << " " << infiles[1];
		}
		else
		{
			trimCommand << "cutadapt -" << trimtype[i] << " " << adapter.seq(adapter_type, revcompliment) << " -n " << ntimes << " -e " << errtol
					<< " -O " << overlap << " -m " << minlength << " -q " << minqual << "," << minqual << " --max-n " << pmissing
					<< " -o " << outnames[0] << " " << infiles[0];
		}
		std::cerr << trimCommand.str() << "\n";
		rv = (system(trimCommand.str().c_str()));
		if (rv)
			return -1;

		// set new input name
		trimFastqOut(infiles, outnames);
	}

	return rv;
}

int filterComplexity (std::string infiles [3], const double dustcutoff, const std::string* outdir, const std::string* lib, const bool pe) {
	std::fstream infile;
	std::fstream outfile;
	std::string outname(*outdir + "/" + *lib);
	static std::vector<std::string> badreads;
	unsigned int nreads = 0;

	// allocate space for storing low complexity reads
	if (!getFILE(infile, infiles[0].c_str(), "in"))
		return -1;
	nreads = std::count(std::istreambuf_iterator<char>(infile), std::istreambuf_iterator<char>(), '\n')/4;
	infile.clear();
	infile.seekg(0, std::ios::beg);
	if (badreads.size() > 0)
		badreads.clear();
	if (badreads.capacity() < nreads)
		badreads.reserve(nreads);

	// create temporary file for partially filtered read1 sequences
	outname += pe ? "_complex_R1_tmp.fastq" : "_complex_R1.fastq";
	std::string fin_tmp(outname); // for use with paired end reads
	if(!getFILE(outfile,outname.c_str(), "out"))
		return -1;

	// find low complexity reads in file1
	std::cerr << "Finding low complexity reads in " << infiles[0] << " ...\n";
	if (findLowComplexity(infile, outfile, &badreads, dustcutoff, true, &infiles[0]))
		return -1;
	infile.close();
	outfile. close();

	// done filtering complexity for single-end file
	if (!pe)
	{
		infiles[0].clear();
		infiles[0] = outname;
		return 0;
	}

	// allocate space for flagging reverse low complexity reads
	static std::vector<std::string> badreads2;
	if (badreads2.size() > 0)
		badreads2.clear();
	if (badreads2.capacity() < nreads)
		badreads2.reserve(nreads);

	// find low complexity reads in file
	if (!getFILE(infile, infiles[1].c_str(), "in"))
		return -1;

	outname.clear();
	outname = *outdir + "/" + *lib + "_complex_R2.fastq";
	if (!getFILE(outfile, outname.c_str(), "out"))
		return -1;

	std::cerr << "Finding low complexity reads in " << infiles[1] << " ...\n";
	std::sort(badreads.begin(), badreads.end()); // flagged read vector must be sorted prior to calling findLowComplexity
	if (findLowComplexity(infile, outfile, &badreads, dustcutoff, true, &infiles[1], &badreads2))
		return -1;
	infile.close();
	outfile.close();

	// record read2 output name
	infiles[1].clear();
	infiles[1] = outname;

	// remove forward mates of discarded low complexity reverse reads
	if (!getFILE(infile, fin_tmp.c_str(), "in"))
		return -1;

	outname.clear();
	outname = *outdir + "/" + *lib + "_complex_R1.fastq";
	if (!getFILE(outfile, outname.c_str(), "out"))
		return -1;

	std::cerr << "Ensuring that complex reads are properly paired in forward and reverse read files ...\n";
	// combine flagged forward and reverse reads and sort
	badreads.insert(badreads.end(), badreads2.begin(), badreads2.end());
	std::sort(badreads.begin(), badreads.end());

	if (findLowComplexity(infile, outfile, &badreads, dustcutoff, false, &fin_tmp))
		return -1;

	infile.close();
	outfile.close();
	// record read1 output name
	infiles[0].clear();
	infiles[0]=outname;

	// delete temporary R1 file
	if (remove(fin_tmp.c_str()))
		std::cerr << "Warning: Failure to delete temporary file " << fin_tmp << "\n";

	return 0;
}

int pear::mergeReads (std::string infiles [3], std::string* outdir, std::string* lib, int minlength, double missing, int phredbase, int nthreads) {
	double pval = 0.01; // p-value for statistical test (-p)
	int min_overlap = 6; // minimum overlap for assembling reads (-v)
	int maxlength = 0; // max length of assemblies, 0 disables (-m)
	int qual_cutoff = 0; // quality score threshold for trimming low quality part of a read (-q)
	int test = 1; // statistical test method (-g)
	int score = 2; // scoring method (-s)
	int qualcap = 41; // quality score upper bound (-c)

	// prepare command stream
	static std::stringstream mergeCommand;
	if (!mergeCommand.str().empty())
	{
		mergeCommand.clear();
		mergeCommand.str(std::string());
	}

	std::string outname = *outdir + "/" + *lib;

	mergeCommand << "pear -f " << infiles[0] << " -r " << infiles[1] << " -o " << outname << " -p " << pval << " -v " << min_overlap
			<< " -m " << maxlength << " -n " << minlength << " -q " << qual_cutoff << " -u " << missing << " -g " << test << " -s " << score
			<< " -b " << phredbase << " -c " << qualcap << " -j " << nthreads;

	std::cerr << mergeCommand.str() << "\n";

	// run PEAR
	if (system(mergeCommand.str().c_str()))
		return -1;

	// set outfile names
	mergeFastqOut(infiles, outdir, lib);

	return 0;
}

void pear::mergeFastqOut (std::string files [3], const std::string* outdir, const std::string* lib) {
	clearFileNames(files);
	files[0] += *outdir + "/" + *lib + ".unassembled.forward.fastq";
	files[1] += *outdir + "/" + *lib + ".unassembled.reverse.fastq";
	files[2] += *outdir + "/" + *lib + ".assembled.fastq";
}

int bowtie2::rmvContamination (std::string infiles [3], const int nfiles, const std::string& contamfile, const std::string* outdir, const std::string* lib, int phredbase, int nthreads, const bool pe) {
	std::stringstream contamCommand;
	static std::string phred (phredbase == 33 ? "--phred33" : "--phred64");
	const std::string preset("--fast");
	std::string out_prefix (*outdir + "/" + *lib);
	const std::string pe_suffix("_paired.fastq");
	const std::string se_suffix("_unpaired.fastq");
	std::string mapfile(*outdir + "/" + *lib + ".sam");

	// check for empty input fastq files
	int isempty = 0;
	if ((isempty = emptyFiles(infiles, nfiles)) != 0 && isempty < 2)
		return -1;

	// ready stream for new command
	if (!contamCommand.str().empty())
	{
		contamCommand.clear();
		contamCommand.str(std::string());
	}

	// map reads to contamination file using bowtie2
	if (pe)
	{
		contamCommand << "bowtie2 -q --quiet --no-head -p " <<  nthreads << " " << phred << " " << preset
				<< " --al-conc " << out_prefix + pe_suffix << " --un " << out_prefix + se_suffix << " -x " << contamfile
				<< " -1 " << infiles[0] << " -2 " << infiles[1] << " -U " << infiles[2] << " -S " << mapfile;
	}
	else
	{
		contamCommand << "bowtie2 -q --quiet --no-head -p " << nthreads << " " << phred << " " << preset << " --un " << out_prefix + se_suffix
				<< " -x " << contamfile << " -U " << infiles[0] << " -S " << mapfile;
	}
	std::cerr << contamCommand.str() << "\n";
	if (system(contamCommand.str().c_str()))
		return -1;

	// Extract noncontaminant reads from mapping results
	if (pe)
	{
		std::string concordant_mapped(out_prefix + "_paired.1.fastq");
		if (doPairedUnmapped(concordant_mapped, infiles, outdir, lib))
		{
			std::cerr << "Error encountered during doPairedUnmapped for processing paired read files\n";
			return -1;
		}
	}

	std::string ufile(out_prefix + "_unpaired.fastq");
	if (doSingleUnmapped(ufile, outdir, lib))
	{
		std::cerr << "Error encountered during doSingleUnmapped for processing unpaired read files\n";
		return -1;
	}

	// set output file names
	contamFastqOut(infiles, outdir, lib, pe);

	return 0;
}

void bowtie2::contamFastqOut (std::string files [3], const std::string* outdir, const std::string* lib, const bool pe) {
	clearFileNames(files);

	if (pe)
	{
		files[0] = *outdir + "/" + *lib + "_nocontam_R1.fastq";
		files[1] = *outdir + "/" + *lib + "_nocontam_R2.fastq";
		files[2] = *outdir + "/" + *lib + "_nocontam_u.fastq";
	}
	else
		files[0] = *outdir + "/" + *lib + "_nocontam_u.fastq";
}

int bowtie2::doSingleUnmapped (const std::string& mapfile, const std::string* outdir, const std::string* lib) {
	std::string newname(*outdir + "/" + *lib + "_nocontam_u.fastq");

	std::cerr << "Finding noncontaminant unpaired reads\n";

	if (mapfile.empty())
	{
		std::cerr << "No unmapped read file name in call to doSingleUnmapped\n";
		return -1;
	}

	if (!bfs::exists(mapfile))
	{
		std::cerr << "Unmapped, unpaired read file " << mapfile << " does not exist\n";
		return -1;
	}
	else
		return rename(mapfile.c_str(), newname.c_str());
}

int bowtie2::doPairedUnmapped (const std::string& mapfile, const std::string fqfiles [3], const std::string* outdir, const std::string* lib) {
	std::fstream fin;
	static std::vector<std::string> mapped;

	if (mapfile.empty())
	{
		std::cerr << "No name for concordantly mapped, paired reads file in call to getPairedUnmapped\n";
		return -1;
	}

	// prepare vector for storing flagged reads
	if (!getFILE(fin, mapfile.c_str(), "in"))
		return -1;
	unsigned int nreads = std::count(std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>(),'\n')/4;
	fin.clear();
	fin.seekg(0, std::ios::beg);
	if (mapped.size() > 0)
		mapped.clear();
	if (mapped.capacity() < nreads)
		mapped.reserve(nreads);

	std::cerr << "Finding noncontaminant paired reads\n";

	// insert mapped read IDs into vector of flagged reads
	if (nreads > 0)
	{
		if(flagMapped(fin, &mapped, mapfile))
			return -1;
		std::sort(mapped.begin(), mapped.end()); // sort vector for binary search
	}
	fin.close();

	// write unmapped PE reads to output file
	std::fstream fout;
	std::stringstream outname;
	for (int i=0; i<2; ++i)
	{
		if (!getFILE(fin, fqfiles[i].c_str(), "in"))
			return -1;
		outname << *outdir << "/" << *lib << "_nocontam_R" << i+1 << ".fastq";
		if(!getFILE(fout, outname.str().c_str(), "out"))
			return -1;
		outputUnmappedPE(fin, fout, &mapped, fqfiles[i]);
		fin.close();
		fout.close();
		outname.str(std::string());
	}

	return 0;
}

int bowtie2::outputUnmappedPE (std::fstream& infile, std::fstream& outfile, std::vector<std::string>* badreads, const std::string& in_name) {
	// badreads must be sorted

	// check for empty input file
	if (infile.peek() == std::ifstream::traits_type::eof())
	{
		std::cerr << in_name << " is an empty file in call to outputUnmappedPE\n";
		return -1;
	}

	// find noncontaminant reads and output them
	std::string fqline;
	std::string read_id;
	int i, n;
	bool iscontam;
	n=0;
	while (getline(infile,fqline))
	{
		if (fqline.empty())
		{
			std::cerr << in_name << " is an incomplete file\n";
			return -1;
		}
		if (n%4 == 0)
		{
			iscontam = false;
			// get unique part of read header
			if (readID(fqline, read_id).empty())
			{
				std::cerr << "Could not extract unique part of read header for file " << in_name << "\n";
				return -1;
			}

			// check if read is a contaminant
			iscontam = badreads->size() > 0 && std::binary_search(badreads->begin(), badreads->end(), read_id) ? true : false;

			// output read
			if (!iscontam)
				outfile << fqline << "\n";
			for (i=0; i<3; ++i)
			{
				if (getline(infile,fqline))
				{
					if (fqline.empty())
					{
						std::cerr << in_name << " is an incomplete file\n";
						return -1;
					}
					if (!iscontam)
						outfile << fqline << "\n";
				}
				else
				{
					std::cerr << in_name << " is an incomplete file\n";
					return -1;
				}
				++n;
			}
		}
		++n;
	}

	return 0;
}

int bowtie2::flagMapped (std::fstream& fs, std::vector<std::string>* reads, const std::string& filename) {
	std::string fqline;
	std::string id; // is static
	int n = 0;
	while(getline(fs,fqline))
	{
		if (fqline.empty())
		{
			std::cerr << filename << " is an incomplete file\n";
			return -1;
		}
		if (n%4 == 0)
		{
			if (readID(fqline, id).empty())
			{
				std::cerr << "Could not extract unique part of read header in call to flagMapped on " << filename  << "\n";
				return -1;
			}
			reads->push_back(id);
		}
		++n;
	}

	return 0;
}

int fastqc::qualityEvaluation (const std::string infiles [3], const std::string* outdir, const bool pe) {
	static std::stringstream evalCommand;
	std::fstream fin;

	for (int i=0; i<3; ++i)
	{
		if (!pe && i>0)
			break;

		// check if input file is empty
		if(!getFILE(fin, infiles[i].c_str(), "in"))
			return -1;
		if (fin.peek() == std::ifstream::traits_type::eof())
		{
			std::cerr << infiles[i] << " is an empty file... Skipping evaluation\n";
			fin.close();
			continue;
		}
		fin.close();

		// ready command stream for new input
		if (!evalCommand.str().empty())
		{
			evalCommand.clear();
			evalCommand.str(std::string());
		}

		// generate command
		evalCommand << "fastqc -o " << *outdir << " " << infiles[i];

		// run fastQC
		std::cerr << evalCommand.str() << "\n";
		if (system(evalCommand.str().c_str()))
		{
			std::cerr << "FastQC evaluation failed for file " << infiles[i] << "\n";
			return -1;
		}
	}

	return 0;
}

int findLowComplexity (std::fstream& infile, std::fstream& outfile, std::vector<std::string>* flagreads, const double cutoff,
		const bool calcdust, const std::string* in_name, std::vector<std::string>* r2flag) {
	int rv = 0;

	// check for empty input file
	if (infile.peek() == std::ifstream::traits_type::eof())
	{
		std::cerr << *in_name << " is an empty file in call to findLowComplexity\n";
		return -1;
	}

	// find and flag low complexity reads
	std::string fqline, seq;
	int n = -1;
	while(getline(infile, fqline))
	{
		++n;
		if (fqline.empty())
		{
			rv = -1;
			break;
		}
		if (n%4 == 0)
		{
			// check if read has already been flagged as low complexity
			if (flagreads->size() > 0 && isLowComplexityRead(fqline, flagreads))
				continue;
			// get read sequence
			if (getline(infile, seq))
			{
				++n;
				if (seq.empty())
				{
					rv = -1;
					break;
				}
				// calculate DUST score and flag if low complexity
				if (calcdust && seqstat::calcDustScore(seq) > cutoff)
				{
					if (addLowComplexityRead(fqline, r2flag ? r2flag : flagreads))
					{
						std::cerr << "Couldn't flag low complexity read in file " << in_name << "\n";
						return -1;
					}
					continue;
				}
				// write read information to output file
				outfile << fqline << "\n" << seq << "\n";
				for (int i=0; i<2; ++i)
				{
					if (getline(infile, fqline))
					{
						++n;
						if (!fqline.empty())
							outfile << fqline << "\n";
						else
						{
							rv = -1;
							break;
						}
					}
					else
					{
						std::cerr << "Missing lines in file " << *in_name << "\n";
						return -1;
					}
				}
				if (rv)
					break;
			}
			else
			{
				std::cerr << "Missing lines in file " << *in_name << "\n";
				return -1;
			}
		}
	}

	if (rv)
		std::cerr << "Unexpected empty line in file " << *in_name << "\n";

	return rv;
}

int addLowComplexityRead (const std::string& read_header, std::vector<std::string>* flagged_reads) {
	if (read_header.empty())
	{
		std::cerr << "Unable to flag low complexity read because its header is empty\n";
		return -1;
	}

	// get unique part of read name
	std::string id;
	if (readID(read_header, id).empty())
	{
		std::cerr << "Failed to extract unique part of read header in call to addLowComplexityRead\n";
		return -1;
	}

	// insert read id into vector of low complexity reads
	flagged_reads->push_back(id);

	return 0;
}

bool isLowComplexityRead (const std::string& read_header, const std::vector<std::string>* flagged_reads) {
	// input vector should be sorted

	if (read_header.empty())
	{
		std::cerr << "Warning: Empty header passed to isLowComplexityRead\n";
		return false;
	}

	// get unique part of read identifier
	std::string read_id;
	if (readID(read_header, read_id).empty())
	{
		std::cerr << "Warning: Unable to extract unique part of read header in call to isLowComplexityRead\n";
		return false;
	}

	// determine if read is present in the flagged reads vector
	return (std::binary_search(flagged_reads->begin(), flagged_reads->end(), read_id));
}

const std::string& readID (const std::string& read_header, std::string& read_id) {
	if (!read_id.empty())
		read_id.clear();

	if (read_header.empty())
	{
		std::cerr << "No read header in call to readID\n";
		return read_id;
	}

	// allocate enough space for unique part of read identifier
	read_id.reserve(read_header.length());

	// extract unique part of read identifier and remove the part pertaining to read 1/2
	std::string::const_iterator iter = read_header.begin();
	while (iter != read_header.end() && *iter != '/' && *iter != ' ' && *iter != '\t')
	{
		read_id.push_back(*iter);
		++iter;
	}

	return read_id;
}

void cutadapt::trimFastqOut (std::string a [3], std::string b [3]) {
	for (int i=0; i<3; ++i)
	{
		if (!a[i].empty())
			a[i].clear();
		if (!b[i].empty())
			a[i] = b[i];
	}
}

void cutadapt::setTrimOutput (std::string output [3], const std::string* dir, const std::string* lib, const char round, const bool pe) {
	clearFileNames(output);
	output[0] = *dir + "/" + *lib + "_trimmed" + round + "_R1.fastq";
	if (pe)
		output[1] = *dir + "/" + *lib + "_trimmed" + round + "_R2.fastq";
}

int setFastqNames1 (std::string names [3], std::string file, const bool pe) {
	clearFileNames(names);

	if (pe)
	{
		if (file.find("_R1.") != std::string::npos)
		{
			names[0] = file;
			names[1] = file.replace(file.rfind("_R1"), 3, "_R2");
		}
		else
		{
			std::cerr << "No 'R1' in the input file name " << file << ": check input file names\n";
			return -1;
		}
	}
	else
		names[0] = file;

	return 0;
}

bool isForwardRead (const bfs::directory_iterator* iter) {
	return (*iter)->path().filename().string().find("_R1.") == std::string::npos ? false : true;
}

int bowtie2::bowtie2Build (std::string file) {
	std::stringstream bowtie_command;
	bowtie_command << "bowtie2-build " << file << " " << file;
	if (system(bowtie_command.str().c_str()))
	{
		std::cerr << "Problem indexing " << file << " with bowtie2-build\n";
		return -1;
	}
	return 0;
}

std::string libName (const bfs::directory_iterator* fqiter, bool pe) {
	static std::string name;

	if (pe)
	{
		std::size_t endpos;
		if (!name.empty())
			name.clear();
		name = (*fqiter)->path().stem().string();
		if ((endpos = name.rfind("_R1")) != std::string::npos)
			return name.substr(0,endpos);
		else if ((endpos = name.rfind("_R2")) != std::string::npos)
			return name.substr(0,endpos);
		else
			return name;
	}

	return (*fqiter)->path().stem().string();
}

int emptyFiles (const std::string fnames [3], const int nfiles) {
	/*
	 * return values:
	 * -1 = missing file names or nonexistant file (error)
	 * 0 = no empty files
	 * 1 = empty single end or at least one paired end file
	 * 2 = empty merged read file
	 */

	static std::ifstream fs;

	for (int i=0; i<nfiles; ++i)
	{
		if (fnames[i].empty())
		{
			std::cerr << "Missing file name when checking for empty files\n";
			return -1;
		}
		if (bfs::exists(fnames[i]))
		{
			fs.open(fnames[i].c_str());
			if (is_empty(fs))
			{
				std::cerr << fnames[i] << " is empty\n";
				if (i < 2)
					return 1;
				else
					return 2;
			}
			fs.close();
		}
		else
		{
			std::cerr << fnames[i] << " does not exist\n";
			return -1;
		}
	}
	return 0;
}

void clearFileNames (std::string names [3]) {
	for (int i=0; i<3; ++i)
	{
		if (!names[i].empty())
			names[i].clear();
	}
}

Adapter::Adapter ()
	:	_truseq ("AGATCGGAAGAGC"),
	 	_truseqrc ("GCTCTTCCGATCT"),
	 	_nextera ("CTGTCTCTTATACACATCT"),
	 	_nexterarc("AGATGTGTATAAGAGACAG"),
	 	_smallrna("TGGAATTCTCGG"),
	 	_smallrnarc("CCGAGAATTCCA")
{
	truseq = std::make_pair(_truseq, _truseqrc);
	nextera = std::make_pair(_nextera, _nexterarc);
	smallrna = std::make_pair(_smallrna, _smallrnarc);
}

std::string Adapter::seq (std::string type, bool rc) {
	if (type == "truseq")
		return rc ? truseq.second : truseq.first;
	else if (type == "nextera")
		return rc ? nextera.second : nextera.first;
	else if (type == "smallrna")
		return rc ? smallrna.second : smallrna.first;
	else
		std::cerr << "Invalid adapter type in call to Adapter::seq()\n";
	return "";
}
