/*
 * vcfCleaner.h
 */

#ifndef VCFCLEANER_H_
#define VCFCLEANER_H_

# include <string>
# include <vector>

// structures

struct region {
	region();
	void clear();
	void makeNew (std::string &seqid, unsigned int pos);
	std::fstream* write (std::fstream& outstream);

	std::string contig;
	unsigned int start;
	unsigned int end;
	int entries; // number of times site appears in vcf
};

struct vcfrecord {
	vcfrecord(size_t flagcapacity = 0);
	void newEntry (std::string* vcfline, std::string* id, unsigned int* position, char altallele, std::string* f = NULL);

	std::string entry;
	std::string flags;
	std::string contig;
	unsigned int pos;
	char alt;
};

// functions

void maininfo ();

void gatkinfo (int &biallelic, int &allsites, int &allsites_vcf, unsigned int &maxcov, unsigned int &mincov, unsigned int &minind_cov,
	unsigned int &minind, double &rms_mapq, double &mqRankSum, double &posbias, double &strandbias, double &baseqbias, double &qual,
	double &varqual_depth, double &hetexcess);

int gatkvcf (int argc, char** argv, std::fstream &invcf, std::fstream &outvcf, std::fstream &passpos, std::fstream &failpos);

int parseGATKargs (int argc, char** argv, std::fstream &invcf, std::fstream &outvcf, std::fstream &passpos, std::fstream &failpos,
	int &biallelic, int &allsites, int &allsites_vcf, unsigned int &maxcov, unsigned int &mincov, unsigned int &minind_cov,
	unsigned int &minind, double &rms_mapq, double &mqRankSum, double &posbias, double &strandbias, double &baseqbias,
	double &qual, double &varqual_depth, double &hetexcess);

size_t nCoveredInd (std::vector<std::string> &vcfvec, unsigned int min_indcov, int* rv);

int isMultiSNP (std::vector<std::string> &vcfvec);

region* updateRegion (vcfrecord& site, region& goodpos, std::fstream& outstream);

std::fstream* writeBads (std::string& flags, std::string &contig, const unsigned int* pos, std::fstream& outstream);

void checkGatkInfo(std::vector<std::string> &info, int n, std::string* flags, const unsigned int &mincov, const unsigned int &maxcov,
	const double &rms_mapq, const double &mqRankSum, const double &posbias, const double &strandbias, const double &baseqbias,
	const double & varqual_depth, const double &hetexcess);

void printUserArgs (const char* invcf_name, std::string &outvcf_name, std::string &goodpos_name, std::string &badpos_name,
	int &biallelic, int &allsites, int &allsites_vcf, unsigned int &maxcov, unsigned int &mincov, unsigned int &minind_cov,
	unsigned int &minind, double &rms_mapq, double &mqRankSum, double &posbias, double &strandbias, double &baseqbias,
	double &qual, double &varqual_depth, double &hetexcess);

#endif /* VCFCLEANER_H_ */
