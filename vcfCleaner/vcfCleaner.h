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
	std::ofstream* write (std::ofstream& outstream);

	std::string contig;
	unsigned int start;
	unsigned int end;
	int entries; // number of times site appears in vcf
};

struct vcfrecord {
	vcfrecord(size_t flagcapacity = 0);
	void newEntry (std::string* vcfline, std::string* id, unsigned int* position,
			char altallele, double allelef, std::string* f = NULL);

	std::string entry;
	std::string flags;
	std::string contig;
	unsigned int pos;
	char alt;
	double af;
};

// functions

void maininfo ();

void gatkinfo (int &biallelic, int &allsites, int &allsites_vcf, unsigned int &maxcov, unsigned int &mincov, unsigned int &minind_cov,
	unsigned int &minind, unsigned int &mingeno, double &rms_mapq, double &mqRankSum, double &posbias, double &strandbias, double &baseqbias, double &qual,
	double &varqual_depth, double &hetexcess, int &varonly, double &maf, int& rmvIndels, int &verbose);

int gatkvcf (int argc, char** argv, std::ifstream &invcf, std::ofstream &outvcf, std::ofstream &passpos, std::ofstream &failpos);

int parseGATKargs (int argc, char** argv, std::ifstream &invcf, std::ofstream &outvcf, std::ofstream &passpos, std::ofstream &failpos,
	int &biallelic, int &allsites, int &allsites_vcf, unsigned int &maxcov, unsigned int &mincov, unsigned int &minind_cov,
	unsigned int &minind, unsigned int &mingeno, double &rms_mapq, double &mqRankSum, double &posbias, double &strandbias, double &baseqbias,
	double& qual, double &varqual_depth, double &hetexcess, int &varonly, double &maf, int &rmvIndels, int &verbose, int& infmt);

int parseFormat (std::vector<std::string> &vcfvec, int* index);

int extractIndInfo (std::vector<std::string> &vcfvec, size_t* indcounts, unsigned int min_indcov, int verbose);

//size_t nCoveredInd (std::vector<std::string> &vcfvec, unsigned int min_indcov, int* rv);

void siteType (std::vector<std::string> & vcfvec, int* isMultiAllelic, int* isIndel);

region* updateRegion (vcfrecord& site, region& goodpos, std::ofstream& outstream);

std::ofstream* writeBads (std::string& flags, std::string &contig, const unsigned int* pos, std::ofstream& outstream);

void checkGatkInfo(std::vector<std::string> &info, int n, std::string* flags, const unsigned int &mincov, const unsigned int &maxcov,
	const double &rms_mapq, const double &mqRankSum, const double &posbias, const double &strandbias, const double &baseqbias,
	const double & varqual_depth, const double &hetexcess);

double getMaf (const std::vector<std::string> &vcfvec, int verbose);

void printUserArgs (const char* invcf_name, std::string &outvcf_name, std::string &goodpos_name, std::string &badpos_name,
	int &biallelic, int &allsites, int &allsites_vcf, unsigned int &maxcov, unsigned int &mincov, unsigned int &minind_cov,
	unsigned int &minind, unsigned int &mingeno, double &rms_mapq, double &mqRankSum, double &posbias, double &strandbias, double &baseqbias,
	double &qual, double &varqual_depth, double &hetexcess, int &varonly, double &maf, int &rmvIndels);

#endif /* VCFCLEANER_H_ */
