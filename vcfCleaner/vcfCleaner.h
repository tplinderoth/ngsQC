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
	vcfrecord();
	void reserveFields (unsigned int nfields);
	void reserveFlags (unsigned int nflags);
	void clearFlags ();
	void setSite ();
	void setFlag (int i);
	std::string contig;
	unsigned int pos;
	std::string ref;
	std::string alt;
	double af;
	int qcfail; // number of set flags
	std::vector<std::string> fields; // tokenized vcf line
	std::vector< std::pair<std::string,int> > flags; // QC-fail flags
};

// functions

void maininfo ();

void gatkinfo (int &biallelic, int &allsites, int &allsites_vcf, unsigned int &maxcov, unsigned int &mincov, unsigned int &minind_cov,
	unsigned int &minind, unsigned int &mingeno, double &rms_mapq, double &mqRankSum, double &posbias, double &strandbias, double &baseqbias, double &qual,
	double &varqual_depth, double &hetexcess, int &varonly, double &maf, int &rmvIndels, int &siteOnly, int &printFilter, int &verbose);

int gatkvcf (int argc, char** argv, std::ifstream &invcf, std::ofstream &outvcf, std::ofstream &passpos, std::ofstream &failpos);

int parseGATKargs (int argc, char** argv, std::ifstream &invcf, std::ofstream &outvcf, std::ofstream &passpos, std::ofstream &failpos,
	int &biallelic, int &allsites, int &allsites_vcf, unsigned int &maxcov, unsigned int &mincov, unsigned int &minind_cov,
	unsigned int &minind, unsigned int &mingeno, double &rms_mapq, double &mqRankSum, double &posbias, double &strandbias, double &baseqbias,
	double& qual, double &varqual_depth, double &hetexcess, int &varonly, double &maf, int &rmvIndels, int &siteOnly, int &printFilter, int &verbose, int& infmt);

int parseFormat (std::vector<std::string> &vcfvec, int* index);

int extractIndInfo (std::vector<std::string> &vcfvec, size_t* indcounts, unsigned int min_indcov, int verbose);

//size_t nCoveredInd (std::vector<std::string> &vcfvec, unsigned int min_indcov, int* rv);

void siteType (std::vector<std::string> & vcfvec, int* isMultiAllelic, int* isIndel);

region* updateRegion (vcfrecord* site, region& goodpos, std::ofstream& outstream);

std::ofstream* writeBads (vcfrecord* vcfrec, std::ofstream& outstream);

void checkGatkInfo(vcfrecord* vcfrec, const unsigned int &mincov, const unsigned int &maxcov,
	const double &rms_mapq, const double &mqRankSum, const double &posbias, const double &strandbias, const double &baseqbias,
	const double & varqual_depth, const double &hetexcess);

double getMaf (const std::vector<std::string> &vcfvec, int verbose);

int recordSite (vcfrecord* site, region &keep, std::ofstream &vcfstream, std::ofstream &goodstream, std::ofstream &failstream,
		const int &allsites_vcf, const int &printFilter, const int &siteOnly, const int &varonly, const double &mafcutoff, const int &allsites);

void writeVcf (vcfrecord* vcfrec, std::ofstream &os, const int &siteOnly, const int &printFilter);

void printUserArgs (const char* invcf_name, std::string &outvcf_name, std::string &goodpos_name, std::string &badpos_name,
	int biallelic, int allsites, int allsites_vcf, unsigned int maxcov, unsigned int mincov, unsigned int minind_cov,
	unsigned int minind, unsigned int mingeno, double rms_mapq, double mqRankSum, double posbias, double strandbias, double baseqbias,
	double qual, double varqual_depth, double hetexcess, int varonly, double maf, int rmvIndels, int siteOnly, int printFilter);

#endif /* VCFCLEANER_H_ */
