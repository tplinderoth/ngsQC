/*
 * vcfCleaner.cpp
 *
 * TODO:
 * 1) dynamically set infosize (number of INFO IDs) based on VCF header
 * 2) implement accurate multiallelic filtering, e.g. maf filtering
 * 3) Boost does not play well with bgzipped files (segfaults or leaves off the last line), so need to
 * switch over to using htslib to parse VCF. Gzipped files are fine.
 */

#include "vcfCleaner.h"
#include <cstring>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

int main (int argc, char** argv) {
	int rv = 0;

	std::ifstream invcf; // input vcf
	std::ofstream outvcf; // filtered output vcf
	std::ofstream pass; // quality passed list of sites
	std::ofstream fail; // failed sites

	if (argc < 2) {
		maininfo();
	} else if (strcmp(argv[1], "gatk") == 0) {
		rv = gatkvcf (argc, argv, invcf, outvcf, pass, fail);
	} else {
		std::cerr << "\n" << argv[1] << " is an invalid VCF type\n";
		rv = 1;
		maininfo();
	}

	// close file streams
	if (invcf.is_open()) invcf.close();
	if (outvcf.is_open()) outvcf.close();
	if (pass.is_open()) pass.close();
	if (fail.is_open()) fail.close();

	return rv;
}

void maininfo () {
	int w=12;
	std::cerr << "\nvcfCleaner [VCF type] [arguments]\n"
	<< "\nVCF types:\n"
	<< std::setw(w) << std::left << "gatk" << "vcf produced by GATK HaplotypeCaller/GenotypeGVCFs workflow\n"
	<< "\n";
}

void gatkinfo (int &biallelic, int &allsites, int &allsites_vcf, unsigned int &maxcov, unsigned int &mincov, unsigned int &minind_cov,
	unsigned int &minind, unsigned int &mingeno, double &rms_mapq, double &mqRankSum, double &posbias, double &strandbias, double &baseqbias, double &qual,
	double &varqual_depth, double &hetexcess, int &varonly, double &maf, int &verbose) {

	int w1=20;
	int w2=8;

	std::cerr << "\nThis program is for filtering VCFs produced by GATK's HaplotypeCaller/GenotypeGVCFs workflow.\n"
	<< "Takes uncompressed and gzipped/bgzipped VCFs, however bgzipped reading can be buggy so use at your own risk.\n"
	<< "Indels are ignored and not included in files containing filtered sites.\n"
	<< "\nvcfCleaner gatk [arguments]\n\n"
	<< std::setw(w1) << std::left << "-vcf" << std::setw(w2) << std::left << "FILE" << "Input VCF file to filter ('-' for STDIN)\n"
	<< std::setw(w1) << std::left << "-out" << std::setw(w2) << std::left << "STRING" << "Output file name (prefix)\n"
	<< std::setw(w1) << std::left << "-biallelic" << std::setw(w2) << std::left << "0|1" << "0: keep multiallelic sites, 1: keep only biallelic sites, A [" << biallelic << "]\n"
	<< std::setw(w1) << std::left << "-allsites" << std::setw(w2) << std::left << "0|1" << "0: process only SNPs, 1: process all sites (except indels) [" << allsites << "]\n"
	<< std::setw(w1) << std::left << "-allsites_vcf" << std::setw(w2) << std::left << "0|1" << "Output VCF contains QC-passed 0: SNPs only, 1: all sites [" << allsites_vcf << "]\n"
	<< std::setw(w1) << std::left << "-maxcov" << std::setw(w2) << std::left << "INT" << "Max total site depth, D [" << maxcov << "]\n"
	<< std::setw(w1) << std::left << "-mincov" << std::setw(w2) << std::left << "INT" << "Min total site depth, C [" << mincov << "]\n"
	<< std::setw(w1) << std::left << "-minind_cov" << std::setw(w2) << std::left << "INT" << "Min individual coverage [" << minind_cov << "]\n"
	<< std::setw(w1) << std::left << "-minind" << std::setw(w2) << std::left << "INT" << "Min number of individuals covered by at least -minind_cov reads, U [" << minind << "]\n"
	<< std::setw(w1) << std::left << "-mingeno" << std::setw(w2) << std::left << "INT" << "Min number of individuals with called genotype, G [" << mingeno << "]\n"
	<< std::setw(w1) << std::left << "-rms_mapq" << std::setw(w2) << std::left << "FLOAT" << "Min RMS mapping quality, R [" << rms_mapq << "]\n"
	<< std::setw(w1) << std::left << "-mapqRankSum" << std::setw(w2) << std::left << "FLOAT" << "Max absolute Wilcoxon rank sum test Z-score of alt vs. ref read map quality, M [" << mqRankSum << "]\n"
	<< std::setw(w1) << std::left << "-posbias" << std::setw(w2) << std::left << "FLOAT" << "Max absolute Wilcoxon rank sum test Z-score of alt vs. ref read position bias, P [" << posbias << "]\n"
	<< std::setw(w1) << std::left << "-strandbias" << std::setw(w2) << std::left << "FLOAT" << "Max Phred-scaled Fisher's exact test p-value of strand bias, S [" << strandbias << "]\n"
	<< std::setw(w1) << std::left << "-baseqbias" << std::setw(w2) << std::left << "FLOAT" << "Max absolute Wilcoxon rank sum test Z-score of alt vs ref base qualities, B [" << baseqbias << "]\n"
	<< std::setw(w1) << std::left << "-qual" << std::setw(w2) << std::left << "FLOAT" << "Min Phred-scaled quality score of ALT assertion, Q [" << qual << "]\n"
	<< std::setw(w1) << std::left << "-varq_depth" << std::setw(w2) << std::left << "FLOAT" << "Min variant Phred-scaled confidence/quality by depth, V [" << varqual_depth << "]\n"
	<< std::setw(w1) << std::left << "-hetexcess" << std::setw(w2) << std::left << "FLOAT" << "Max Phred-scaled p-value for exact test of excess heterozygosity, H [" << hetexcess << "]\n"
	<< std::setw(w1) << std::left << "-varonly" << std::setw(w2) << std::left << "0|1" << "The INFO AF for SNP-only output VCF must be in range 0=[0,1], or 1=(0,1) [" << varonly << "]\n"
	<< std::setw(w1) << std::left << "-maf" << std::setw(w2) << std::left << "FLOAT" << "Minor allele frequency lower bound for SNP-only VCF [" << maf << "]\n"
	<< std::setw(w1) << std::left << "-verbose" << std::setw(w2) << std::left << "0|1|2" << "Level of warnings to issue: 0 = suppress all, 1 = site-level, 2 = individual-level [" << verbose << "]\n"
	<< "\nOther site QC fail flags:\n"
	<< "F, Unkown reference allele\n"
	<< "N, No data for any individuals\n"
	<< "\n";
}

int gatkvcf (int argc, char** argv, std::ifstream &invcf, std::ofstream &outvcf, std::ofstream &passpos, std::ofstream &failpos) {
	int rv = 0;

	// filtering parameters
	int biallelic = 1; // 0: keep multiallelic sites, 1: keep only biallelic sites
	int allsites = 0; // 0: process only SNPs, 1: process all sites (except indels)
	int allsites_vcf = 0; // output VCF contains QC-passed 0: SNPs only or 1: all sites
	unsigned int maxcov = 1000000; // maximum total site coverage (DP)
	unsigned int mincov = 2; // minimum total site coverage (DP)
	unsigned int minind_cov = 1; // minimum individual coverage (DP)
	unsigned int minind = 1; // minimum number of individuals with data
	unsigned int mingeno = 1; // minimum number of individuals with called genotype
	double rms_mapq = 40.0; // min RMS mapping quality (MQ)
	double mqRankSum = 12.5; // max absolute Wilcoxon rank sum test z-score of alt vs ref read map quality (MQRankSum)
	double posbias = 8.0; // max absolute Wilcoxon rank sum test z-score of alt vs ref read position bias  (ReadPosRankSum)
	double strandbias = 60.0; // max Phred-scale p-value using Fisher's exact test of strand bias (FS)
	double baseqbias = 21.0; // max absolute Wilcoxon rank sum test z-score of alt vs. ref base qualities (BaseQRankSum)
	double varqual_depth = 2.0; // min variant Phred-scaled confidence/quality by depth (QD)
	double qual = 30; // min Phred-scaled quality score of ALT assertion (QUAL field)
	double hetexcess = 40.0; // max Phred-scaled p-value for exact test of excess heterozygosity (ExcessHet)
	int varonly = 1; // The INFO for SNP-only VCF must have AF in 1=range (0,1), or 0=range [0,1]
	double mafcutoff = 0.0; // minor allele frequency lower bound for SNP-only VCF
	int verbose = 2; // amount of warnings outputted, 0=none, 1=site level, 2=individual level
	int infmt = 0; // 0=uncompressed vcf, 1=gzipped input, 2=uncompressed standard input

	if ((rv=parseGATKargs(argc, argv, invcf, outvcf, passpos, failpos, biallelic, allsites, allsites_vcf, maxcov, mincov, minind_cov,
			minind, mingeno, rms_mapq, mqRankSum, posbias, strandbias, baseqbias, qual, varqual_depth, hetexcess, varonly, mafcutoff, verbose, infmt))) {
		if (rv > 0)
			return 0;
		else if (rv < 0)
			return -1;
	}

	// set up zipped file reading if necessary
	std::streambuf *inbuf = NULL;
	boost::iostreams::filtering_streambuf<boost::iostreams::input> zipbuf;
	zipbuf.push(boost::iostreams::gzip_decompressor());
	if (infmt == 0) {
		inbuf = invcf.rdbuf();
	}
	else if (infmt == 1) {
		zipbuf.push(invcf);
		inbuf = &zipbuf;
	} else if (infmt == 2) {
		inbuf = std::cin.rdbuf();
	}
	std::istream instream(inbuf);

	// process VCF file
	std::string vcfline;
	std::vector<std::string> vcfvec;
	std::vector<std::string>::iterator iter;

	// print headers and set up vector to hold tokens
	while (getline(instream, vcfline)) {
		outvcf << vcfline << "\n";
		if (vcfline[0] != '#' || vcfline[1] != '#')
			break;
	}

	std::stringstream ss(vcfline);
	std::string tok;
	while (ss >> tok) {
		vcfvec.push_back("");
	}
	vcfvec.resize(vcfvec.size());
	unsigned int nind = vcfvec.size()-9;

	std::cerr << "Number individuals in VCF: " << nind << "\n\n";

	// filter sites
	/*
	 * flags: N=No data for site, C=min coverage, D=max coverage, U=min number of individuals with data, A=biallelic,
	 * R=RMS map quality, M=map quality bias, P=position bias, S=strand bias, B=base quality bias, V=variant quality,
	 * Q=QUAL, H=excess heterozygosity, F=Unknown reference allele, G=min number of genotyped individuals
	*/

	int infosize = 20; // can dynamically set this
	std::vector<std::string> infovec;
	infovec.resize(infosize);

	std::stringstream infostream;
	std::string contig = "";

	std::string badflags;
	badflags.reserve(16);
	region keep; // good sites
	vcfrecord vcfinfo(badflags.capacity()); // vcf entry
	std::pair <std::string,unsigned int> prev ("",0);

	int i=0;
	size_t indcounts [2]; // 0=number genotyped individuals, 1=number of individuals with min coverage

	unsigned int pos = 0;
	double maf = 0.0;

	while (getline(instream, vcfline)) {
	// assume vcf fields are [0]=CHROM, [1]=POS, [2]=ID, [3]=REF, [4]=ALT, [5]=QUAL, [6]=FILTER, [7]=INFO, [8]=FORMAT

		// fill up site information vector
		ss.clear();
		ss.str(vcfline);
		iter=vcfvec.begin();
		while (ss >> *iter) {
			++iter;
		}
		if (iter < vcfvec.end()) {
			std::cerr << "Input VCF appears truncated\n";
			return -1;
		}

		// process site
		badflags.clear();
		pos = atoi(vcfvec[1].c_str());

		// check for completely missing data
		if (vcfvec[7] == ".") {
			badflags.push_back('N');
		} else {
			// check QUAL
			if (atoi(vcfvec[5].c_str()) < qual) {
				badflags.push_back('Q');
			}
			// check the type of site
			if (vcfvec[3].size() > 1 || vcfvec[4].size() > 1 || vcfvec[4] == "*") {
				if (biallelic && isMultiSNP(vcfvec)) {
					// multi-allelic
					badflags.push_back('A');
				} else {
					// indel, *=deletion
					badflags.push_back('I');
				}
			} else {
				if (vcfvec[3] == "A" || vcfvec[3] == "C" || vcfvec[3] == "G" || vcfvec[3] == "T") {

					if (allsites || vcfvec[4] != ".") {

						// extract individual information
						if (extractIndInfo(vcfvec, indcounts, minind_cov, verbose)) {
							return -1;
						}

						// check number of genotyped individuals
						//std::cerr << indcounts[0] << "\n"; //debug
						if (indcounts[0] < mingeno) badflags.push_back('G');

						// check number of individuals with data
						//std::cerr << indcounts[1] << "\n"; //debug
						if (indcounts[1] < minind) badflags.push_back('U');

						// examine INFO field
						infostream.str(std::string());
						infostream.clear();
						infostream.str(vcfvec[7]);
						i=0;

						while (std::getline(infostream, infovec[i], ';')) {
							++i;
						}

						checkGatkInfo(infovec, i, &badflags, mincov, maxcov, rms_mapq, mqRankSum, posbias, strandbias,
								baseqbias, varqual_depth, hetexcess);

						// get allele frequency info for potential snp
						if (vcfvec[4] != ".") {
							maf = getMaf(vcfvec, verbose);
						} else {
							maf = 0.0;
						}
					}

				} else {
					// unknown reference allele
					badflags.push_back('F');
				}
			}
		}

		if (!vcfinfo.entry.empty()) {

			// check for multiple entries of position
			if ((vcfinfo.pos == pos && vcfinfo.contig == vcfvec[0]) || (vcfinfo.pos == prev.second && vcfinfo.contig == prev.first)) {
				if (biallelic) vcfinfo.flags.push_back('A');
			}

			if (vcfinfo.flags.empty()) {
				// update good sites
				if(!updateRegion(vcfinfo, keep, passpos)) {
					return -1;
				}

				// output to VCF
				if (allsites_vcf) {
					// print site whether it's variable or not
					outvcf << vcfinfo.entry << "\n";
				} else if ( vcfinfo.alt != '.') {
					// potential snp
					if (!varonly || (varonly && vcfinfo.af > mafcutoff)) {
						outvcf << vcfinfo.entry << "\n";
					}
				}

			} else {
				// process sites that failed QC
				if (!writeBads(vcfinfo.flags, vcfinfo.contig, &vcfinfo.pos, failpos)) {
					return -1;
				}
			}
		}

		prev.first = vcfinfo.contig;
		prev.second = vcfinfo.pos;
		vcfinfo.newEntry(&vcfline, &vcfvec[0], &pos, (vcfvec[4].c_str())[0], maf, &badflags);
	}

	// write last staged sites
	if (vcfinfo.flags.empty()) {
		if (!updateRegion(vcfinfo, keep, passpos)) {
			return -1;
		}

		if (allsites_vcf) {
			// print site whether it's variable or not
			outvcf << vcfinfo.entry << "\n";
		} else if ( vcfinfo.alt != '.') {
			// potential snp
			if (!varonly || (varonly && vcfinfo.af > mafcutoff)) {
				outvcf << vcfinfo.entry << "\n";
			}
		}

	} else {
		if (!writeBads(vcfinfo.flags, vcfinfo.contig, &vcfinfo.pos, failpos)) {
			return -1;
		}
	}


	if (keep.entries > 0) {
		if (!keep.write(passpos)) {
			return -1;
		}
	}

	std::cerr << "Finished processing VCF\n";

	return rv;
}

double getMaf (const std::vector<std::string> &vcfvec, int verbose) {
	// only returns frequency of first-listed ALT allele
	static char f [30];
	f[0] = '\0';
	int j=0;
	static std::string info;
	info.clear();
	info = (vcfvec[7]);
	for (unsigned int i=0; i<info.size(); ++i) {
		if (i+2 < info.size() && info[i]=='A' && info[i+1]=='F' && info[i+2]=='=') {
			i += 3;
			while (i < info.size() && info[i] != ';' && info[i] != ',') {
				f[j] = info[i];
				++j;
				++i;
			}
			f[j] = '\0';
			break;
		}
	}

	if (!f[0]) {
		if (verbose) std::cerr << "WARNING: No AF found in INFO field for " << vcfvec[0] << " " << vcfvec[1] << "\n";
		return -1.0;
	}

	double freq = atof(f);
	return (freq < 0.5 ? freq : 1.0-freq);
}

int isMultiSNP (std::vector<std::string> &vcfvec) {
	unsigned int i;

	for (i=0; i<vcfvec[4].size(); ++i) {
		if (vcfvec[4][i] == ',') return 1;
	}

	for (i=0; i<vcfvec[3].size(); ++i) {
		if (vcfvec[3][i] == ',') return 1;
	}

	return 0;
}

region* updateRegion (vcfrecord &site, region& goodpos, std::ofstream& outstream) {
	if (goodpos.start == 0) {
		goodpos.makeNew(site.contig, site.pos);
	} else if (site.pos - goodpos.end == 1) {
		++goodpos.end;
	} else if (site.pos == goodpos.end && site.contig == goodpos.contig) {
		++goodpos.entries;
	} else if (site.pos - goodpos.end > 1 || site.contig != goodpos.contig) {
			if(!goodpos.write(outstream)) {
				return NULL;
		}
		goodpos.makeNew(site.contig, site.pos);
	} else {
		std::cerr << "Positions in VCF do not appear to be in increasing order\n";
		return NULL;
	}

	return &goodpos;
}

std::ofstream* writeBads (std::string& flags, std::string &contig, const unsigned int* pos, std::ofstream& outstream) {
	if (!outstream.is_open()) {
		std::cerr << "Unable to write qc-failed site information to unopened outstream\n";
		return NULL;
	}

	for (unsigned int i=0; i<flags.size(); ++i) {
		// don't output sites that are indels
		if (flags[i] == 'I') return &outstream;
	}

	//std::cerr << contig << "\t" << *pos << "\t" << flags << "\n"; //debug
	outstream << contig << "\t" << *pos << "\t" << flags << "\n";

	return &outstream;
}

int parseFormat (std::vector<std::string> &vcfvec, int* index) {
	// determine where the GT and DP info is in the FORMAT field
	int rv = 0;
	static const int nflags = 2; // 0=GT, 1=DP
	for (int i=0; i<nflags; ++i) index[i] = 0;

	int gtfound = 0, dpfound = 0;

	for (unsigned int i=0; i<vcfvec[8].size(); ++i) {
		if (vcfvec[8][i] == 'G' && vcfvec[8][i+1] == 'T') {
			gtfound = 1;
		} else if (vcfvec[8][i] == 'D' && vcfvec[8][i+1] == 'P') {
			dpfound = 1;
		}

		if (vcfvec[8][i] == ':') {
			if (!gtfound) ++index[0];
			if (!dpfound) ++index[1];
		}

		if (gtfound + dpfound == nflags) break;
	}

	if (!gtfound) {
		rv = -1;
		std::cerr << "ERROR: No GT found in FORMAT field for " << vcfvec[0] << " " << vcfvec[1] << "\n";
	}

	if (!dpfound) {
		rv = -1;
		std::cerr << "ERROR: No DP found in FORMAT field for " << vcfvec[0] << " " << vcfvec[1] << "\n";
	}

	return rv;
}

int extractIndInfo (std::vector<std::string> &vcfvec, size_t* indcounts, unsigned int min_indcov, int verbose) {
	int rv = 0;
	static const int nflags = 2;

	// Find location of info in format field
	int formatidx [nflags];
	if (parseFormat(vcfvec, formatidx)) {
		return -1;
	}

	unsigned int i = 0, j = 0;
	static std::vector<std::string>::iterator iter;
	static char dp [30];
	static char gt [10];

	for (i=0; i<nflags; ++i) indcounts[i] = 0;
	int c=0;
	int ind=1;

	static std::string missing = ""; // flags which subfields are missing
	static int miss [nflags]; // 0=GT, 1=DP
	for (i=0; i<nflags; ++i) miss[i] = 0;

	for (iter=vcfvec.begin()+9; iter < vcfvec.end(); ++iter) {
		dp[0] = '\0';
		gt[0] = '\0';
		c = 0;
		i = 0;
		int nparsed=0;
		int validgeno=0;

		for (i=0; i<iter->size(); ++i) {
			if (c == formatidx[0]) {
				// get genotype
				j=0;
				while (i < iter->size()) {
					if ((*iter)[i] == ':') break;
					gt[j] = (*iter)[i];
					if (j > 0 && (gt[j] == '/' || gt[j] == '|')) ++validgeno;
					++i;
					++j;
				}
				gt[j] = '\0';
				++c;
				++nparsed;
			} else if (c == formatidx[1]) {
				// get DP
				j=0;
				while (i < iter->size()) {
					if ((*iter)[i] == ':') break;
					dp[j] = (*iter)[i];
					++i;
					++j;
				}
				dp[j] = '\0';
				++c;
				++nparsed;
			} else if ((*iter)[i] == ':') {
				++c;
			}

			if (nparsed >= nflags) break;
		}

		// Check for missing genotype
		// Alleles can be separated by '/' (unphased) or '|' (phased)
		// The first subfield must always be genotype according to VCF documentation
		if (validgeno) {

			// check if genotype is called
			if (strcmp(gt, "./.") != 0 && strcmp(gt, ".|.") != 0) {
				++indcounts[0];
			}

			// check if coverage requirement is met
			if (dp[0]) {
				if (static_cast<unsigned int>(atoi(dp)) >= min_indcov) ++indcounts[1];
			} else {
				if (verbose == 2) {
					std::cerr << "WARNING: Couldn't find DP value for " << vcfvec[0] << " " << vcfvec[1] << " individual " << ind << "\n";
				} else if (verbose == 1 && !miss[1]) {
					missing += " DP";
					miss[1] = 1;
				}
			}

		} else {
			std::cerr << "ERROR: Couldn't find genotype information for " << vcfvec[0] << " " << vcfvec[1] << " individual " << ind << "\n";
			return -1;
		}

		/*
		std::cerr << gt << "\t" << dp << "\n"; // debug
		std::cerr << indcounts[0] << "\t" << indcounts[1] << "\n"; // debug
		 */

		++ind;
	}

	// notify of missing subfields for a site
	if (verbose == 1 && !missing.empty()) {
		std::cerr << "WARNING: Missing" << missing << " values for " << vcfvec[0] << " " << vcfvec[1] << "\n";
	}

	return rv;
}

void checkGatkInfo(std::vector<std::string> &info, int n, std::string* flags, const unsigned int &mincov, const unsigned int &maxcov,
		const double &rms_mapq, const double &mqRankSum, const double &posbias, const double &strandbias, const double &baseqbias,
		const double & varqual_depth, const double &hetexcess) {

	static char id [20];
	static char val [40];
	char* idptr = id;
	double v = 0;
	double absv = 0;

	for (int i=0; i<n; ++i) {

		int k=0;
		idptr = id;
		for (unsigned int j=0; j<info[i].length(); ++j) {
			if (info[i][j] == '=') {
				idptr[k] = '\0';
				idptr = val;
				k=0;
				continue;
			}
			idptr[k] = info[i][j];
			++k;
		}
		idptr[k] = '\0';
		v = atof(val);
		absv = fabs(v);

		if (strcmp(id, "DP") == 0) {
			if (v < mincov) flags->push_back('C');
			if (v > maxcov) flags->push_back('D');
		}

		else if (strcmp(id, "MQ") == 0 && v < rms_mapq) {
			flags->push_back('R');
		}

		else if (strcmp(id, "MQRankSum") == 0 && absv > mqRankSum) {
			flags->push_back('M');
		}

		else if (strcmp(id, "ReadPosRankSum") == 0 && absv > posbias) {
			flags->push_back('P');
		}

		else if (strcmp(id, "FS") == 0 && v > strandbias) {
			flags->push_back('S');
		}

		else if (strcmp(id, "BaseQRankSum") == 0 && fabs(v) > baseqbias) {
			flags->push_back('B');
		}

		else if (strcmp(id, "QD") == 0 && v < varqual_depth) {
			flags->push_back('V');
		}

		else if (strcmp(id, "ExcessHet") == 0 && v > hetexcess) {
			flags->push_back('H');
		}

	}
}

int parseGATKargs (int argc, char** argv, std::ifstream &invcf, std::ofstream &outvcf, std::ofstream &passpos, std::ofstream &failpos,
	int &biallelic, int &allsites, int &allsites_vcf, unsigned int &maxcov, unsigned int &mincov, unsigned int &minind_cov,
	unsigned int &minind, unsigned int &mingeno, double &rms_mapq, double &mqRankSum, double &posbias, double &strandbias, double &baseqbias,
	double& qual, double &varqual_depth, double &hetexcess, int &varonly, double &maf, int &verbose, int &infmt) {

	int rv = 0;

	if (argc < 6) {
		gatkinfo(biallelic, allsites, allsites_vcf, maxcov, mincov, minind_cov, minind, mingeno, rms_mapq, mqRankSum,
				posbias, strandbias, baseqbias, qual, varqual_depth, hetexcess, varonly, maf, verbose);
		return 1;
	}

	int argpos = 2;
	const char* invcf_name = NULL;
	std::string outvcf_name;
	std::string goodpos_name;
	std::string badpos_name;

	while (argpos < argc) {
		if (strcmp(argv[argpos], "-vcf") == 0) {
			// input VCF
			invcf_name = argv[argpos+1];
			if (strcmp(invcf_name, "-") == 0) {
				infmt = 2; // reading from standard input
			} else {
				invcf.open(argv[argpos+1], std::ios_base::in | std::ios_base::binary);
				if (! invcf) {
					std::cerr << "Unable to open input VCF " << argv[argpos+1] << "\n";
					return -1;
				}
				// check for gzipped file by reading magic numbers
				unsigned char magic [2] = {0};
				invcf.read(reinterpret_cast<char*>(magic), sizeof(magic));
				infmt = (magic[0] == 0x1f && magic[1] == 0x8b) ? 1 : 0;
				invcf.seekg(0, std::ios_base::beg);
			}
		}

		else if (strcmp(argv[argpos], "-out") == 0) {
			std::string outprefix(argv[argpos+1]);
			// output VCF
			outvcf_name = outprefix + ".vcf";
			outvcf.open(outvcf_name.c_str(), std::ios::out);
			if (! outvcf) {
				std::cerr << "Unable to open output VCF " << outvcf_name << "\n";
				return -1;
			}
			// output QC pass positions
			goodpos_name = outprefix + "_pass.pos";
			passpos.open(goodpos_name.c_str(), std::ios::out);
			if (! passpos) {
				std::cerr << "Unable to open QC-passed position file " << goodpos_name << "\n";
				return -1;
			}
			// output QC failed positions
			badpos_name = outprefix + "_fail.pos";
			failpos.open(badpos_name.c_str(), std::ios::out);
			if (! failpos) {
				std::cerr << "Unable to open QC-failed position file " << badpos_name << "\n";
				return -1;
			}
		}

		else if (strcmp(argv[argpos],"-biallelic") == 0) {
			// biallelic SNP filter
			biallelic = atoi(argv[argpos+1]);
			switch (biallelic) {
				case 0:
					std::cerr << "retaining multiallelic sites (-biallelic 0) is still in development...\n";
					return 1;
				case 1:
					break;
				default:
					std::cerr << "-biallelic must be 0 to keep all SNPs or 1 to keep only biallelic SNPs\n";
					return -1;
			}
		}

		else if (strcmp(argv[argpos], "-allsites") == 0) {
			// process all input sites
			allsites = atoi(argv[argpos+1]);
			switch (allsites) {
				case 0:
					break;
				case 1:
					break;
				default:
					std::cerr << "-allsites must be 0 to process only SNPs or 1 to process all sites\n";
					return -1;
			}

		}

		else if (strcmp(argv[argpos], "-allsites_vcf") == 0) {
			// output vcf type
			allsites_vcf = atoi(argv[argpos+1]);
			switch (allsites_vcf) {
				case 0:
					break;
				case 1:
					break;
				default:
					std::cerr << "-allsites_vcf must be 0 to have the output VCF contain SNPs only or 1 to contain all sites\n";
					return -1;
			}
		}

		else if (strcmp(argv[argpos], "-maxcov") == 0) {
			// max total site depth
			maxcov = atoi(argv[argpos+1]);
			if (maxcov <= 0) {
				std::cerr << "-maxcov must be > 0\n";
				return -1;
			}
		}

		else if (strcmp(argv[argpos], "-mincov") == 0) {
			// min total site depth
			mincov = atoi(argv[argpos+1]);
			if (mincov < 0) {
				std::cerr << "-mincov must be >= 0\n";
				return -1;
			}
		}

		else if (strcmp(argv[argpos], "-minind_cov") == 0) {
			// min individual depth
			minind_cov = atoi(argv[argpos+1]);
			if (minind_cov < 0) {
				std::cerr << "-minind_cov must be >= 0\n";
				return -1;
			}
		}

		else if (strcmp(argv[argpos],"-minind") == 0) {
			// min number of individuals with data
			minind = atoi(argv[argpos+1]);
			if (minind < 0) {
				std::cerr << "-minind must be >= 0\n";
				return -1;
			}
		}

		else if (strcmp(argv[argpos],"-mingeno") == 0) {
			// min number of individuals with called genotype
			mingeno = atoi(argv[argpos+1]);
			if (mingeno < 0) {
				std::cerr << "-mingeno must be >= 0\n";
				return -1;
			}
		}

		else if (strcmp(argv[argpos], "-rms_mapq") == 0) {
			// min RMS mapping quality
			rms_mapq = atof(argv[argpos+1]);
			if (rms_mapq < 0.0) {
				std::cerr << "-rms_mapq must be >= 0\n";
				return -1;
			}
		}

		else if (strcmp(argv[argpos],"-mapqRankSum") == 0) {
			// max map quality bias
			mqRankSum = atof(argv[argpos+1]);
			if (mqRankSum < 0.0) {
				std::cerr << "-mapqRankSum must be >= 0\n";
				return -1;
			}
		}

		else if (strcmp(argv[argpos], "-posbias") == 0) {
			// max allele position bias
			posbias = atof(argv[argpos+1]);
			if (posbias < 0) {
				std::cerr << "-posbias must be >= 0\n";
				return -1;
			}
		}

		else if (strcmp(argv[argpos], "-strandbias") == 0) {
			// max strand bias
			strandbias = atof(argv[argpos+1]);
			if (strandbias < 0) {
				std::cerr << "-strandbias must be >= 0\n";
				return -1;
			}
		}

		else if (strcmp(argv[argpos], "-baseqbias") == 0) {
			// max base quality bias
			baseqbias = atof(argv[argpos+1]);
			if (baseqbias < 0) {
				std::cerr << "-baseqbias must be >= 0\n";
				return -1;
			}
		}

		else if (strcmp(argv[argpos], "-qual") == 0) {
			// min QUAL field value
			qual = atof(argv[argpos+1]);
			if (qual < 0) {
				std::cerr << "-qual must be >= 0\n";
				return -1;
			}
		}

		else if (strcmp(argv[argpos], "-varq_depth") == 0) {
			// min variant quality divided by depth of non-homo-ref individuals
			varqual_depth = atof(argv[argpos+1]);
			if (varqual_depth < 0) {
				std::cerr << "-varq_depth must be >= 0\n";
				return -1;
			}
		}

		else if (strcmp(argv[argpos],"-hetexcess") == 0) {
			// max excess heterozygosity
			hetexcess = atof(argv[argpos+1]);
			if (hetexcess < 0) {
				std::cerr << "-hetexcess must be >= 0\n";
				return -1;
			}
		}

		else if (strcmp(argv[argpos], "-varonly") == 0) {
			// AF criteria for outputting variable sites
			varonly = atoi(argv[argpos+1]);
			switch (varonly) {
				case 0 :
					break;
				case 1 :
					break;
				default:
					std::cerr << "Invalid -varonly value, (should be 0 or 1)\n";
					return -1;
			}
		}

		else if (strcmp(argv[argpos], "-maf") == 0) {
			// MAF lower bound for variable site VCFs
			maf = atof(argv[argpos+1]);
			if (maf < 0.0 || maf > 0.5) {
				std::cerr << "-maf out of [0, 0.5] range\n";
				return -1;
			}
		}

		else if (strcmp(argv[argpos], "-verbose") == 0) {
			// level of warnings to issue
			verbose = atoi(argv[argpos+1]);
			switch (verbose) {
				case 0 :
					break;
				case 1 :
					break;
				case 2 :
					break;
				default:
					std::cerr << "Invalid -verbose level (should be 0, 1, or 2)\n";
					return -1;
			}
		}

		else {
			std::cerr << "Unknown argument " << argv[argpos] << "\n";
			return -1;
		}

		argpos += 2;
	}

	if (!allsites && allsites_vcf) {
		std::cerr << "Warning: -allsites_vcf is 1 but only SNPs will be printed to VCF because -allsites is 0\n";
	}

	printUserArgs(invcf_name, outvcf_name, goodpos_name, badpos_name, biallelic, allsites, allsites_vcf, maxcov,
			mincov, minind_cov, minind, mingeno, rms_mapq, mqRankSum, posbias, strandbias, baseqbias, qual, varqual_depth, hetexcess,
			varonly, maf);

	return rv;
}

void printUserArgs (const char* invcf_name, std::string &outvcf_name, std::string &goodpos_name, std::string &badpos_name,
	int &biallelic, int &allsites, int &allsites_vcf, unsigned int &maxcov, unsigned int &mincov, unsigned int &minind_cov,
	unsigned int &minind, unsigned int &mingeno, double &rms_mapq, double &mqRankSum, double &posbias, double &strandbias, double &baseqbias,
	double &qual, double &varqual_depth, double &hetexcess, int &varonly, double &maf) {

	int w=20;
	std::cerr << "\n"
	<< std::setw(w) << std::left << "Reading VCF: " << invcf_name << "\n"
	<< std::setw(w) << std::left << "Filtered VCF: " << outvcf_name << "\n"
	<< std::setw(w) << std::left << "Passed positions: " << goodpos_name << "\n"
	<< std::setw(w) << std::left << "Failed positions: " << badpos_name << "\n"
	<< std::setw(w) << std::left << "allsites: " << allsites << "\n"
	<< std::setw(w) << std::left << "allsites_vcf: " << allsites_vcf << "\n"
	<< std::setw(w) << std::left << "biallelic: " << biallelic << "\n"
	<< std::setw(w) << std::left << "maxcov: " << maxcov << "\n"
	<< std::setw(w) << std::left << "mincov: " << mincov << "\n"
	<< std::setw(w) << std::left << "minind_cov: " << minind_cov << "\n"
	<< std::setw(w) << std::left << "minind: " << minind << "\n"
	<< std::setw(w) << std::left << "mingeno: " << mingeno << "\n"
	<< std::setw(w) << std::left << "rms_mapq: " << rms_mapq << "\n"
	<< std::setw(w) << std::left << "mapqRankSum: " << mqRankSum << "\n"
	<< std::setw(w) << std::left << "posbias: " << posbias << "\n"
	<< std::setw(w) << std::left << "strandbias: "  << strandbias << "\n"
	<< std::setw(w) << std::left << "baseqbias: " << baseqbias << "\n"
	<< std::setw(w) << std::left << "qual: " << qual << "\n"
	<< std::setw(w) << std::left << "varq_depth: " << varqual_depth << "\n"
	<< std::setw(w) << std::left << "hetexcess: " << hetexcess << "\n"
	<< std::setw(w) << std::left << "varonly: " << varonly << "\n"
	<< std::setw(w) << std::left << "maf: " << maf << "\n"
	<< "\n";
}

region::region ()
	: start(0),
	  end(0),
	  entries(0)
{}

void region::clear () {
	contig.clear();
	start = 0;
	end = 0;
	entries = 0;
}

void region::makeNew (std::string &seqid, unsigned int pos) {
	contig = seqid;
	start = pos;
	end = pos;
	entries = 1;
}

std::ofstream* region::write (std::ofstream& outstream) {
	if (! outstream.is_open()) {
		std::cerr << "Unable to write qc-passed region information to unopened outstream\n";
		return NULL;
	}
	outstream << contig << ":" << start;
	if (start != end) outstream << "-" << end;
	outstream << "\n";
	return &outstream;
}

vcfrecord::vcfrecord (size_t flagcapacity)
	: entry(""),
	  flags(""),
	  contig(""),
	  pos(0),
	  alt('\0'),
	  af(0)
{
	flags.reserve(flagcapacity);
}

void vcfrecord::newEntry (std::string* vcfline, std::string* id, unsigned int* position,
		char altallele, double allelef, std::string* f) {
	entry = *vcfline;
	contig = *id;
	pos = *position;
	alt = altallele;
	af = allelef;
	if (f) {
		flags = *f;
	} else {
		flags.clear();
	}
}
