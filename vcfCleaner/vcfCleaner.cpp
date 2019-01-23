/*
 * vcfCleaner.cpp
 *
 * TODO:
 * - dynamically set infosize (number of INFO IDs) based on VCF header
 */

#include "vcfCleaner.h"
#include <cstring>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>

int main (int argc, char** argv) {
	int rv = 0;

	std::fstream invcf; // input vcf
	std::fstream outvcf; // filtered output vcf
	std::fstream pass; // quality passed list of sites
	std::fstream fail; // failed sites

	if (argc < 2) {
		maininfo();
	} else if (strcmp(argv[1], "gatk") == 0) {
		gatkvcf (argc, argv, invcf, outvcf, pass, fail);
	} else {
		std::cerr << "\n" << argv[1] << " is an invalid VCF type\n";
		maininfo();
	}

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
	unsigned int &minind, double &rms_mapq, double &mqRankSum, double &posbias, double &strandbias, double &baseqbias,
	double &varqual_depth, double &hetexcess) {

	int w1=20;
	int w2=8;

	std::cerr << "\nThis is the subroutine used to filter VCFs produced by the GATK's HaplotypeCaller/GenotypeGVCFs workflow\n"
	<< "Indels are ignored and not included in files containing filtered sites\n"
	<< "\nvcfCleaner gatk [arguments]\n\n"
	<< std::setw(w1) << std::left << "-vcf" << std::setw(w2) << std::left << "FILE" << "Input VCF file to filter\n"
	<< std::setw(w1) << std::left << "-out" << std::setw(w2) << std::left << "STRING" << "Output file name (prefix)\n"
	<< std::setw(w1) << std::left << "-biallelic" << std::setw(w2) << std::left << "0|1" << "0: keep multiallelic sites, 1: keep only biallelic sites, A [" << biallelic << "]\n"
	<< std::setw(w1) << std::left << "-allsites" << std::setw(w2) << std::left << "0|1" << "0: process only SNPs, 1: process all sites (except indels) [" << allsites << "]\n"
	<< std::setw(w1) << std::left << "-allsites_vcf" << std::setw(w2) << std::left << "0|1" << "Output VCF contains QC-passed 0: SNPs only, 1: all sites [" << allsites_vcf << "]\n"
	<< std::setw(w1) << std::left << "-maxcov" << std::setw(w2) << std::left << "INT" << "Max total site depth, D [" << maxcov << "]\n"
	<< std::setw(w1) << std::left << "-mincov" << std::setw(w2) << std::left << "INT" << "Min total site depth, C [" << mincov << "]\n"
	<< std::setw(w1) << std::left << "-minind_cov" << std::setw(w2) << std::left << "INT" << "Min individual coverage [" << minind_cov << "]\n"
	<< std::setw(w1) << std::left << "-minind" << std::setw(w2) << std::left << "INT" << "Min number of individuals covered by at least -minind_cov reads, U [" << minind << "]\n"
	<< std::setw(w1) << std::left << "-rms_mapq" << std::setw(w2) << std::left << "FLOAT" << "Min RMS mapping quality, R [" << rms_mapq << "]\n"
	<< std::setw(w1) << std::left << "-mapqRankSum" << std::setw(w2) << std::left << "FLOAT" << "Max absolute Wilcoxon rank sum test Z-score of alt vs. ref read map quality, M [" << mqRankSum << "]\n"
	<< std::setw(w1) << std::left << "-posbias" << std::setw(w2) << std::left << "FLOAT" << "Max absolute Wilcoxon rank sum test Z-score of alt vs. ref read position bias, P [" << posbias << "]\n"
	<< std::setw(w1) << std::left << "-strandbias" << std::setw(w2) << std::left << "FLOAT" << "Max Phred-scaled Fisher's exact test p-value of strand bias, S [" << strandbias << "]\n"
	<< std::setw(w1) << std::left << "-baseqbias" << std::setw(w2) << std::left << "FLOAT" << "Max absolute Wilcoxon rank sum test Z-score of alt vs ref base qualities, B [" << baseqbias << "]\n"
	<< std::setw(w1) << std::left << "-varq_depth" << std::setw(w2) << std::left << "FLOAT" << "Min variant Phred-scaled confidence/quality by depth, V [" << varqual_depth << "]\n"
	<< std::setw(w1) << std::left << "-hetexcess" << std::setw(w2) << std::left << "FLOAT" << "Max Phred-scaled p-value for exact test of excess heterozygosity, H [" << hetexcess << "]\n"
	<< "\nOther site QC fail flags:\n"
	<< "F, Unkown reference allele\n"
	<< "\n";
}

int gatkvcf (int argc, char** argv, std::fstream &invcf, std::fstream &outvcf, std::fstream &passpos, std::fstream &failpos) {
	int rv = 0;

	// filtering parameters
	int biallelic = 1; // 0: keep multiallelic sites, 1: keep only biallelic sites
	int allsites = 0; // 0: process only SNPs, 1: process all sites (except indels)
	int allsites_vcf = 0; // output VCF contains QC-passed 0: SNPs only or 1: all sites
	unsigned int maxcov = 1000000; // maximum total site coverage (DP)
	unsigned int mincov = 2; // minimum total site coverage (DP)
	unsigned int minind_cov = 1; // minimum individual coverage (DP)
	unsigned int minind = 1; // minimum number of individuals with data
	double rms_mapq = 40.0; // min RMS mapping quality (MQ)
	double mqRankSum = 12.5; // max absolute Wilcoxon rank sum test z-score of alt vs ref read map quality (MQRankSum)
	double posbias = 60.0; // max absolute Wilcoxon rank sum test z-score of alt vs ref read position bias  (ReadPosRankSum)
	double strandbias = 8.0; // max Phred-scale p-value using Fisher's exact test of strand bias (FS)
	double baseqbias = 21.0; // max absolute Wilcoxon rank sum test z-score of alt vs. ref base qualities (BaseQRankSum)
	double varqual_depth = 2.0; // min variant Phred-scaled confidence/quality by depth (QD)
	double hetexcess = 40.0; // max Phred-scaled p-value for exact test of excess heterozygosity (ExcessHet)

	if ((rv=parseGATKargs(argc, argv, invcf, outvcf, passpos, failpos, biallelic, allsites, allsites_vcf, maxcov, mincov, minind_cov,
			minind, rms_mapq, mqRankSum, posbias, strandbias, baseqbias, varqual_depth, hetexcess))) {
		if (rv > 0)
			return 0;
		else if (rv < 0)
			return -1;
	}

	// process VCF file
	std::string vcfline;
	std::vector<std::string> vcfvec;
	std::vector<std::string>::iterator iter;

	// skip headers and set up vector to hold tokens
	getline(invcf, vcfline);
	while (!vcfline.empty()) {
		outvcf << vcfline << "\n";
		if (vcfline[0] != '#' || vcfline[1] != '#')
			break;
		getline(invcf, vcfline);
	}
	std::stringstream ss(vcfline);
	std::string tok;
	while (ss >> tok) {
		vcfvec.push_back(tok);
	}
	vcfvec.resize(vcfvec.size());
	unsigned int nind = vcfvec.size()-9;

	std::cerr << "Number individuals in VCF: " << nind << "\n\n";

	// filter sites
	/*
	 * flags: N=No data for site, C=min coverage, D=max coverage, U=min number of individuals with data, A=biallelic,
	 * R=RMS map quality, M=map quality bias, P=position bias, S=strand bias, B=base quality bias, V=variant quality,
	 * H=excess heterozygosity, F=Unknown reference allele
	*/

	int infosize = 20; // can dynamically set this
	std::vector<std::string> infovec;
	infovec.resize(infosize);

	std::stringstream infostream;
	std::string contig = "";

	std::string badflags;
	badflags.reserve(14);
	region keep; // good sites
	vcfrecord vcfinfo(badflags.capacity()); // vcf entry
	std::pair <std::string,unsigned int> prev ("",0);

	int i=0;
	size_t ncov = 0;
	unsigned int pos = 0;

	while (!invcf.eof()) {
	// assume vcf fields are [0]=CHROM, [1]=POS, [2]=ID, [3]=REF, [4]=ALT, [5]=QUAL, [6]=FILTER, [7]=INFO, [8]=FORMAT

		// get a new site
		getline(invcf,vcfline);
		if (vcfline.empty()) break;
		ss.str(std::string());
		ss.clear();
		ss.str(vcfline);
		iter = vcfvec.begin();
		while (ss >> *iter) {
			++iter;
		}
		badflags.clear();
		pos = atoi(vcfvec[1].c_str());

		// check for completely missing data
		if (vcfvec[7] == ".") {
			badflags.push_back('N');

		} else {
			if (vcfvec[3].size() == 1 && vcfvec[4].size() == 1) {
				if (vcfvec[3] == "A" || vcfvec[3] == "C" || vcfvec[3] == "G" || vcfvec[3] == "T") {
					if (allsites || vcfvec[4] != ".") {

						// determine number of individuals with data
						ncov = nCoveredInd(vcfvec, minind_cov, &rv);
						if (rv) {
							break;
						} else {
							if (ncov < minind) badflags.push_back('U');
						}

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
					}
				} else {
					// unknown reference allele
					badflags.push_back('F');
				}
			} else {
				if (isMultiSNP(vcfvec)) {
					// multi-allelic
					badflags.push_back('A');
				} else {
					// indel
					badflags.push_back('I');
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
					rv = -1;
					break;
				}

				// output to VCF
				if (allsites_vcf || vcfinfo.alt != '.') {
					outvcf << vcfinfo.entry << "\n";
				}

			} else {
				// process sites that failed QC
				if (!writeBads(vcfinfo.flags, vcfinfo.contig, &vcfinfo.pos, failpos)) {
					rv = -1;
					break;
				}
			}
		}

		prev.first = vcfinfo.contig;
		prev.second = vcfinfo.pos;
		vcfinfo.newEntry(&vcfline, &vcfvec[0], &pos, (vcfvec[4].c_str())[0], &badflags);
	}

	// write last staged sites
	if (vcfinfo.flags.empty()) {
		if (!updateRegion(vcfinfo, keep, passpos)) {
			rv = -1;
		}

		if (allsites_vcf || vcfinfo.alt != '.') {
			outvcf << vcfinfo.entry << "\n";
		}

	} else {
		if (!writeBads(vcfinfo.flags, vcfinfo.contig, &vcfinfo.pos, failpos)) {
			rv = -1;
		}
	}

	if (!keep.write(passpos)) {
		rv = -1;
	}

	// close file streams
	if (invcf.is_open()) invcf.close();
	if (outvcf.is_open()) outvcf.close();
	if (passpos.is_open()) passpos.close();
	if (failpos.is_open()) failpos.close();

	std::cerr << "Finished processing VCF\n";

	return rv;
}

int isMultiSNP (std::vector<std::string> &vcfvec) {
	unsigned int i;

	for (i=0; i<vcfvec[4].size(); ++i) {
		if (vcfvec[4][i] == ':') return 1;
	}

	for (i=0; i<vcfvec[3].size(); ++i) {
		if (vcfvec[3][i] == ':') return 1;
	}

	return 0;
}

region* updateRegion (vcfrecord &site, region& goodpos, std::fstream& outstream) {
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

	/*
	 if (goodpos.start == 0) {
		goodpos.makeNew(site.contig, site.pos);
	} else if (site.pos - goodpos.end == 1) {
		++goodpos.end;
	} else if (site.pos - goodpos.end > 1 || site.contig != goodpos.contig) {
		if (biallelic && goodpos.entries > 1) --goodpos.end;
		if (goodpos.end >= goodpos.start) {
			if(!goodpos.write(outstream)) return NULL;
		}
		goodpos.makeNew(site.contig, site.pos);
	} else if (site.pos == goodpos.end && site.contig == goodpos.contig) {
		++goodpos.entries;
		if (biallelic && site.flags.find('A') == std::string::npos) site.flags.push_back('A');
	} else {
		std::cerr << "Positions in VCF do not appear to be in increasing order\n";
		return NULL;
	}
	 *
	 */

	return &goodpos;
}

std::fstream* writeBads (std::string& flags, std::string &contig, const unsigned int* pos, std::fstream& outstream) {
	if (!outstream.is_open()) {
		std::cerr << "Unable to write qc-failed site information to unopened outstream\n";
		return NULL;
	}

	for (unsigned int i=0; i<flags.size(); ++i) {
		// don't output sites that are indels
		if (flags[i] == 'I') return &outstream;
	}

	outstream << contig << "\t" << *pos << "\t" << flags << "\n";

	return &outstream;
}

size_t nCoveredInd (std::vector<std::string> &vcfvec, unsigned int min_indcov, int* rv) {
	// determine where the DP info is in the FORMAT field
	size_t n = 0;
	int dpidx = 0;
	int dpfound = 0;
	unsigned int i = 0;
	for (i=0; i<vcfvec[8].size(); ++i) {
		if (vcfvec[8][i] == 'D' && vcfvec[8][i+1] == 'P') {
			dpfound = 1;
			break;
		}
		if (vcfvec[8][i] == ':') ++dpidx;
	}

	if (!dpfound) {
		*rv = -1;
		std::cerr << "No DP found in FORMAT field for " << vcfvec[0] << " " << vcfvec[1] << "\n";
		return 0;
	}

	// count number of covered individuals
	static std::vector<std::string>::iterator iter;
	static char dp [10];
	int c=0, j=0;

	for (iter=vcfvec.begin()+9; iter < vcfvec.end(); ++iter) {
		// find DP start
		c = 0;
		i = 0;
		while (c < dpidx) {
			++i;
			if ((*iter)[i] == ':') {
				++c;
				++i;
			}
		}

		// get dp value
		j=0;
		while (i < iter->size() && (*iter)[i] != ':') {
			dp[j] = (*iter)[i];
			++i;
			++j;
		}
		dp[j] = '\0';

		// check if coverage requirement is met
		if (static_cast<unsigned int>(atoi(dp)) >= min_indcov) ++n;
	}

	return n;
}

void checkGatkInfo(std::vector<std::string> &info, int n, std::string* flags, const unsigned int &mincov, const unsigned int &maxcov,
		const double &rms_mapq, const double &mqRankSum, const double &posbias, const double &strandbias, const double &baseqbias,
		const double & varqual_depth, const double &hetexcess) {

	static char id [20];
	static char val [20];
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

int parseGATKargs (int argc, char** argv, std::fstream &invcf, std::fstream &outvcf, std::fstream &passpos, std::fstream &failpos,
	int &biallelic, int &allsites, int &allsites_vcf, unsigned int &maxcov, unsigned int &mincov, unsigned int &minind_cov,
	unsigned int &minind, double &rms_mapq, double &mqRankSum, double &posbias, double &strandbias, double &baseqbias,
	double &varqual_depth, double &hetexcess) {

	int rv = 0;

	if (argc < 6) {
		gatkinfo(biallelic, allsites, allsites_vcf, maxcov, mincov, minind_cov, minind, rms_mapq, mqRankSum,
				posbias, strandbias, baseqbias, varqual_depth, hetexcess);
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
			invcf.open(argv[argpos+1], std::ios::in);
			if (! invcf) {
				std::cerr << "Unable to open input VCF " << argv[argpos+1] << "\n";
				return -1;
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
					break;
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
			mincov, minind_cov, minind, rms_mapq, mqRankSum, posbias, strandbias, baseqbias, varqual_depth, hetexcess);

	return rv;
}

void printUserArgs (const char* invcf_name, std::string &outvcf_name, std::string &goodpos_name, std::string &badpos_name,
	int &biallelic, int &allsites, int &allsites_vcf, unsigned int &maxcov, unsigned int &mincov, unsigned int &minind_cov,
	unsigned int &minind, double &rms_mapq, double &mqRankSum, double &posbias, double &strandbias, double &baseqbias,
	double &varqual_depth, double &hetexcess) {

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
	<< std::setw(w) << std::left << "rms_mapq: " << rms_mapq << "\n"
	<< std::setw(w) << std::left << "mapqRankSum: " << mqRankSum << "\n"
	<< std::setw(w) << std::left << "posbias: " << posbias << "\n"
	<< std::setw(w) << std::left << "strandbias: "  << strandbias << "\n"
	<< std::setw(w) << std::left << "baseqbias: " << baseqbias << "\n"
	<< std::setw(w) << std::left << "varq_depth: " << varqual_depth << "\n"
	<< std::setw(w) << std::left << "hetexcess: " << hetexcess << "\n"
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

std::fstream* region::write (std::fstream& outstream) {
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
	  alt('\0')
{
	flags.reserve(flagcapacity);
}

void vcfrecord::newEntry (std::string* vcfline, std::string* id, unsigned int* position, char altallele, std::string* f) {
	entry = *vcfline;
	contig = *id;
	pos = *position;
	alt = altallele;
	if (f) {
		flags = *f;
	} else {
		flags.clear();
	}
}
