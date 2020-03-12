/*
 * vcfCleaner.cpp
 *
 * TODO:
 * 1) dynamically set infosize (number of INFO IDs) based on VCF header
 * 2) implement accurate multiallelic filtering, e.g. maf filtering
 * 3) Boost does not play well with bgzipped files (segfaults or leaves off the last line), so need to
 * switch over to using htslib to parse VCF. Gzipped files are fine.
 * 4) Need to update VCF header to reflect FILTER annotations
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
	double &varqual_depth, double &hetexcess, int &varonly, double &maf, int &rmvIndels, int &siteOnly, int &printFilter, int &verbose) {

	int w1=20;
	int w2=8;

	std::cerr << "\nThis program is for filtering VCFs produced by GATK's HaplotypeCaller/GenotypeGVCFs workflow.\n"
	<< "Takes uncompressed and gzipped/bgzipped VCFs, however bgzipped reading can be buggy so use at your own risk.\n"
	<< "For multiallelic sites only overall site depth, map quality, and variant quality filters are supported.\n"
	<< "\nvcfCleaner gatk [arguments]\n\n"
	<< std::setw(w1) << std::left << "-vcf" << std::setw(w2) << std::left << "FILE" << "Input VCF file to filter ('-' for STDIN)\n"
	<< std::setw(w1) << std::left << "-out" << std::setw(w2) << std::left << "STRING" << "Output file name (prefix)\n"
	<< std::setw(w1) << std::left << "-biallelic" << std::setw(w2) << std::left << "" << "Discard multiallelic variants\n"
	<< std::setw(w1) << std::left << "-allsites" << std::setw(w2) << std::left << "" << "Process all sites (variable and nonvariable)\n"
	<< std::setw(w1) << std::left << "-allsites_vcf" << std::setw(w2) << std::left << "" << "Output VCF contains both variable and nonvariable sites\n"
	<< std::setw(w1) << std::left << "-siteOnly" << std::setw(w2) << std::left << "" << "Do not print individual information to output VCF\n"
	<< std::setw(w1) << std::left << "-printFilter" << std::setw(w2) << std::left << "0|1|2" << "0: Discard QC-failed sites. 1: Overwrite, or, 2: append to FILTER field [" << printFilter << "]\n"
	<< std::setw(w1) << std::left << "-maxcov" << std::setw(w2) << std::left << "INT" << "Max total site depth [" << maxcov << "]\n"
	<< std::setw(w1) << std::left << "-mincov" << std::setw(w2) << std::left << "INT" << "Min total site depth [" << mincov << "]\n"
	<< std::setw(w1) << std::left << "-minind_cov" << std::setw(w2) << std::left << "INT" << "Min individual coverage [" << minind_cov << "]\n"
	<< std::setw(w1) << std::left << "-minind" << std::setw(w2) << std::left << "INT" << "Min number of individuals covered by at least -minind_cov reads [" << minind << "]\n"
	<< std::setw(w1) << std::left << "-mingeno" << std::setw(w2) << std::left << "INT" << "Min number of individuals with called genotype [" << mingeno << "]\n"
	<< std::setw(w1) << std::left << "-rms_mapq" << std::setw(w2) << std::left << "FLOAT" << "Min RMS mapping quality [" << rms_mapq << "]\n"
	<< std::setw(w1) << std::left << "-mapqRankSum" << std::setw(w2) << std::left << "FLOAT" << "Max absolute Wilcoxon rank sum test Z-score of alt vs. ref read map quality [" << mqRankSum << "]\n"
	<< std::setw(w1) << std::left << "-posbias" << std::setw(w2) << std::left << "FLOAT" << "Max absolute Wilcoxon rank sum test Z-score of alt vs. ref read position bias [" << posbias << "]\n"
	<< std::setw(w1) << std::left << "-strandbias" << std::setw(w2) << std::left << "FLOAT" << "Max Phred-scaled Fisher's exact test p-value of strand bias [" << strandbias << "]\n"
	<< std::setw(w1) << std::left << "-baseqbias" << std::setw(w2) << std::left << "FLOAT" << "Max absolute Wilcoxon rank sum test Z-score of alt vs ref base qualities [" << baseqbias << "]\n"
	<< std::setw(w1) << std::left << "-qual" << std::setw(w2) << std::left << "FLOAT" << "Min Phred-scaled quality score of ALT assertion [" << qual << "]\n"
	<< std::setw(w1) << std::left << "-varq_depth" << std::setw(w2) << std::left << "FLOAT" << "Min variant Phred-scaled confidence/quality by depth [" << varqual_depth << "]\n"
	<< std::setw(w1) << std::left << "-hetexcess" << std::setw(w2) << std::left << "FLOAT" << "Max Phred-scaled p-value for exact test of excess heterozygosity [" << hetexcess << "]\n"
	<< std::setw(w1) << std::left << "-rmvIndels" << std::setw(w2) << std::left << "" << "Discard indels (still count towards multiallelic sites)\n"
	<< std::setw(w1) << std::left << "-varonly" << std::setw(w2) << std::left << "" << "The INFO AF for SNP-only output VCF must be in range (0,1) (as opposed to [0,1])\n"
	<< std::setw(w1) << std::left << "-maf" << std::setw(w2) << std::left << "FLOAT" << "Minor allele frequency lower bound for SNP-only VCF [" << maf << "]\n"
	<< std::setw(w1) << std::left << "-verbose" << std::setw(w2) << std::left << "0|1|2" << "Level of warnings to issue: 0 = suppress all, 1 = site-level, 2 = individual-level [" << verbose << "]\n"
	<< "\n";
}

int gatkvcf (int argc, char** argv, std::ifstream &invcf, std::ofstream &outvcf, std::ofstream &passpos, std::ofstream &failpos) {
	int rv = 0;

	// filtering parameters
	int biallelic = 0; // 0: keep multiallelic sites, 1: keep only biallelic sites
	int allsites = 0; // 0: process only SNPs, 1: process all sites
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
	int varonly = 0; // The INFO for SNP-only VCF must have AF in 1=range (0,1), or 0=range [0,1]
	double mafcutoff = 0.0; // minor allele frequency lower bound for SNP-only VCF
	int verbose = 2; // amount of warnings outputted, 0=none, 1=site level, 2=individual level
	int infmt = 0; // 0=uncompressed vcf, 1=gzipped input, 2=uncompressed standard input
	int rmvIndels = 0; // 0 = process indels, 1 = discard sites with indels (indels still count towards identifying multiallelic sites)
	int siteOnly = 0; // 0 = individual information is printed, 1 = Do not print individual info to VCF
	int printFilter = 0; // 0 = sites that fail QC are not printed to VCF, 1 = overwrite FILTER field info with failed filters, 2 = append failed filters to INFO field

	if ((rv=parseGATKargs(argc, argv, invcf, outvcf, passpos, failpos, biallelic, allsites, allsites_vcf, maxcov, mincov, minind_cov,
			minind, mingeno, rms_mapq, mqRankSum, posbias, strandbias, baseqbias, qual, varqual_depth, hetexcess, varonly, mafcutoff,
			rmvIndels, siteOnly, printFilter, verbose, infmt))) {
		if (rv > 0)
			return 0;
		else if (rv < 0)
			return -1;
	}

	// determine if individual info needs to be checked
	int parseIndividuals = (minind == 0 && mingeno == 0) ? 0 : 1;

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
	vcfrecord vcftoke1; // vcf entry for odd lines
	vcfrecord vcftoke2; // vcf entry for even lines
	vcfrecord* vcfptr = &vcftoke1;
	vcfrecord* prev = NULL;
	std::string vcfline;

	// print headers and set number of vcf fields
	while (getline(instream, vcfline)) {
		outvcf << vcfline << "\n";
		if (vcfline[0] != '#' || vcfline[1] != '#')
			break;
	}

	std::stringstream ss(vcfline);
	std::string tok;
	while (ss >> tok) {
		vcftoke1.fields.push_back("");
	}
	vcftoke1.fields.resize(vcftoke1.fields.size());
	vcftoke2.fields.resize(vcftoke1.fields.size());
	unsigned int nind = vcftoke1.fields.size()-9;

	std::cerr << "Number individuals in VCF: " << nind << "\n\n";

	// filter sites
	/*
	 * flags: N=No data for site, C=min coverage, D=max coverage, U=min number of individuals with data, A=biallelic,
	 * R=RMS map quality, M=map quality bias, P=position bias, S=strand bias, B=base quality bias, V=variant quality,
	 * Q=QUAL, H=excess heterozygosity, F=Unknown reference allele, G=min number of genotyped individuals
	*/

	unsigned int nflags = 17;
	vcftoke1.flags.resize(nflags);
	vcftoke2.flags.resize(nflags);

	region keep; // good sites (stores info to get printed to region file)

	size_t indcounts [2]; // 0=number genotyped individuals, 1=number of individuals with min coverage
	int isMultiAllelic = 0;
	int isIndel = 0;
	unsigned int nlines = 0; // number of lines processed

	while (getline(instream, vcfline)) {
	// assume vcf fields are [0]=CHROM, [1]=POS, [2]=ID, [3]=REF, [4]=ALT, [5]=QUAL, [6]=FILTER, [7]=INFO, [8]=FORMAT
		++nlines;
		if (nlines % 2 == 0) {
			// processing even line
			vcfptr = &vcftoke2;
		} else {
			// processing odd line
			vcfptr = &vcftoke1;
		}

		// fill up site information vector
		ss.clear();
		ss.str(vcfline);
		std::vector<std::string>::iterator iter=vcfptr->fields.begin();
		while (ss >> *iter) {
			++iter;
		}
		if (iter < vcfptr->fields.end()) {
			std::cerr << "Input VCF line " << nlines << " appears truncated\n";
			return -1;
		}
		vcfptr->setSite();
		vcfptr->clearFlags();

		// process site

		if (vcfptr->fields[7] == ".") {
			// check for completely missing data
			vcfptr->setFlag(0);
		} else if ( vcfptr->ref.size() == 1 && (vcfptr->ref=="A" || vcfptr->ref=="C" || vcfptr->ref=="G" || vcfptr->ref=="T" || vcfptr->ref=="N") == 0 ) {
			// check for unknown reference allele
			vcfptr->setFlag(1);
		} else  {
			if (allsites || vcfptr->alt != ".") {
				siteType(vcfptr->fields, &isMultiAllelic, &isIndel);

				// check for multiallelic sites
				if (biallelic && isMultiAllelic) {
					vcfptr->setFlag(2);
				}

				// check for indels
				if (rmvIndels && isIndel) {
					vcfptr->setFlag(3);
				}

				// check QUAL
				if (atoi(vcfptr->fields[5].c_str()) < qual) {
					vcfptr->setFlag(4);
				}

				if (parseIndividuals) {
					// extract individual information
					if (extractIndInfo(vcfptr->fields, indcounts, minind_cov, verbose)) {
						return -1;
					}

					// check number of genotyped individuals
					//std::cerr << indcounts[0] << "\n"; //debug
					if (indcounts[0] < mingeno) vcfptr->setFlag(8);

					// check number of individuals with data
					//std::cerr << indcounts[1] << "\n"; //debug
					if (indcounts[1] < minind) vcfptr->setFlag(7);
				}

				// examine INFO field
				checkGatkInfo(vcfptr, mincov, maxcov, rms_mapq, mqRankSum, posbias, strandbias, baseqbias, varqual_depth, hetexcess);

				// get allele frequency info for potential snp
				if (vcfptr->alt != ".") {
					vcfptr->af = getMaf(vcfptr->fields, verbose);
				} else {
					vcfptr->af = 0.0;
				}
			}

		}

		// write output
		if (prev) {
			// check for multiple entries of position
			if (vcfptr->pos == prev->pos && vcfptr->contig == prev->contig) {
				if (biallelic) prev->setFlag(2);
			}

			if(recordSite (prev, keep, outvcf, passpos, failpos, allsites_vcf, printFilter, siteOnly, varonly, mafcutoff, allsites)) {
				return -1;
			}

		}

		prev = vcfptr;
	}

	// check for multiple entries of last position
	if (vcfptr->pos == prev->pos && vcfptr->contig == prev->contig) {
		if (biallelic) vcfptr->setFlag(2);
	}

	// write last staged sites
	if(recordSite (vcfptr, keep, outvcf, passpos, failpos, allsites_vcf, printFilter, siteOnly, varonly, mafcutoff, allsites)) {
		return -1;
	}

	if (keep.entries > 0) {
		if (!keep.write(passpos)) {
			return -1;
		}
	}

	std::cerr << "Finished processing VCF\n";

	return rv;
}

int recordSite (vcfrecord* site, region &keep, std::ofstream &vcfstream, std::ofstream &goodstream, std::ofstream &failstream,
		const int &allsites_vcf, const int &printFilter, const int &siteOnly, const int &varonly, const double &mafcutoff, const int &allsites) {
	// update passed and failed sites output
	if (allsites || (site->alt != "." && (!varonly || (varonly && site->af > mafcutoff)))) {
		if (site->qcfail == 0) {
			// update good sites
			if (!updateRegion(site, keep, goodstream)) {
				return -1;
			}
		} else if (!printFilter) {
			// write site that fails QC (only needed if annotating VCF FILTERS field)
			if (!writeBads(site, failstream)) {
				return -1;
			}
		}
	}

	// write to output VCF
	if (allsites_vcf) {
		// output all sites VCF
		if (printFilter) {
			// write VCF with FILTER annotation
			writeVcf(site, vcfstream, siteOnly, printFilter);
		} else if (site->qcfail == 0) {
			// write VCF only for sites that pass all filters
			writeVcf(site, vcfstream, siteOnly, printFilter);
		}
	} else if (site->alt != ".") {
		// output variant VCF
		if (!varonly || (varonly && site->af > mafcutoff) ) {
			// check whether MAF condition is met
			if (printFilter) {
				// variant VCF with FILTER annotation
				writeVcf(site, vcfstream, siteOnly, printFilter);
			} else if (site->qcfail == 0) {
				// write only variants sites that pass filters
				writeVcf(site, vcfstream, siteOnly, printFilter);
			}
		}
	}

	return 0;
}

void writeVcf (vcfrecord* vcfrec, std::ofstream &os, const int &siteOnly, const int &printFilter) {
	unsigned int i = 0;
	static std::vector< std::pair<std::string,int> >::iterator it;
	for (i = 0; i<6; ++i) os << vcfrec->fields[i] << "\t";

	// print FILTER field
	if (vcfrec->qcfail == 0 && vcfrec->fields[6] == ".") {
		os << "PASS\t";
	} else {
		if (printFilter == 2 && vcfrec->fields[6] != ".") os << vcfrec->fields[6] << ",";
		int j = 0;
		for (it = vcfrec->flags.begin(); it != vcfrec->flags.end(); ++it) {
			if (it->second) {
				os << it->first;
				++j;
				if (j < vcfrec->qcfail) {
					os << ",";
				} else {
					os << "\t";
					break;
				}
			}
		}
	}

	// print INFO
	os << vcfrec->fields[7];

	// print individual genotype information
	if (!siteOnly) {
		for (i = 8; i < vcfrec->fields.size(); ++i) {
			os << "\t" << vcfrec->fields[i];
		}
	}

	os << "\n";
}

double getMaf (const std::vector<std::string> &vcfvec, int verbose) {
	// returns largest MAF (when there are multiple alleles)
	static char f [30];
	double maf = 0.0;
	double maxmaf = 0.0;
	double freq;
	f[0] = '\0';
	int j=0;
	static std::string info;
	info.clear();
	info = vcfvec[7];
	for (unsigned int i=0; i<info.size(); ++i) {
		if (i+2 < info.size() && info[i]=='A' && info[i+1]=='F' && info[i+2]=='=') {
			i += 3;
			while (i < info.size() && info[i] != ';') {
				f[j] = info[i];
				++j;
				++i;
				if (info[i] == ',') {
					f[j] = '\0';
					freq = atof(f);
					maf = (freq < 0.5) ? freq : 1.0-freq;
					if (maf > maxmaf) maxmaf = maf;
					++i;
					j = 0;
				}
			}
			f[j] = '\0';
			freq = atof(f);
			maf = (freq < 0.5) ? freq : 1.0-freq;
			if (maf > maxmaf) maxmaf = maf;
			break;
		}
	}

	if (!f[0]) {
		if (verbose) std::cerr << "WARNING: No AF found in INFO field for " << vcfvec[0] << " " << vcfvec[1] << "\n";
		return -1.0;
	}

	return maxmaf;
}

void siteType (std::vector<std::string> & vcfvec, int* isMultiAllelic, int* isIndel) {

	*isMultiAllelic = 0;
	*isIndel = 0;

	if (vcfvec[3].size() > 1 && vcfvec[4] != ".") {
		*isIndel = 1;
	}

	if (vcfvec[4] == "*") {
		*isIndel = 1; // deletion
	}

	if (vcfvec[4].size() > 1) {
		static std::string allele;
		allele.clear();
		for (unsigned int i=0; i<vcfvec[4].size(); ++i) {
			if (vcfvec[4][i] == ',') {
				*isMultiAllelic = 1;
				if (allele.size() > 1 || allele == "*") {
					*isIndel = 1;
				}
				allele.clear();
			} else {
				allele.push_back(vcfvec[4][i]);
			}
		}
	}
}

region* updateRegion (vcfrecord* site, region& goodpos, std::ofstream& outstream) {
	if (goodpos.start == 0) {
		goodpos.makeNew(site->contig, site->pos);
	} else if (site->pos - goodpos.end == 1) {
		++goodpos.end;
	} else if (site->pos == goodpos.end && site->contig == goodpos.contig) {
		++goodpos.entries;
	} else if (site->pos - goodpos.end > 1 || site->contig != goodpos.contig) {
			if(!goodpos.write(outstream)) {
				return NULL;
		}
		goodpos.makeNew(site->contig, site->pos);
	} else {
		std::cerr << "Positions in VCF do not appear to be in increasing order\n";
		return NULL;
	}

	return &goodpos;
}

std::ofstream* writeBads (vcfrecord* vcfrec, std::ofstream& outstream) {
	if (!outstream.is_open()) {
		std::cerr << "Unable to write qc-failed site information to unopened outstream\n";
		return NULL;
	}

	outstream << vcfrec->contig << "\t" << vcfrec->pos << "\t";

	// print FILTER field
	static std::vector< std::pair<std::string,int> >::iterator it;
	int j = 0;
	for (it = vcfrec->flags.begin(); it != vcfrec->flags.end(); ++it) {
		if (it->second) {
			outstream << it->first;
			++j;
			if (j < vcfrec->qcfail) {
				outstream << ",";
			} else {
				outstream << "\n";
				break;
			}
		}
	}

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

void checkGatkInfo(vcfrecord* vcfrec, const unsigned int &mincov, const unsigned int &maxcov,
		const double &rms_mapq, const double &mqRankSum, const double &posbias, const double &strandbias, const double &baseqbias,
		const double & varqual_depth, const double &hetexcess) {

	const int infosize = 20; // can dynamically set this
	static std::vector<std::string> info; // stores tokenized VCF INFO subfields
	if (info.size() == 0) info.resize(infosize);

	static std::stringstream infostream;
	infostream.clear();
	infostream.str(vcfrec->fields[7]);

	int n = 0;
	while (std::getline(infostream, info[n], ';')) {
		++n;
	}

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
			if (v < mincov) vcfrec->setFlag(5);
			if (v > maxcov) vcfrec->setFlag(6);
		}

		else if (strcmp(id, "MQ") == 0 && v < rms_mapq) {
			vcfrec->setFlag(9);
		}

		else if (strcmp(id, "MQRankSum") == 0 && absv > mqRankSum) {
			vcfrec->setFlag(10);
		}

		else if (strcmp(id, "ReadPosRankSum") == 0 && absv > posbias) {
			vcfrec->setFlag(11);
		}

		else if (strcmp(id, "FS") == 0 && v > strandbias) {
			vcfrec->setFlag(12);
		}

		else if (strcmp(id, "BaseQRankSum") == 0 && absv > baseqbias) {
			vcfrec->setFlag(13);
		}

		else if (strcmp(id, "QD") == 0 && v < varqual_depth) {
			vcfrec->setFlag(14);
		}

		else if (strcmp(id, "ExcessHet") == 0 && v > hetexcess) {
			vcfrec->setFlag(15);
		}

	}
}

int parseGATKargs (int argc, char** argv, std::ifstream &invcf, std::ofstream &outvcf, std::ofstream &passpos, std::ofstream &failpos,
	int &biallelic, int &allsites, int &allsites_vcf, unsigned int &maxcov, unsigned int &mincov, unsigned int &minind_cov,
	unsigned int &minind, unsigned int &mingeno, double &rms_mapq, double &mqRankSum, double &posbias, double &strandbias, double &baseqbias,
	double& qual, double &varqual_depth, double &hetexcess, int &varonly, double &maf, int &rmvIndels, int &siteOnly, int &printFilter, int &verbose, int &infmt) {

	int rv = 0;

	if (argc < 6) {
		gatkinfo(biallelic, allsites, allsites_vcf, maxcov, mincov, minind_cov, minind, mingeno, rms_mapq, mqRankSum,
				posbias, strandbias, baseqbias, qual, varqual_depth, hetexcess, varonly, maf, rmvIndels, siteOnly, printFilter, verbose);
		return 1;
	}

	int argpos = 2;
	const char* invcf_name = NULL;
	std::string outvcf_name;
	std::string goodpos_name;
	std::string badpos_name;
	std::string outprefix;

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
			++argpos;
		}

		else if (strcmp(argv[argpos], "-out") == 0) {
			outprefix = argv[argpos+1];
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
			++argpos;
		}

		else if (strcmp(argv[argpos],"-biallelic") == 0) {
			// biallelic SNP filter
			biallelic = 1;
		}

		else if (strcmp(argv[argpos], "-allsites") == 0) {
			// process all input sites
			allsites = 1;
		}

		else if (strcmp(argv[argpos], "-allsites_vcf") == 0) {
			// output vcf type
			allsites_vcf = 1;
		}

		else if (strcmp(argv[argpos], "-maxcov") == 0) {
			// max total site depth
			maxcov = atoi(argv[argpos+1]);
			if (maxcov <= 0) {
				std::cerr << "-maxcov must be > 0\n";
				return -1;
			}
			++argpos;
		}

		else if (strcmp(argv[argpos], "-mincov") == 0) {
			// min total site depth
			mincov = atoi(argv[argpos+1]);
			if (mincov < 0) {
				std::cerr << "-mincov must be >= 0\n";
				return -1;
			}
			++argpos;
		}

		else if (strcmp(argv[argpos], "-minind_cov") == 0) {
			// min individual depth
			minind_cov = atoi(argv[argpos+1]);
			if (minind_cov < 0) {
				std::cerr << "-minind_cov must be >= 0\n";
				return -1;
			}
			++argpos;
		}

		else if (strcmp(argv[argpos],"-minind") == 0) {
			// min number of individuals with data
			minind = atoi(argv[argpos+1]);
			if (minind < 0) {
				std::cerr << "-minind must be >= 0\n";
				return -1;
			}
			++argpos;
		}

		else if (strcmp(argv[argpos],"-mingeno") == 0) {
			// min number of individuals with called genotype
			mingeno = atoi(argv[argpos+1]);
			if (mingeno < 0) {
				std::cerr << "-mingeno must be >= 0\n";
				return -1;
			}
			++argpos;
		}

		else if (strcmp(argv[argpos], "-rms_mapq") == 0) {
			// min RMS mapping quality
			rms_mapq = atof(argv[argpos+1]);
			if (rms_mapq < 0.0) {
				std::cerr << "-rms_mapq must be >= 0\n";
				return -1;
			}
			++argpos;
		}

		else if (strcmp(argv[argpos],"-mapqRankSum") == 0) {
			// max map quality bias
			mqRankSum = atof(argv[argpos+1]);
			if (mqRankSum < 0.0) {
				std::cerr << "-mapqRankSum must be >= 0\n";
				return -1;
			}
			++argpos;
		}

		else if (strcmp(argv[argpos], "-posbias") == 0) {
			// max allele position bias
			posbias = atof(argv[argpos+1]);
			if (posbias < 0) {
				std::cerr << "-posbias must be >= 0\n";
				return -1;
			}
			++argpos;
		}

		else if (strcmp(argv[argpos], "-strandbias") == 0) {
			// max strand bias
			strandbias = atof(argv[argpos+1]);
			if (strandbias < 0) {
				std::cerr << "-strandbias must be >= 0\n";
				return -1;
			}
			++argpos;
		}

		else if (strcmp(argv[argpos], "-baseqbias") == 0) {
			// max base quality bias
			baseqbias = atof(argv[argpos+1]);
			if (baseqbias < 0) {
				std::cerr << "-baseqbias must be >= 0\n";
				return -1;
			}
			++argpos;
		}

		else if (strcmp(argv[argpos], "-qual") == 0) {
			// min QUAL field value
			qual = atof(argv[argpos+1]);
			if (qual < 0) {
				std::cerr << "-qual must be >= 0\n";
				return -1;
			}
			++argpos;
		}

		else if (strcmp(argv[argpos], "-varq_depth") == 0) {
			// min variant quality divided by depth of non-homo-ref individuals
			varqual_depth = atof(argv[argpos+1]);
			if (varqual_depth < 0) {
				std::cerr << "-varq_depth must be >= 0\n";
				return -1;
			}
			++argpos;
		}

		else if (strcmp(argv[argpos],"-hetexcess") == 0) {
			// max excess heterozygosity
			hetexcess = atof(argv[argpos+1]);
			if (hetexcess < 0) {
				std::cerr << "-hetexcess must be >= 0\n";
				return -1;
			}
			++argpos;
		}

		else if (strcmp(argv[argpos], "-varonly") == 0) {
			// AF criteria for outputting variable sites
			varonly = 1;
		}

		else if (strcmp(argv[argpos], "-maf") == 0) {
			// MAF lower bound for variable site VCFs
			maf = atof(argv[argpos+1]);
			if (maf < 0.0 || maf > 0.5) {
				std::cerr << "-maf out of [0, 0.5] range\n";
				return -1;
			}
			++argpos;
		}

		else if (strcmp(argv[argpos], "-rmvIndels") == 0) {
			// decide whether to discard sites with indels
			rmvIndels = 1;
		}

		else if (strcmp(argv[argpos], "-siteOnly") == 0) {
			// whether to print individual-level information
			siteOnly = 1;
		}

		else if (strcmp(argv[argpos], "-printFilter") == 0) {
			// whether to filter out QC-failed sites or annotate them in the FILTER field
			printFilter = atoi(argv[argpos+1]);
			switch (printFilter) {
				case 0:
					break;
				case 1:
					break;
				case 2:
					break;
				default:
					std::cerr << "Invalid argument to -printFilter: " << printFilter << "\n";
					return -1;
			}
			++argpos;
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
			++argpos;
		}

		else {
			std::cerr << "Unknown argument " << argv[argpos] << "\n";
			return -1;
		}

		++argpos;
	}

	// output QC failed positions
	if (!printFilter) {
		badpos_name = outprefix + "_fail.pos";
		failpos.open(badpos_name.c_str(), std::ios::out);
		if (! failpos) {
			std::cerr << "Unable to open QC-failed position file " << badpos_name << "\n";
			return -1;
		}
	}

	if (!allsites && allsites_vcf) {
		std::cerr << "Warning: -allsites_vcf is 1 but only SNPs will be printed to VCF because -allsites is 0\n";
	}

	printUserArgs(invcf_name, outvcf_name, goodpos_name, badpos_name, biallelic, allsites, allsites_vcf, maxcov,
			mincov, minind_cov, minind, mingeno, rms_mapq, mqRankSum, posbias, strandbias, baseqbias, qual, varqual_depth, hetexcess,
			varonly, maf, rmvIndels, siteOnly, printFilter);

	return rv;
}

void printUserArgs (const char* invcf_name, std::string &outvcf_name, std::string &goodpos_name, std::string &badpos_name,
	int biallelic, int allsites, int allsites_vcf, unsigned int maxcov, unsigned int mincov, unsigned int minind_cov,
	unsigned int minind, unsigned int mingeno, double rms_mapq, double mqRankSum, double posbias, double strandbias, double baseqbias,
	double qual, double varqual_depth, double hetexcess, int varonly, double maf, int rmvIndels, int siteOnly, int printFilter) {

	int w=20;
	std::cerr << "\n"
	<< std::setw(w) << std::left << "Reading VCF: " << invcf_name << "\n"
	<< std::setw(w) << std::left << "Filtered VCF: " << outvcf_name << "\n"
	<< std::setw(w) << std::left << "Passed positions: " << goodpos_name << "\n";
	if (!printFilter) {
		std::cerr << std::setw(w) << std::left << "Failed positions: " << badpos_name << "\n";
	}
	if (allsites) {
		std::cerr << "-allsites\n";
	}
	if (allsites_vcf) {
		std::cerr << "-allsites_vcf\n";
	}
	if (siteOnly) {
		std::cerr << "-siteOnly\n";
	}
	std::cerr << std::setw(w) << std::left << "-printFilter: " << printFilter << "\n";
	if (biallelic) {
		std::cerr << "-biallelic\n";
	}
	if (varonly) {
		std::cerr << "-varonly\n";
	}
	if (rmvIndels) {
		std::cerr << "-rmvIndels\n";
	}
	std::cerr << std::setw(w) << std::left << "-maxcov: " << maxcov << "\n"
	<< std::setw(w) << std::left << "-mincov: " << mincov << "\n"
	<< std::setw(w) << std::left << "-minind_cov: " << minind_cov << "\n"
	<< std::setw(w) << std::left << "-minind: " << minind << "\n"
	<< std::setw(w) << std::left << "-mingeno: " << mingeno << "\n"
	<< std::setw(w) << std::left << "-rms_mapq: " << rms_mapq << "\n"
	<< std::setw(w) << std::left << "-mapqRankSum: " << mqRankSum << "\n"
	<< std::setw(w) << std::left << "-posbias: " << posbias << "\n"
	<< std::setw(w) << std::left << "-strandbias: "  << strandbias << "\n"
	<< std::setw(w) << std::left << "-baseqbias: " << baseqbias << "\n"
	<< std::setw(w) << std::left << "-qual: " << qual << "\n"
	<< std::setw(w) << std::left << "-varq_depth: " << varqual_depth << "\n"
	<< std::setw(w) << std::left << "-hetexcess: " << hetexcess << "\n"
	<< std::setw(w) << std::left << "-maf: " << maf << "\n"
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

vcfrecord::vcfrecord ()
	: contig(""),
	  pos(0),
	  ref(""),
	  alt(""),
	  af(0),
	  qcfail(0)
{
	flags.resize(16);
	flags[0].first = "Missing";
	flags[1].first = "UknownRef";
	flags[2].first = "MultiAllele";
	flags[3].first = "Indel";
	flags[4].first = "LowQual";
	flags[5].first = "LowDP";
	flags[6].first = "HighDP";
	flags[7].first = "LowIndDP";
	flags[8].first = "LowGeno";
	flags[9].first = "LowMQ";
	flags[10].first = "MapQualBias";
	flags[11].first = "PosBias";
	flags[12].first = "FS";
	flags[13].first = "BaseQualBias";
	flags[14].first = "LowQD";
	flags[15].first = "ExcessHet";

}

void vcfrecord::reserveFields (unsigned int nfields) {
	fields.reserve(nfields);
}

void vcfrecord::reserveFlags (unsigned int nflags) {
	flags.reserve(nflags);
}

void vcfrecord::clearFlags () {
	static std::vector< std::pair<std::string,int> >::iterator it;
	for (it = flags.begin(); it != flags.end(); ++it) {
		it->second = 0;
	}
	qcfail = 0;
}

void vcfrecord::setSite () {
	contig = fields[0];
	pos = atoi(fields[1].c_str());
	ref = fields[3];
	alt = fields[4];
}


void vcfrecord::setFlag (int i) {
	if (flags[i].second == 0) {
		flags[i].second = 1;
		++qcfail;
	}
}
