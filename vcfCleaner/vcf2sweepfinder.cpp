/*
 * vcf2sweepfinder.cpp
 * This program converts VCF to native SweepFinder format
 * Assumes biallelic sites
 */

# include <iostream>
# include <iomanip>
# include <fstream>
# include <sstream>
# include <cstdlib>
# include <cstring>
# include <map>
# include <vector>

typedef std::map< std::string, std::vector<unsigned long> > fasta_index;
typedef std::pair< std::string, std::vector<char> > faseq;

void vcf2sweepfinderInfo (int poltype, int minderived, int minminor) {
	int w=12;

	std::cerr << "\n"
	<< "vcf2sweepfinder [vcf file] [options]\n"
	<< "\noptions\n"
	<< std::setw(w) << std::left << "-polarize" << std::setw(w) << std::left << "0|1|2" << "SNP polarization method [" << poltype << "]\n"
	<< std::setw(w) << std::left << "-anc" << std::setw(w) << std::left << "FILE" << "ancestral state file\n"
	<< std::setw(w) << std::left << "-minderived" << std::setw(w) << std::left << "INT" << "minimum number of derived alleles [" << minderived << "]\n"
	<< std::setw(w) << std::left << "-minminor" << std::setw(w) << std::left << "INT" << "minimum number of minor alleles at folded sites [" << minminor << "]\n"
	<< "\npolarization method\n"
	<< "0: " << "Do not polarize (all sites are folded with output x column referring to ALT allele counts\n"
	<< "1: " << "Use allele from the VCF 'AA' INFO field as ancestral\n"
	<< "2: " << "Use alleles supplied in a FASTA format ancestral state file as ancestral\n"
	<< "\nOutput columns\n"
	<< "(1) location: position of the SNP\n"
	<< "(2) x: frequency of the SNP\n"
	<< "(3) n: the number of valid sequences at a SNP\n"
	<< "(4) folded: SNP is unfolded (0) or folded (1)\n"
	<< "\n";
}

int getFastaIndex (const char* fasta, fasta_index &faidx) {
	char fastaidx [strlen(fasta) + 5];
	strcat(strcpy(fastaidx, fasta), ".fai");
	int nfields = 5;


	std::ifstream index_stream(fastaidx, std::ifstream::in);
	if (!index_stream) {
		std::cerr << "Unable to open FASTA index " << fastaidx << "\n";
		return -1;
	}


	std::stringstream faistream;
	std::string failine;
	getline(index_stream, failine);

	while(!failine.empty()) {
		faistream.str(std::string());
		faistream.clear();
		faistream.str(failine);
		std::string seqid;
		faistream >> seqid;
		std::vector<unsigned long> fainfo;
		faidx.insert( std::pair< std::string, std::vector<unsigned long> >(seqid, fainfo) );
		faidx[seqid].resize(nfields);
		int i=0;
		while (i < nfields && faistream >> faidx[seqid][i]) {
			++i;
		}
		getline(index_stream, failine);
	}

	return 0;
}

void setFasta (std::ifstream& fafile, fasta_index &faidx, faseq &seq, std::string seqid) {
	// reserve space to store longest fasta sequence
	if (seq.second.capacity() < 1) {
		unsigned int maxlen = 0;
		for (fasta_index::iterator it = faidx.begin(); it != faidx.end(); ++it) {
			unsigned int seqlen = it->second[0];
			if (seqlen > maxlen) maxlen = seqlen;
		}
		seq.second.resize(maxlen, 'N');
	}

	// extract sequence
	seq.first = seqid;
	std::string faline;
	std::stringstream fastream;
	unsigned long offset = faidx[seqid][1];
	fafile.seekg(offset, std::ios_base::beg);

	unsigned int base = 0;
	unsigned int seqlen = faidx[seqid][0];
	while (base < seqlen) {
		fastream.str(std::string());
		fastream.clear();
		getline(fafile, faline);
		if (!faline.empty()) {
			fastream.str(faline);
			while (fastream >> seq.second[base]) {
				++base;
			}
		} else {
			break;
		}
	}

	//std::cerr << seq.first << " length: " << base << "\n";
	//std::cerr << seq.first << ": " << seq.second[base-4] << seq.second[base-3] << seq.second[base-2] << seq.second[base-1] << "\n";

	if (base != seqlen) {
		std::cerr << seqid << " fasta sequence length (" << base << ") does not match index\n";
		seq.second.clear();
	}
}

char getInfoAnc (std::vector<std::string>& vcfentry) {
	// extracts the ancestal allele from the VCF INFO field

	char anc = 'N';
	std::string ancstr;
	std::string* info = &vcfentry[7];
	static std::string::iterator it;

	for (it = info->begin(); it != info->end(); ++it) {
		if (it+3 < info->end()) {
			if (*it == 'A' && *(it+1) == 'A' && *(it+2) == '=') {
				it += 3;
				while (it != info->end() && *it != '|' && *it != ';') {
					ancstr.push_back(*it);
					++it;
				}
				break;
			}
		}
	}

	if (ancstr.empty()) {
		std::cerr << "No ancestral allele found for " << vcfentry[0] << " " << vcfentry[1] << "\n";
		anc = 'N';
	} else if (ancstr.size() > 1) {
		std::cerr << "Invalid ancestral allele '" << ancstr << "' for " << vcfentry[0] << " " << vcfentry[1] << "\n";
		anc = 'N';
	} else {
		anc = ancstr[0];
	}

	return anc;
}


int parseArgs (int argc, char** argv, std::ifstream &vcf, int &anctype, std::ifstream &ancfile, fasta_index &faidx, int &minderived, int &minminor) {

	if (argc < 2){
		vcf2sweepfinderInfo(anctype, minderived, minminor);
		return 1;
	}

	// parse VCF
	vcf.open(argv[1]);
	if (!vcf) {
		std::cerr << "Unable to open VCF " << argv[1] << "\n";
		return -1;
	}

	// parse options
	int argpos = 2;
	int i;

	while (argpos < argc) {
		if (strcmp(argv[argpos], "-polarize") == 0) {
			anctype = atoi(argv[argpos+1]);

			switch (anctype) {
				case 0 :
					std::cerr << "Assuming all sites are folded\n";
					break;
				case 1 :
					std::cerr << "Taking ancestral state from VCF\n";
					break;
				case 2 :
					std::cerr << "Taking ancestral state from FASTA\n";
					for(i=2; i<argc; ++i){
						if (strcmp(argv[i], "-anc") == 0) {
							ancfile.open(argv[i+1]);
							if (!ancfile) {
								std::cerr << "Unable to open ancestral state FASTA " << argv[i+1] << "\n";
								return -1;
							}
							if (getFastaIndex(argv[i+1], faidx)) {
								return -1;
							}
							break;
						}
					}
					if (!ancfile) {
						std::cerr << "Must provide an ancestral state FASTA with -polarize 2\n";
						return -1;
					}
					break;
				default :
					std::cerr << "Unknown polarization method " << anctype << "\n";
					return -1;
			}
		} else if (strcmp(argv[argpos], "-minderived") == 0) {
			minderived = atoi(argv[argpos+1]);
			if (minderived < 0) {
				std::cerr << "-minderived must be >= 0\n";
				return -1;
			}
		} else if (strcmp(argv[argpos], "-minminor") == 0) {
			minminor = atoi(argv[argpos+1]);
			if (minminor < 0) {
				std::cerr << "-minminor must be >= 0\n";
				return -1;
			}
		}
		argpos += 2;
	}

	return 0;
}

unsigned int countAlleles (std::vector<std::string> vcfvec, char anc, int* fold, size_t* n) {
	// counts number of viable sequences and alternate/ancestral alleles
	// returns the ancestral allele count if able to polarize, otherwise returns the alt allele count

	int nref_total = 0;
	int nalt_total = 0;
	static std::vector<std::string>::iterator it;
	static std::string::iterator strit;
	static std::string geno;

	// assumes genotype is always first in each sample's entry as per VCF standard
	for (it = vcfvec.begin()+9; it != vcfvec.end(); ++it) {
		strit = it->begin();
		geno.clear();
		int validfmt = 1;

		// parse genotype
		while (strit != it->end() && *strit != ':') {
			if (*strit == '.') { // missing genotype
				geno.clear();
				break;
			}
			geno.push_back(*strit);
			++strit;
		}

		// check genotype format and count alleles
		int nref = 0;
		int nalt = 0;
		if (geno.length() == 3 && (geno[1] == '/' || geno[1] == '|')) {
			for (int i=0; i < 3; i += 2) {
				if (geno[i] == '0') {
					++nref;
				} else if (geno[i] == '1') {
					++nalt;
				} else {
					validfmt = 0;
					break;
				}
			}
		} else {
			validfmt = 0;
		}

		// add to running totals
		if (validfmt) {
			nref_total += nref;
			nalt_total += nalt;
		} else {
			std::cerr << "Invalid genotype '" << geno << "' at " << vcfvec[0] << " " << vcfvec[1] << "\n";
			*n = 0;
			return 0;
		}
	}

	// check polarization and calculate counts
	*n = nalt_total + nref_total;
	unsigned int allele_count;
	char ref = toupper(vcfvec[3][0]);
	char alt = toupper(vcfvec[4][0]);
	if (anc == ref) {
		allele_count = nalt_total;
		*fold = 0;
	} else if (anc == alt) {
		allele_count = nref_total;
		*fold = 0;
	} else {
		allele_count = nalt_total;
		*fold = 1;
	}

	return allele_count;
}

unsigned int extractAlleleCounts (std::vector<std::string> vcfvec, char anc, int* fold, size_t* n) {
	// extracts allele counts from the VCF INFO field

	std::string* info = &vcfvec[7];
	static std::string acbuffer;
	static std::string anbuffer;
	static std::string::iterator it;
	static std::stringstream ss;

	acbuffer.clear();
	anbuffer.clear();

	for (it = info->begin(); it != info->end(); ++it) {
		if (it+3 < info->end()) {
			if (*it == 'A' && *(it+1) == 'C' && *(it+2) == '=') {
				it +=3;
				while (it != info->end() && *it != ';') {
					acbuffer.push_back(*it);
					++it;
				}
			}
			if (*it == 'A' && *(it+1) == 'N' && *(it+2) == '=') {
				it += 3;
				while (it != info->end() && *it != ';') {
					anbuffer.push_back(*it);
					++it;
				}
			}
		}
		if (!acbuffer.empty() && !anbuffer.empty()) break;
	}

	unsigned int altcount;
	ss.str(std::string());
	ss.clear();
	ss.str(acbuffer);
	ss >> altcount;

	ss.str(std::string());
	ss.clear();
	ss.str(anbuffer);
	ss >> *n;

	// check polarization and calculate counts
	unsigned int allele_count;
	char ref = toupper(vcfvec[3][0]);
	char alt = toupper(vcfvec[4][0]);
	if (anc == ref) {
		allele_count = altcount;
		*fold = 0;
	} else if (anc == alt) {
		allele_count = *n - altcount;;
		*fold = 0;
	} else {
		allele_count = altcount;
		*fold = 1;
	}

	return allele_count;
}

int isBiallelic (std::vector<std::string> vcfvec) {
	 for (int i=3; i<5; ++i) {
		 if (vcfvec[i].size() == 1) {
			 switch (vcfvec[i][0]) {
			 	 case ('A') :
			 			 break;
			 	 case ('C') :
			 			 break;
			 	 case ('G') :
			 			 break;
			 	 case ('T') :
			 			 break;
			 	 case ('N') :
			 			 break;
			 	 default :
			 		 return 0;
			 }
		 } else {
			 return 0;
		 }
	 }

	return 1;
}

std::string& parseVcfHeader (std::ifstream &vcf, std::string &vcfline, int* flags) {
	getline(vcf, vcfline);
	while (!vcfline.empty() && (vcfline[0] == '#' && vcfline[1] == '#')) {
		if (vcfline.substr(0,13) == "##INFO=<ID=AC") flags[0] = 1;
		if (vcfline.substr(0,13) == "##INFO=<ID=AN") flags[1] = 1;
		getline(vcf, vcfline);
	}
	return vcfline;
}

int processVcf (std::ifstream &vcf, const int anctype, std::ifstream &ancfile, fasta_index &faidx, unsigned int minderived, unsigned int minminor) {

	std::string vcfline;
	std::vector<std::string> vcfvec;
	std::vector<std::string>::iterator iter;

	faseq fastaseq;
	fastaseq.first = "";

	// skip VCF header and set up vector to hold tokens
	int countmeth = 0;
	int headerflag [10] = {}; // 0=AC, 1=AN
	if(!parseVcfHeader(vcf, vcfline, headerflag).empty()){
		if (headerflag[0] && headerflag[1]) countmeth = 1;
	} else {
		std::cerr << "VCF is incomplete\n";
		return -1;
	}

	// parse sample info line
	std::stringstream ss(vcfline);
	std::string tok;
	while (ss >> tok) {
		vcfvec.push_back(tok);
	}
	vcfvec.resize(vcfvec.size());
	unsigned int nind = vcfvec.size()-9;

	std::cerr << "Number individuals in VCF: " << nind << "\n";
	if (countmeth) {
		std::cerr << "Extracting allele frequencies from INFO field\n";
	} else {
		std::cerr << "Counting alleles\n";
	}

	// write header
	std::cout << "position\tx\tn\tfolded\n";

	// process site information
	size_t n;
	int fold;
	char ancallele = 'N';
	unsigned int pos;
	unsigned int altfreq = 0;
	unsigned int minorfreq = 0;

	int posn = 0;
	while (!vcf.eof()) {

		// get a new site and vectorize it
		getline(vcf,vcfline);
		if (vcfline.empty()) break;
		ss.str(std::string());
		ss.clear();
		ss.str(vcfline);
		iter = vcfvec.begin();
		while (ss >> *iter) {
			++iter;
		}
		std::istringstream posstream (vcfvec[1]);
		posstream >> pos;

		// skip non-biallelic sites
		if (! isBiallelic(vcfvec)) {
			continue;
		}

		// get ancestral allele for unfolded
		if (anctype == 1) {
			ancallele = toupper(getInfoAnc(vcfvec));
		}
		else if (anctype == 2) {
			// get new fasta sequence if needed
			if (vcfvec[0] != fastaseq.first) {
				setFasta(ancfile, faidx, fastaseq, vcfvec[0]);
				if (fastaseq.second.size() == 0) return -1;
			}
			if (pos > fastaseq.second.size()) {
				std::cerr << vcfvec[0] << " " << vcfvec[1] << " exceeds FASTA sequence length\n";
				return -1;
			}
			ancallele = toupper(fastaseq.second[pos-1]);
		} else {
			ancallele = 'N';
		}
		//std::cerr << vcfvec[0] << " " << vcfvec[1] << ": " << ancallele << "\n";

		// get allele counts
		altfreq = countmeth ? extractAlleleCounts(vcfvec, ancallele, &fold, &n) : countAlleles(vcfvec, ancallele, &fold, &n);

		// write output
		if (n > 0) {
			if (fold) {
				minorfreq = (float)altfreq/n < 0.5 ? altfreq : n-altfreq;
				if (minorfreq < minminor) continue;
			} else {
				if (altfreq < minderived) continue;
			}
			std::cout << pos << "\t" << altfreq << "\t" << n << "\t" << fold << "\n";
		}

		if (posn >= 4) break;
		++posn;
	}

	return 0;
}

int main (int argc, char** argv) {
	int rv = 0;

	std::ifstream vcf;
	std::ifstream ancfile;
	int anctype = 0;
	fasta_index faidx;
	int minderived = 1;
	int minminor = 1;

	if (!parseArgs(argc, argv, vcf, anctype, ancfile, faidx, minderived, minminor)) {
		rv = processVcf(vcf, anctype, ancfile, faidx, minderived, minminor);
	}

	if (vcf.is_open()) vcf.close();
	if (ancfile.is_open()) ancfile.close();

	return rv;
}
