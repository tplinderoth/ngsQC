/*
 * vcf2eigenstrat.cpp
 * converts VCF into HAPMIX input
 */

# include <fstream>
# include <iostream>
# include <iomanip>
# include <sstream>
# include <cstdlib>
# include <cstring>
# include <vector>
# include <map>

void maininfo () {
	int w=6;
	std::cerr << "\n"
	<< "parseArgs <makeGeno|makeSnp>\n\n"
	<< "makeGeno: " << std::setw(w) << std::left << "generate eigenstrat format genotype file\n"
	<< "makeSnp: " << std::setw(w) << std::left << "generate eigenstrat format SNP file\n"
	<< "\n";
}

void makeGenoInfo (const int genotype) {
	std::cerr << "\nvcf2eigenstrat makeGeno [options] <vcf file | '-' for VCF from STDIN>\n\n"
	<< "Output: each row is a SNP that contains either two columns per individual corresponding\n"
	<< "to the 1st and 2nd phased haplotype allele respectively if phased or one column per individual\n"
	<< "corresponding to the (0,1,2,9) genotype if unphased.\n"
	<< "\nOptions:\n"
	<< "--genotype [0|1]: If providing phased VCF, --genotype 1 will output unphased genotypes as 0,1,2,9 [" << genotype << "]\n"
	<< "\n";
}

void makeSnpInfo () {
	int w = 20;
	std::cerr << "\nvcf2eigenstrat makeSnp <rate input> <vcf file | '-' for VCF from STDIN>\n"
	<< "\nrate input:\n"
	<< "--ratefile [FILE]: " << std::setw(w) << std::left << "File with columns (1) chromosome, (2) position, (3) genetic position in cM\n"
	<< "--rate [FLOAT]: " << std::setw(w) << std::left << "recombination rate in cM/Mb for generating a flat genetic map\n"
	<< "\nOutput:\n"
	<< "(1) SNP name\n"
	<< "(2) chromosome\n"
	<< "(3) genetic position in Morgans\n"
	<< "(4) physical position in bases\n"
	<< "(5) reference allele\n"
	<< "(6) alternate allele\n"
	<< "\n";
}

int parseArgs (int argc, char** argv, std::istream& vcf, std::ifstream& vcffile, std::ifstream& ratefile, double& rate, int& genostyle) {

	int rv = 0;

	if (argc < 2 || strcmp(argv[1], "help") == 0) {
		maininfo();
		return rv;
	}

	// parse VCF file
	if (argc > 2) {
		if (strcmp(argv[argc-1], "-") == 0) {
			vcf.rdbuf(std::cin.rdbuf());
		} else {
			vcffile.open(argv[argc-1]);
			if (!vcf) {
				std::cerr << "Unable to open VCF file " << argv[argc-1] << "\n";
				rv = -1;
				return rv;
			}
		}
	}

	// parse function-specific arguments
	int i = 2;
	if (strcmp(argv[1],"makeGeno") == 0) {
		if (argc < 3) {
			makeGenoInfo(genostyle);
			rv = 0;
			return rv;
		}
		rv = 1;
		while (i < argc-1) {
			if (strcmp(argv[i],"--genotype") == 0) {
				genostyle = atoi(argv[i+1]);
				switch (genostyle) {
					case 0 :
						break;
					case 1 :
						break;
					default :
						std::cerr << "--genotype must be either 1 for unphased genotype or 0 for phased genotype\n";
						rv = -1;
						break;
				}
				if (rv == -1) break;
			} else {
				std::cerr << "Unknown argument " << argv[i] << "\n";
				rv = -1;
				break;
			}
			i += 2;
		}
	} else if (strcmp(argv[1],"makeSnp") == 0) {
		if (argc < 5) {
			makeSnpInfo();
			rv = 0;
			return rv;
		}
		rv = 2;
		while (i < argc-1) {
			if (strcmp(argv[i], "--ratefile") == 0) {
				ratefile.open(argv[i+1]);
				if (!ratefile) {
					std::cerr << "Unable to open recombination rate file " << argv[i+1] << "\n";
					rv = -1;
					break;
				}
			} else if (strcmp(argv[i], "--rate") == 0) {
				rate = atof(argv[i+1]);
				if (rate <= 0) {
					std::cerr << "Recombination rate must be positive\n";
					rv = -1;
					break;
				}
			} else {
				std::cerr << "Unknown argument: " << argv[i] << "\n";
				rv = -1;
				break;
			}
			i += 2;
		}
	} else {
		std::cerr << "Unknown function " << argv[1] << "\n";
		rv = -2;
	}

	return rv;
}

int vcf2geno (std::istream& vcf, const int genotype) {
	int rv = 0;
	std::string vcfline;
	std::vector<std::string> vcfvec;
	std::vector<std::string>::iterator iter;

	// skip VCF header and set up vector to hold tokens
	getline(vcf, vcfline);
	while (!vcfline.empty()) {
		if (vcfline[0] != '#' || vcfline[1] != '#')
			break;
		getline(vcf, vcfline);
	}

	std::stringstream ss(vcfline);
	std::string tok;
	while (ss >> tok) {
		vcfvec.push_back(tok);
	}
	vcfvec.resize(vcfvec.size());
	unsigned int nfields = vcfvec.size();
	unsigned int nind = nfields-9;

	std::cerr << "Number individuals in VCF: " << nind << "\n";

	// process VCF sites
	unsigned int nfields_site = 0;

	while (!vcf.eof()) {

		// get a new VCF line
		getline(vcf,vcfline);
		if (vcfline.empty()) break;
		ss.str(std::string());
		ss.clear();
		ss.str(vcfline);
		iter = vcfvec.begin();
		nfields_site = 0;

		// print genotypes
		while (ss >> *iter) {
			if (nfields_site > 8) {
				if ((*iter)[1] == '|' && !genotype) { // phased haplotypes
					if (((*iter)[0] == '0' || (*iter)[0] == '1') && ((*iter)[2] == '0' || (*iter)[2] == '1')) {
						std::cout << (*iter)[0] << (*iter)[2];
					} else {
						std::cerr << "Invalid genotype " << iter->substr(0,3) << " at " << vcfvec[0] << " " << vcfvec[1] << "\n";
						rv = -1;
						break;
					}
				} else if ((*iter)[1] == '/' || ((*iter)[1] == '|' && genotype)) { // unphased genotypes
					if ((*iter)[0] == '.') {
						std::cout << "9";
					} else if (((*iter)[0] == '0' || (*iter)[0] == '1') && ((*iter)[2] == '0' || (*iter)[2] == '1')) {
						std::cout << atoi(iter->substr(0,1).c_str()) + atoi(iter->substr(2,1).c_str());
					} else {
						std::cerr << "Invalid genotype " << iter->substr(0,3) << " at " << vcfvec[0] << " " << vcfvec[1] << "\n";
						rv = -1;
						break;
					}
				} else {
						std::cerr << "Invalid genotype " << iter->substr(0,3) << " at " << vcfvec[0] << " " << vcfvec[1] << "\n";
						rv = -1;
						break;
					}
				}
			++iter;
			++nfields_site;
		}
		std::cout << "\n";

		if (nfields != nfields_site) {
			std::cout << vcfvec[0] << " " << vcfvec[1] << " VCF line is incomplete\n";
			rv = -1;
		}

		if (rv) break;
	}

	return rv;
}

int vcf2Snp (std::istream& vcf, std::ifstream& ratefile, double rate) {
	int rv = 0;
	std::map<std::string, double> ratemap;
	std::stringstream ss;

	// create map of recombination rates if necessary
	if (ratefile) {
		std::string rateline;
		std::string pos[2];
		double rr;
		getline(ratefile, rateline); // skip header
		while(getline(ratefile, rateline)) {
			ss.clear();
			ss.str(rateline);
			for (int i=0; i<2; ++i) ss >> pos[i];
			ss >> rr;
			ratemap.insert( std::pair<std::string, double>(pos[0] + '_' + pos[1], rr) );
		}
	}

	// add recombination rates to VCF entries
	const int nfields = 5;
	std::string vcfline;
	std::string vcftok[nfields];
	std::vector<std::string> vcfvec;

	// skip VCF header and set up vector to hold tokens
	getline(vcf, vcfline);
	while (!vcfline.empty()) {
		if (vcfline[0] != '#' || vcfline[1] != '#')
			break;
		getline(vcf, vcfline);
	}

	// process VCF sites
	std::map<std::string, double>::iterator it;
	double morgans;
	while (!vcf.eof()) {

		// get a new VCF line
		getline(vcf,vcfline);
		if (vcfline.empty()) break;
		ss.clear();
		ss.str(vcfline);
		int i = 0;

		// extract and print info
		while (i < nfields && ss >> vcftok[i]) ++i;
		std::string snp = vcftok[0] + '_' + vcftok[1];
		std::string chr = vcftok[0].substr(3, vcftok[0].size());
		if  (!ratemap.empty()) {
			if ((it = ratemap.find(snp)) != ratemap.end()) {
				morgans = it->second/100.0; // convert centimorgans to morgans
			} else {
				std::cerr << "WARNING: No rate for VCF position " << vcftok[0] << " " << vcftok[1] << "\n";
			}
		} else {
			morgans = rate * atof(vcftok[1].c_str())/100000000.0; //
		}
		std::cout << snp << "\t" << chr << "\t" << std::setprecision(16) << morgans << "\t" << vcftok[1] << "\t" << vcftok[3] << "\t" << vcftok[4] << "\n";
		if (rv) break;
	}

	return rv;
}

int main (int argc, char** argv) {
	int rv = 0;
	std::ifstream vcffile;
	std::istream vcfstream(vcffile.rdbuf());
	std::ifstream ratefile;
	double rate = 0;
	int genostyle = 0;

	if (!(rv = parseArgs(argc, argv, vcfstream, vcffile, ratefile, rate, genostyle))) {
		return 0;
	} else if (rv == -1) {
		rv = -1;
	} else if (rv == 1) {
		rv = vcf2geno(vcfstream, genostyle);
	} else if (rv == 2) {
		rv = vcf2Snp(vcfstream, ratefile, rate);
	} else {
		rv = -1;
	}

	if (vcffile) vcffile.close();
	if (ratefile) ratefile.close();

	return rv;
}
