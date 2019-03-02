/*
 * vcf2bimbam.cpp
 * this program converts VCF to BIMBAM format with posterior expected genotypes for use with GEMMA
 * only works for biallelic sites
 */

# include <fstream>
# include <iostream>
# include <iomanip>
# include <vector>
# include <cstring>
# include <sstream>
# include <cstdlib>
# include <math.h>

void vcf2bimbamInfo () {
	int w=6;
	std::cerr << "\n"
	<< "vcf2bimbam [input vcf] [genotype prior]\n"
	<< "\npriors:\n"
	<< std::setw(w) << std::left << "hwe" << "Hardy-Weinberg Equilibrium probability of the genotypes given the ALT allele frequency\n"
	<< std::setw(w) << std::left << "unif" << "Uniform probability for genotypes\n\n"
	<< "BIMBAM output:\n"
	<< "(1) SNP, (2) ALT, (3) REF, (4..N) Posterior expected genotype\n"
	<< "Missing genotypes encoded as 'NA'\n"
	<< "\n";
}

int parseArgs (int argc, char** argv, std::ifstream &vcf, int &genoprior) {
	int rv = 0;

	if (argc < 3) {
		vcf2bimbamInfo();
		return 1;
	}

	vcf.open(argv[1]);
	if (!vcf) {
		std::cerr << "Unable to open input VCF " << argv[1] << "\n";
		return -1;
	}

	if (strcmp(argv[2], "unif") == 0) {
		genoprior = 0;
	} else if (strcmp(argv[2], "hwe") == 0) {
		genoprior = 1;
	} else {
		std::cerr << "Invalid genotype prior " << argv[2] << "\n";
		return -1;
	}

	return rv;
}

int plIndex (std::vector<std::string> &vcfvec) {
	// determine where the PL info is in the FORMAT field

	int index = 0;
	int plfound = 0;

	for (unsigned int i=0; i<vcfvec[8].size(); ++i) {
		if (vcfvec[8][i] == 'P' && vcfvec[8][i+1] == 'L') {
			plfound = 1;
		}
		if (vcfvec[8][i] == ':') ++index;
		if (plfound) break;
	}

	if (!plfound) {
		index = -1;
		std::cerr << "ERROR: No PL found in FORMAT field for " << vcfvec[0] << " " << vcfvec[1] << "\n";
	}

	return index;
}

double calcGeno (double* gl, double* prior) {
	static double d [3];
	double pdata = 0.0;
	int i;
	for (i=0; i<3; ++i) {
		d[i] = gl[i] * prior[i];
		pdata += d[i];
	}

	if (pdata == 0) {
		std::cerr << "Data has probability zero\n";
		return -9.0; // use better exception handling here
	}

	double postgeno = 0.0;
	for (i=1; i<3; ++i) {
		postgeno += i * (d[i]/pdata);
	}

	//std::cerr << gl[0] << "," << gl[1] << "," << gl[2] << ": " << postgeno << "\n"; // debug

	return postgeno;
}

int expectGeno (std::vector<std::string> &vcfvec, double* geno, unsigned int nind, double* genoprior) {
	int rv = 0;

	// find location of PL in format field
	int plidx;
	if ((plidx = plIndex(vcfvec)) < 0) {
		return -1;
	}

	// get genotype likelihoods
	static double pl [3];
	static char plstr [10];
	static std::vector<std::string>::iterator inditer;
	unsigned int n=0;

	// loop through all individuals
	for (inditer=vcfvec.begin()+9; inditer < vcfvec.end(); ++inditer) {
		if (n >= nind) {
			std::cerr << "Number of individuals at " << vcfvec[0] << " " << vcfvec[1] << " is greater than number of individuals in header\n";
			return -1;
		}

		// loop through individual's information
		int k = 0;
		pl[0] = -9;
		for (unsigned int i=0; i<inditer->size(); ++i) {
			if (k == plidx) {
				// loop through genotype likelihood information
				plstr[0] = '\0';
				int j = 0, m = 0;
				while (i < inditer->size()) {
					if ((*inditer)[i] == ':') {
						break;
					} else if ((*inditer)[i] == ',') {
						plstr[j] = '\0';
						if (m < 3) {
							pl[m] = pow(10, -atof(plstr)/10); // convert out of Phred-scale
							++m;
							j = 0;
						} else {
							std::cerr << "ERROR: More than 3 PL values for " << vcfvec[0] << " " << vcfvec[1] << " individual " << n << "\n";
							return -1;
						}
					} else {
						plstr[j] = (*inditer)[i];
						++j;
					}
					++i;
				}
				if (m < 3) {
					plstr[j] = '\0';
					pl[m] = pow(10, -atof(plstr)/10); // convert out of phred scale
				} else {
					std::cerr << "ERROR: More than 3 PL values for " << vcfvec[0] << " " << vcfvec[1] << " individual " << n << "\n";
					return -1;
				}
			} else if ((*inditer)[i] == ':') {
				++k;
			}
		}

		// calculate mean genotype
		if (pl[0] == -9.0) {
			// missing genotype likelihood information
			geno[n] = -9;
		} else {
			geno[n] = calcGeno(pl, genoprior);
			if (geno[n] == -9.0) {
				std::cerr << "WARNING: Problem calculating genotype for " << vcfvec[0] << " " << vcfvec[1] << " individual " << n << "\n";
			}
		}

		++n;
	}

	return rv;
}

double getaf (const std::vector<std::string> &vcfvec) {
	static char f [10];
	f[0] = '\0';
	int j=0;
	static std::string info;
	info.clear();
	info = (vcfvec[7]);
	for (unsigned int i=0; i<info.size(); ++i) {
		if (i+2 < info.size() && info[i]=='A' && info[i+1]=='F' && info[i+2]=='=') {
			i += 3;
			while (i < info.size() && info[i] != ';') {
				f[j] = info[i];
				++j;
				++i;
			}
			f[j] = '\0';
			break;
		}
	}

	if (!f[0]) {
		std::cerr << "ERROR: No AF found in INFO field\n";
		return -1.0;
	}

	return atof(f);
}

void hweprior (double *prob, double altf) {
	// prob[0] = p(REF/REF), p[1] = p(REF/ALT), p[2] = p(ALT/ALT)
	double reff = 1.0 - altf;
	prob[0] = reff*reff;
	prob[1] = 2*reff*altf;
	prob[2] = altf*altf;
}

void printGenoInfo (std::vector<std::string> &vcfvec, double* geno, unsigned int n) {
	std::string id(vcfvec[0] + "_" + vcfvec[1]);
	std::cout << id << ", " << vcfvec[4] << ", " << vcfvec[3] << ", ";

	unsigned int i=0;
	for (i=0; i<n-1; ++i) {
		if (geno[i] == -9.0) {
			// missing genotype
			std::cout << "NA, ";
		} else
		{
			std::cout << geno[i] << ", ";
		}
	}

	// print last genotype
	if (geno[i] == -9.0) {
		std::cout << "NA\n";
	} else {
		std::cout << geno[i] << "\n";
	}
}

int processVcf (std::ifstream &vcf, const int priortype) {
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
	unsigned int nind = vcfvec.size()-9;

	std::cerr << "Number individuals in VCF: " << nind << "\n";

	// process sites in VCF
	double af;
	double gprior [3] = {1.0/3, 1.0/3, 1.0/3};
	double * geno = new double[nind];

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

		// get alternate allele frequency from INFO field for prior
		if (priortype == 1) {
			if ((af=getaf(vcfvec)) < 0.0) {
				rv = -1;
				break;
			}
			hweprior(gprior, af);
		}

		// calculate expected genotypes
		if (expectGeno(vcfvec, geno, nind, gprior)) {
			rv = -1;
			break;
		}

		// print expected genotype information
		printGenoInfo(vcfvec, geno, nind);
	}

	delete [] geno;
	return rv;
}

int main(int argc, char** argv) {
	int rv = 0;
	std::ifstream vcf;
	int genoprior; // 0=unif, 1=hwe

	if (!parseArgs(argc, argv, vcf, genoprior)) {
		processVcf(vcf, genoprior);
	}

	if (vcf.is_open()) vcf.close();

	return rv;
}
