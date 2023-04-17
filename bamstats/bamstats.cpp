/*
* bamstats.cpp
*
* Extracts the following stats from BAM file:
* Total depth
* Average map quality
* RMS map quality
* Fraction MQ 0 reads
*
* Compile g++ -O3 -std=c++11 -o bamstats bamstats.cpp
*
* bamstats [options to samtools mpileup] <bam file>
*
* The only additional fields that this script works with is map quality
* produced using samtools mpileup -s/--output-MQ
*
*/

#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <math.h>

void info () {
	std::cerr << "\nbamstats [samtools mpileup options] <bam file>\n\n";
}

int fexists(const char* str)
{
	struct stat buffer;
        return (stat(str, &buffer) == 0);
}


int checkDependencies () {
	if(system("which samtools > /dev/null 2>&1")) {
		std::cerr << "Unable to locate samtools dependency\n";
		return -1;
	}
	return 0;
}

int parseArgs (int argc, char** argv, std::string &bamfile, std::string &cmd) {
	if (argc < 2 || (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))) {
		info();
		return 0;
	}

	if (argc > 2) {
		for (int i = 1; i < argc-1; i++) {
			cmd += static_cast<std::string>(argv[i]) + " ";
		}
	}

	bamfile = argv[argc-1];
	if (!fexists(bamfile.c_str())) {
		std::cerr << "input bam file " << bamfile << " does not exist\n";
		return -1;
	}

	return 0;
}

unsigned int qualstats(const std::string* qstr, double* stats, unsigned int depth) {
	// mapstats[0] = average quality
	// mapstats[1] = RMS quality
	// mapstats[2] = fraction with MQ zero
	static const int offset = 33;
	int q;
	unsigned int n;
	if (depth == 0 || qstr->empty()) {
		for (int i=0; i<3; i++) { stats[i] = -9;}
		n = 0;
	} else {
		for (int i=0; i<3; i++) { stats[i] = 0;}
		n = qstr->length();
		for (unsigned int i = 0; i < n; i++) {
			q = int((*qstr)[i]) - offset;
			stats[0] += q;
			stats[1] += pow(q,2);
			if (q == 0) stats[2]++;
		}
		stats[0] /= n;
		stats[1] = sqrt(stats[1]/n);
		stats[2] /= n;
	}

	return n;
}

int pileupStats (std::string bamfile, std::string cmdoptions) {
	int rv = 0;

	std::string cmd = cmdoptions.empty() ? "samtools mpileup " + bamfile : "samtools mpileup " + cmdoptions + bamfile;
	int coln = (cmd.find("-s") != std::string::npos || cmd.find("--output-MQ") != std::string::npos) ? 7 : 6;

	FILE *fp;
	if ((fp = popen(cmd.c_str(), "r")) == NULL) {
		std::cerr << "Failure reading from pipe: " << cmd << "\n";
		return -1;
	}

	unsigned int buffsize = 1024;
	char buf [buffsize];
	std::string pileline;
	std::string tok;
	double baseq_stats [3];
	double mapq_stats [3];
	unsigned int depth = 0;
	unsigned int nsites = 0;

	//print header
	std::cout << "chr\tpos\tdepth\taverage_baseq\trms_baseq\tfraction_baseq0";
	if (coln > 6) {std::cout << "\taverage_mapq\trms_mapq\tfraction_mapq0";}
	std::cout << "\n";

	while (fgets(buf, buffsize, fp) != NULL) {
		pileline += buf;
		if (pileline[pileline.length() - 1] == '\n') {
			std::stringstream ss(pileline);
			for (int i = 0; i < coln; i++) {
				ss >> tok;
				switch (i) {
					case 0:
						std::cout << tok << "\t";
						break;
					case 1:
						std::cout << tok;
						break;
					case 3:
						depth = std::stoul(tok,NULL,10);
						break;
					case 5:
						if (qualstats(&tok, baseq_stats, depth) != depth) {
							std::cerr << "Base quality string length does not match depth:\n" << pileline;
							rv = -1;
						}
						break;
					case 6:
						if (qualstats(&tok, mapq_stats, depth) != depth) {
							std::cerr << "Map quality string length does not match depth:\n" << pileline;
							rv = -1;
						}
						break;
				}
			}
			if (rv) break;
			std::cout << "\t" << depth;
			for (int j = 0; j < 3; j++) {
				if (baseq_stats[j] == -9) {
					std::cout << "\t.";
				} else {
					std::cout << "\t" << baseq_stats[j];
				}
			}
			if (coln > 6) {
				for (int j = 0; j < 3; j++) {
					if (mapq_stats[j] == -9) {
						std::cout << "\t.";
					} else {
 						std::cout << "\t" << mapq_stats[j];
					}
 				}
			}
			std::cout << "\n";
			pileline.clear();
			++nsites;
		}
	}

	if (!rv && nsites < 1) { std::cerr << "WARNING: '" << cmd << "' returned 0 sites\n"; }
	if (pclose(fp) == -1) {rv = -1;}
	return rv;
}

int main (int argc, char** argv) {
	int rv = 0;
	std::string bamfile;
	std::string cmd = "";

	if (!(rv = parseArgs(argc, argv, bamfile, cmd))) {
		if (!bamfile.empty()) {
			if (!(rv = checkDependencies())) {
				rv = pileupStats(bamfile, cmd);
			}
		}
	}

	if (rv) { std::cerr << "--> exiting\n"; }

	return rv;
}
