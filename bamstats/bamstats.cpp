/*
* bamstats.cpp
*
* Extracts the following stats from BAM file:
* Total depth
* Average mapping and base qualities
* RMS base and mapping qualities
* Fraction base quality zero and mapping quality zero reads
* Number of samples with data
*
* Compile g++ -O3 -o bamstats bamstats.cpp
*
* version 1.1.0
*/

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <sys/stat.h>
#include <math.h>

void info (int qoffset, size_t minind) {
	std::string version("1.1.0");
	int w1 = 32;
	std::cerr << "\nbamstats [options] < BAM | -b BAM_LIST >\n\n"
	<< "The BAM list supplied with -b should contain one BAM file per row.\n"
	<< "\nNative options:\n"
	<< std::setw(w1) << std::left << "--minind INT" << std::setw(w1) << std::left << "Only print sites for which at least INT individuals have data. [" << minind << "]\n"
	<< std::setw(w1) << std::left << "--qual_offset INT" << std::setw(w1) << std::left << "Subtract INT from the raw Phred score values. [" << qoffset << "]\n"
	<< "\nSAMtools mpileup options (set to SAMtools defaults):\n"
	<< std::setw(w1) << std::left << "-A, --count-orphans" << "Keep anamolous read pairs.\n"
	<< std::setw(w1) << std::left << "-a" << "Output all positions, including zero depth sites.\n"
	<< std::setw(w1) << std::left << "-aa" << "Output all positions, including zero depth sites and unused reference sequences.\n"
	<< std::setw(w1) << std::left << "-B, --no-BAQ" << "Disable base alignment quality.\n"
	<< std::setw(w1) << std::left << "-C, --adjust-MQ INT" << std::setw(w1) << std::left << "Adjust mapping quality (0: disable).\n"
	<< std::setw(w1) << std::left << "-A, --count-orphans" << "Keep anamolous read pairs.\n"
	<< std::setw(w1) << std::left << "-d, --max-depth INT" << "Maximum per BAM file depth.\n"
	<< std::setw(w1) << std::left << "-E, --redo-BAQ" << "Recalculate base alignment qualities.\n"
	<< std::setw(w1) << std::left << "-f, --fasta-ref FILE" << std::setw(w1) << std::left  << "Indexed reference sequence in fasta format.\n"
	<< std::setw(w1) << std::left << "-G, --exclude-RG FILE" << std::setw(w1) << std::left << "Exclude read groups listed in FILE.\n"
	<< std::setw(w1) << std::left << "-r, --region STRING" << "Region to analyze.\n"
	<< std::setw(w1) << std::left << "-l, --positions FILE" << std::setw(w1) << std::left << "Skip unlisted positions in BED region or \"chr position\" format.\n"
	<< std::setw(w1) << std::left << "-q, --min-MQ INT" << std::setw(w1) << std::left << "Skip alignments with map quality less than INT.\n"
	<< std::setw(w1) << std::left << "-Q, --min-BQ INT" << std::setw(w1) << std::left << "Skip alignments with base quality less than INT.\n"
	<< std::setw(w1) << std::left << "-A, --count-orphans" << "Keep anamolous read pairs.\n"
	<< std::setw(w1) << std::left << "-R, --ignore-RG" << "Ignore RG tags.\n"
	<< std::setw(w1) << std::left << "--rf, --incl-flags STRING|INT" << std::setw(w1) << std::left  << "Only keep reads with any of the mask bits set (STRING is comma delimited).\n"
	<< std::setw(w1) << std::left << "--ff, --excl-flags STRING|INT" << std::setw(w1) << std::left << "Skip reads with any of the mask bits set (STRING is comma delimited).\n"
	<< std::setw(50) << std::left << "-x, --ignore-overlaps, --disable-overlap-removal" << "Disable reads pair overlap detection and removal.\n"
	<< std::setw(w1) << std::left << "-X, --customized-index FILE" << "Use custome index files.\n"
	<< "\nNotes:\n"
	<< "*Options -a and -aa are incomptable with with --minind values greater than zero.\n"
	<< "*This program expects SAMtools executable to be in the users's PATH.\n"
	<< "\nversion " << version << "\n\n";
}

int fexists(const char* str)
{
	struct stat buffer;
        return (stat(str, &buffer) == 0);
}


int checkDependencies () {
	if(system("which samtools > /dev/null 2>&1")) {
		std::cerr << "Unable to locate SAMtools dependency\n";
		return -1;
	}
	return 0;
}

int argAppend (char** v, int c, int i, std::string &argstr, const char argtype) {
	if (i+1 == c) {
		std::cerr << "Malformed or truncated input argument list\n";
		return -1;
	}
	switch (argtype) {
		case 'f':
			if (!fexists(v[i+1])) {
				std::cerr << v[i] << " file " << v[i+1] << " does not exist\n";
				return -1;
			}
		case 'd':
			if (atof(v[i+1]) < 0) {
				std::cerr << v[i] << " cannot take values less than zero\n";
				return -1;
			}
		case 'i':
			if (atoi(v[i+1]) < 0) {
				std::cerr << v[i] << " cannot take values less than zero\n";
				return -1;
			}
		case 'n':
		default:
			break;
	}
	if (!argstr.empty()) argstr += " ";
	argstr += static_cast<std::string>(v[i]) + " " + static_cast<std::string>(v[i+1]);

	return 0;
}

template <typename T> bool setValue (char** v, int c, int i, T &par) {
	if (i+1 == c) {
		std::cerr << "Malformed or truncated input argument list\n";
		return false;
	}

	std::stringstream ss(v[i+1]);
	if (ss >> par) {
		return true;
	} else {
		std::cerr << "Failed to parse user arguments\n";
		return false;
	}
}

int parseArgs (int argc, char** argv, std::string &bam, size_t &minind, int &qoffset, std::string &sam_arg) {
	if (argc < 2 || (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))) {
		info(qoffset, minind);
		return 0;
	}

	bool allsites = false;

	for (int i = 1; i < argc; i++) {
		if (i == argc-1) {
			bam = argv[i];
			if (!fexists(bam.c_str())) {
				std::cerr << "Input bam file " << bam << " does not exist\n";
				return -1;
			}
			break;
		}
		if (strcmp(argv[i],"-b") == 0 || strcmp(argv[i],"--bam-list") == 0) {
			bam = "";
			if (argAppend(argv, argc, i, bam, 'f') != 0) return -1;
			++i;
		} else if (strcmp(argv[i],"--minind") == 0) {
			if (setValue<size_t>(argv, argc, i, minind) == false) return -1;
			if (minind < 0) {
				std::cerr << "--minind cannot be less than zero\n";
				return -1;
			}
			++i;
		} else if (strcmp(argv[i],"--qual_offset") == 0) {
			if (setValue<int>(argv, argc, i, qoffset) == false) return -1;
			++i;
		} else if (strcmp(argv[i], "-A") == 0 || strcmp(argv[i],"--count-orphans") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 'n') != 0) return -1;
		} else if (strcmp(argv[i],"-B") == 0 || strcmp(argv[i],"--no-BAQ") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 'n') != 0) return -1;
		} else if (strcmp(argv[i],"-C") == 0 || strcmp(argv[i],"--adjust-MQ") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 'd') != 0) return -1;
			++i;
		} else if (strcmp(argv[i],"-d") == 0 || strcmp(argv[i],"--max-depth") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 'i') != 0) return -1;
			++i;
		} else if (strcmp(argv[i],"-E") == 0 || strcmp(argv[i],"--redo-BAQ") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 'n') != 0) return -1;
		} else if (strcmp(argv[i],"-f") == 0 || strcmp(argv[i],"--fasta-ref") == 0) {
			 if (argAppend(argv, argc, i, sam_arg, 'f') != 0) return -1;
			++i;
		} else if (strcmp(argv[i],"-G") == 0 || strcmp(argv[i],"--exclude-RG") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 'f') != 0) return -1;
			++i;
		} else if (strcmp(argv[i],"-l") == 0 || strcmp(argv[i],"--positions") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 'f') != 0) return -1;
			++i;

		} else if (strcmp(argv[i],"-q") == 0 || strcmp(argv[i],"--min-MQ") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 'd') != 0) return -1;
 			++i;
		} else if (strcmp(argv[i],"-Q") == 0 || strcmp(argv[i],"--min-BQ") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 'd') != 0) return -1;
			++i;
		} else if (strcmp(argv[i],"-r") == 0 || strcmp(argv[i],"--region") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 's') != 0) return -1;
			++i;
		} else if (strcmp(argv[i],"-R") == 0 || strcmp(argv[i],"--ignore-RG") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 'n') != 0) return -1;
		} else if (strcmp(argv[i],"--rf") == 0 || strcmp(argv[i],"--incl-flags") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 's') != 0) return -1;
			++i;
		} else if (strcmp(argv[i],"--ff") == 0 || strcmp(argv[i],"--excl-flags") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 's') != 0) return -1;
			++i;
		} else if (strcmp(argv[i],"-x") == 0 || strcmp(argv[i],"--ignore-overlaps-removal") == 0 || strcmp(argv[i],"--disable-overlap-removal") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 'n') != 0) return -1;
		} else if (strcmp(argv[i],"-X") == 0 || strcmp(argv[i],"--customized-index") == 0) {
			 if (argAppend(argv, argc, i, sam_arg, 'n') != 0) return -1;
		} else if (strcmp(argv[i],"-a") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 'n') != 0) return -1;
			allsites = true;
		} else if (strcmp(argv[i],"-aa") == 0) {
			if (argAppend(argv, argc, i, sam_arg, 'n') != 0) return -1;
			allsites = true;
		}
	}

	// check for incompatible arguments/values
	if (allsites && minind > 0) {
		std::cerr << "Option -a or -aa is incompatible with a --minind value greater than zero\n";
		return -1;
	}

	return 0;
}

bool updateQualStats(const std::string* qstr, size_t depth, double* stats, const int offset, std::string* err) {
	int q;
	size_t nreads = qstr->length();
	if (depth == 0) {
		if (nreads == 0 || (*qstr)[0] != '*') {
			*err = "unexpected quality values for BAM with zero depth";
			return false;
		}
	} else if (nreads != depth) {
		*err = "quality string length does not match depth";
		return false;
	} else {
		for (unsigned int i = 0; i < depth; i++) {
			q = int((*qstr)[i]) - offset;
			stats[0] += q; // for average quality
			stats[1] += pow(q,2); // for RMS quality
			if (q == 0) stats[2]++; // for fraction of reads with quality of zero
		}
		stats[3] += depth; // stores total number of quality values
	}
	return true;
}

size_t countBams (const std::string bamarg) {
	std::stringstream ss(bamarg);
	std::string tok;
	size_t nbams = 0;
	ss >> tok;
	if (tok == "-b" || tok == "--bam-list") {
		if (ss.eof()) {
			std::cerr << "BAM list argument, " << tok << ", specified with no BAM list\n";
			return 0;
		}
		ss >> tok;
		std::ifstream bamstream;
		bamstream.open(tok.c_str());
		if (!bamstream) {
			std::cerr << "Unable to open BAM list " << tok << "\n";
			return 0;
		}
		std::string bamfile;
		while (getline(bamstream, bamfile)) {
			if (fexists(bamfile.c_str())) {
				++nbams;
			} else {
				std::cerr << "BAM file, " << bamfile << ", does not exist\n";
				return 0;
			}
		}
		bamstream.close();
	} else {
		if (fexists(tok.c_str())) {
			nbams = 1;
		} else {
			std::cerr << "BAM file, " << tok << ", does not exist\n";
			return 0;
		}
	}
	return nbams;
}

void printStats (double* stats, std::ostream& out) {
	if (stats[3] == 0) {
		out << ".\t.\t.";
	} else {
		for (int i = 0; i < 3; ++i) {
			switch (i) {
				case (0): // Average quality
					out << stats[i]/stats[3];
					break;
				case (1): // RMS quality
					out << sqrt(stats[i]/stats[3]);
					break;
				case (2): // proportion zero quality reads
					out << stats[i]/stats[3];
			}
			if (i < 2) out << "\t";
		}
	}
}

void zeroArray(double (&stats)[4]) {
	for (double& val : stats) {
		val = 0.0;
	}
}

int pileupStats (const std::string bam, std::string sam_options, const int qoffset, const size_t minind) {
	int rv = 0;

	// count number of bams being processed
	size_t nbams = countBams(bam);
	if (!nbams) {
		std::cerr << "Found no valid BAM input\n";
		return -1;
	} else {
		std::cerr << "Number of bams to analyze: " << nbams << "\n";
	}

	std::string cmd = "samtools mpileup " + sam_options + " " + bam;
	std::cerr << "\nGenerating pileup with:\n";
	std::cerr << cmd << "\n";

	FILE *fp;
	if ((fp = popen(cmd.c_str(), "r")) == NULL) {
		std::cerr << "Failure reading from pipe: " << cmd << "\n";
		return -1;
	}

	unsigned int buffsize = 4096; // 1024
	char buf [buffsize];
	std::string pileline;
	std::string tok;
	std::string scaffold;
	size_t pos;
	double baseq_stats [4] = {0};
	double mapq_stats [4] = {0};
	unsigned int depth;
	size_t site_depth;
	unsigned int nsites = 0;
	std::string errmsg;

//	print header
	std::cout << "chr\tpos\tdepth\taverage_baseq\trms_baseq\tfraction_baseq0\taverage_mapq\trms_mapq\tfraction_mapq0\tn_covered\n";

	while (fgets(buf, buffsize, fp) != NULL) {
		pileline += buf;
		if (pileline[pileline.length() - 1] == '\n') {
			std::stringstream ss(pileline);
			site_depth = 0;

			// print position information
			for (int i = 0; i < 3; ++i) {
				switch (i) {
					case (0):
						ss >> scaffold;
						break;
					case (1):
						ss >> pos;
						break;
					default:
						ss >> tok;
				}
			}

			// loop through individuals/bams and extract read information
			size_t ncovered = 0;
			for (size_t k = 1; k <= nbams; ++k) {
				// get info for individual/bam
				for (int i = 0; i < 4; i++) {
					if (ss.eof()) {
						std::cerr << "The following pileup line appears truncated:\n" << pileline << "\n";
						return -1;
					}
					ss >> tok;
					switch (i) {
						case 0:
							depth = std::stoul(tok,NULL,10);
							site_depth += depth;
							if (depth > 0) ++ncovered;
							break;
						case 2:
							if (!updateQualStats(&tok, depth, baseq_stats, qoffset, &errmsg)) {
								std::cerr << "When reading base qualities for BAM " << k << ", " << errmsg << "\n";
								rv = -1;
							}
							break;
						case 3:
							if (!updateQualStats(&tok, depth, mapq_stats, qoffset, &errmsg)) {
								std::cerr << "When reading map qualities for BAM " << k << ", " << errmsg << "\n";
								rv = -1;
							}
							break;
					}
					if (rv) {
						std::cerr << "Offending pileup line:\n" << pileline;
						break;
					}
				}
			}

			// print site information
			if (ncovered >= minind) {
				std::cout << scaffold << "\t" << pos << "\t" << site_depth << "\t";
				printStats(baseq_stats, std::cout);
				std::cout << "\t";
				printStats(mapq_stats, std::cout);
				std::cout << "\t" << ncovered << "\n";
			}
			// reset objects for next pileup line
			pileline.clear();
			zeroArray(baseq_stats);
			zeroArray(mapq_stats);
			++nsites;
		}
	}

	if (!rv && nsites < 1) { std::cerr << "WARNING: '" << cmd << "' returned 0 sites\n"; }
	if (pclose(fp) == -1) {rv = -1;}
	return rv;
}

int main (int argc, char** argv) {
	int rv = 0;
	std::string bam;
	size_t minind = 0;
	int qoffset = 33; // quality score offset
	std::string sam_arg = "-s";

	if (!(rv = parseArgs(argc, argv, bam, minind, qoffset, sam_arg))) {
		if (!bam.empty()) {
			if (!(rv = checkDependencies())) {
				rv = pileupStats(bam, sam_arg, qoffset, minind);
			} else rv = -1;
		} else {
			if (argc > 1)  {
				std::cerr << "No BAM input found.\n";
				rv = -1;
			}
		}
	}

	if (rv) { std::cerr << "--> exiting on error\n"; }

	return rv;
}
