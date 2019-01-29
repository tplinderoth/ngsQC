/*
 * expandrf.cpp
 */

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>

void expandinfo () {
	std::cerr << "\n"
	<< "./expandrf [format: bed|pos] [input region file]\n"
	<< "\nformats:\n"
	<< "bed - chr region_start region_end (0-indexed)\n"
	<< "pos - chr position (1-indexed)\n"
	<< "\n";
}

int parseArgs (int argc, char** argv, std::ifstream &rf, std::string& format) {

	int rv = 0;

	if (argc < 3) {
		expandinfo();
		return rv;
	}

	format = argv[1];
	if (format == "bed" || format == "pos") {
	} else {
		std::cerr << "Unknown format " << format << "\n";
		return -1;
	}

	rf.open(argv[2]);
	if (!rf) {
		std::cerr << "Unable to open input region file " << argv[2] << "\n";
		return -1;
	}

	return rv;
}

int doBed (std::ifstream &rf) {
	int rv = 0;
	std::string region;
	std::string chr;
	unsigned int start = 0;
	unsigned int end = 0;
	size_t colon = 0;
	size_t dash = 0;

	getline(rf, region);
	while (!region.empty()) {
		colon = region.find(':');
		dash = region.find('-');
		std::string chr = region.substr(0,colon);
		if (dash != std::string::npos) {
		        start = atoi((region.substr(colon+1, dash)).c_str());
		        end = atoi((region.substr(dash+1, region.size()-dash)).c_str());
		} else {
		        start = atoi((region.substr(colon+1, region.size()-colon)).c_str());
		        end = start;
		}
		--start;

		std::cout << chr << "\t" << start << "\t" << end << "\n";

		getline(rf, region);
	}

	return rv;
}

int doPos (std::ifstream &rf) {
	int rv = 0;
	std::string region;
	std::string chr;
	unsigned int start = 0;
	unsigned int end = 0;
	size_t colon = 0;
	size_t dash = 0;

	getline(rf, region);
	while (!region.empty()) {
		colon = region.find(':');
		dash = region.find('-');
		std::string chr = region.substr(0,colon);
		if (dash != std::string::npos) {
		        start = atoi((region.substr(colon+1, dash)).c_str());
		        end = atoi((region.substr(dash+1, region.size()-dash)).c_str());
		} else {
		        start = atoi((region.substr(colon+1, region.size()-colon)).c_str());
		        end = start;
		}

		for (unsigned int i=start; i <= end; ++i) {
			std::cout << chr << "\t" << i << "\n";
		}

		getline(rf, region);
	}

	return rv;
}

int main (int argc, char** argv) {
	int rv = 0;

	std::string format;
	std::ifstream rf;

	if (!parseArgs(argc, argv, rf, format)) {
		if (format == "bed") {
			rv = doBed(rf);
		} else if (format == "pos") {
			rv = doPos(rf);
		}
	}

	if (rf.is_open()) {
		rf.close();
	}

	return rv;
}
