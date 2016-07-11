/*
 * seqstat.cpp
 *
 *  Created on: Jul 3, 2016
 *      Author: tyler
 */

#include "seqstat.h"
#include <map>


double seqstat::calcDustScore (const std::string& seq) {
/*
 * From Morgulis et al. 2006, Journal of Computational Biology
 * A Fast and Symmetric DUST Implementation to Mask Low-Complexity DNA Sequences
 */
	unsigned int width = 3;
	static std::map<std::string, int> subscores;

	// sequence must be at least 3 nt to calculate score
	if (seq.size() <= width)
		return 0.0;

	// slide width-mer window over sequence and place subsequences into the map
	for (size_t i=0; i<=seq.size()-width; ++i)
		subscores[seq.substr(i, width)]++;

	// calculate score for entire sequence
	double sum = 0;
	int ct;
	std::map<std::string, int>::iterator iter;
	for (iter = subscores.begin(); iter != subscores.end(); ++iter)
	{
		ct = iter->second;
		sum += static_cast<double>(ct*(ct-1))/2.0;
		iter->second = 0;
	}
	int l = seq.size() - (width-1); // number of windows (subsequences) in sequence
	return sum/(l-1);
}


