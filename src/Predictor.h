/*
 * Predictor.h
 *
 *  Created on: Oct 1, 2014
 *      Author: diego
 */

#ifndef PREDICTOR_H_
#define PREDICTOR_H_

#include "pll.h"
#include <vector>

#define PLL_UNDEFINED_SITE 15

namespace seqpred {

class Predictor {
public:
	Predictor(pllInstance * tree, partitionList * partitions, int partitionNumber);
	virtual ~Predictor();

	void predictMissingSequences();
	std::vector<int> getMissingSequences() {
		return missingSequences;
	}
private:
	boolean subtreeIsMissing(nodeptr node);
	nodeptr findMissingDataAncestor();
	std::vector<int> findMissingSequences();
	void mutateSequence ( char * currentSequence, char * ancestralSequence );
	void evolveNode(nodeptr node, char * ancestralSequence);

	pllInstance * tree;
	unsigned int start, end, length;
	partitionList * partitions;
	std::vector<int> missingSequences;
};

} /* namespace partest */

#endif /* PREDICTOR_H_ */
