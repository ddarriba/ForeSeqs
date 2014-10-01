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

	pllInstance * tree;					/** PLL instance */
	partitionList * partitions;			/** PLL list of partitions */
	unsigned int start;					/** Starting position of the partition */
	unsigned int end; 					/** Ending position of the partition */
	unsigned int length;				/** Number of sites (length) of the partition */
	int numStates;						/** Number of states (DNA=4, AA=20) */
	int numRateCategories;				/** Number of gamma rate categories */
	std::vector<int> missingSequences;	/** Vector of taxa with missing sequences */
};

} /* namespace seqpred */

#endif /* PREDICTOR_H_ */
