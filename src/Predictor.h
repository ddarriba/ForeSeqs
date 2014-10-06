/*
 * Predictor.h
 *
 *  Created on: Oct 1, 2014
 *      Author: diego
 */

#ifndef PREDICTOR_H_
#define PREDICTOR_H_

#include "pll.h"
#include "Model.h"

#include <map>
#include <vector>

#define PLL_UNDEFINED_SITE 15

namespace seqpred {

class Predictor {
public:
	/**
	* @brief Construct a new Predictor for a single partition
	*/
	Predictor(pllInstance * tree, partitionList * partitions, int partitionNumber);
	virtual ~Predictor();

	/**
	 * @brief Predict the missing sequences for all taxa
	 */
	void predictMissingSequences();

	/**
	 * @brief Get the number of taxa with missing sequences
	 */
	int getNumberOfMissingSequences(void) const {
		return missingSequences.size();
	}

	/**
	* @brief Get the taxa with missing sequences
	*/
	const std::vector<int> getMissingSequences( void ) const {
		return missingSequences;
	}

	/**
	 * @brief Get the new state according to a cumulative probability vector P
	 */
	const std::map<int, char *> getPredictedSequences( void ) const {
		return predictedSequences;
	}
private:
	/**
	* @brief Check whether all taxa data in a subtree is missing
	*/
	boolean subtreeIsMissing(nodeptr node);

	/**
	* @brief Find the farthest common ancestor with all missing data
	*/
	nodeptr findMissingDataAncestor();

	/**
	* @brief Find all the taxa with missing sequence.
	* This method computes the data retrieved by getMissingSequences
	*/
	std::vector<int> findMissingSequences();

	/**
	* @brief Mutates a sequence following the current model starting from the ancestor sequence
	*/
	void mutateSequence ( char * currentSequence, char * ancestralSequence, double branchLength );

	/**
	* @brief Predict the sequences for a whole subtree
	*/
	void evolveNode(nodeptr node, char * ancestralSequence);

	/**
	* @brief Get the new state according to a cumulative probability vector P
	*/
	char getState(double * P);

	pllInstance * tree;					/** PLL instance */
	partitionList * partitions;			/** PLL list of partitions */
	int partitionNumber;				/** Partition for predicting the sequences */
	unsigned int start;					/** Starting position of the partition */
	unsigned int end; 					/** Ending position of the partition */
	unsigned int length;				/** Number of sites (length) of the partition */
	int numStates;						/** Number of states (DNA=4, AA=20) */
	int numRateCategories;				/** Number of gamma rate categories */
	std::vector<int> missingSequences;	/** Vector of taxa with missing sequences */
	Model curModel;						/** Model for computing the P matrix */

	/** Map containing the final predicted sequences */
	std::map<int, char*> predictedSequences;
};

} /* namespace seqpred */

#endif /* PREDICTOR_H_ */
