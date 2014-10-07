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
	Predictor(pllInstance * tree, partitionList * partitions, pllAlignmentData * phylip, int partitionNumber);
	virtual ~Predictor( void );

	/**
	 * @brief Predict the missing sequences for all taxa
	 */
	void predictMissingSequences( void );

	/**
	 * @brief Get the number of taxa with missing sequences
	 */
	int getNumberOfMissingSequences(void) const {
		return missingSequences.size();
	}

private:
	/**
	* @brief Check whether all taxa data in a subtree is missing
	*/
	boolean subtreeIsMissing(nodeptr node) const;

	/**
	* @brief Find the farthest common ancestor with all missing data
	*/
	nodeptr findMissingDataAncestor( void ) const;

	/**
	* @brief Find all the taxa with missing sequence.
	* This method computes the data retrieved by getMissingSequences
	*/
	std::vector<int> findMissingSequences( void ) const;

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
	char getState(double * P) const;

	/**
	 * @brief Compute the branch length for a node
	 */
	double computeBranchLength(nodeptr node) const;

	pllInstance * tree;					/** PLL instance */
	partitionList * partitions;			/** PLL list of partitions */
	pllAlignmentData * phylip;			/** PLL alignment data */
	int partitionNumber;				/** Partition for predicting the sequences */
	unsigned int start;					/** Starting position of the partition */
	unsigned int end; 					/** Ending position of the partition */
	unsigned int partitionLength;				/** Number of sites (length) of the partition */
	int numStates;						/** Number of states (DNA=4, AA=20) */
	int numRateCategories;				/** Number of gamma rate categories */
	std::vector<int> missingSequences;	/** Vector of taxa with missing sequences */
	Model * curModel;						/** Model for computing the P matrix */
};

} /* namespace seqpred */

#endif /* PREDICTOR_H_ */
