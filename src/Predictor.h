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
	Predictor(const Predictor&);
	virtual ~Predictor( void );

	/**
	 * @brief Predict the missing sequences for all taxa
	 */
	void predictMissingSequences( void );

	/**
	 * @brief Predict all sequences one by one, for testing purposes
	 */
	void predictAllSequences( void );

	/**
	 * @brief Get the number of taxa with missing sequences
	 */
	int getNumberOfMissingSequences(void) const {
		return _missingSequences.size();
	}

	Predictor& operator=(const Predictor&);

private:
	/**
	* @brief Check whether all taxa data in a subtree is missing
	*/
	boolean subtreeIsMissing(const nodeptr node) const;

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
	void mutateSequence ( char * currentSequence, const char * ancestralSequence, double branchLength );

	/**
	* @brief Predict the sequences for a whole subtree
	*/
	void evolveNode(const nodeptr node, const char * ancestralSequence);

	/**
	 * @brief Compute the branch length for a node
	 */
	double computeBranchLength(const nodeptr node) const;

	pllInstance * _pllTree;				/** PLL instance */
	partitionList * _pllPartitions;		/** PLL list of partitions */
	pllAlignmentData * _pllAlignment;	/** PLL alignment data */

	int _partitionNumber;			/** Partition for predicting the sequences */
	int _numberOfStates;			/** Number of different states (4 for NT, 20 for AA) */
	unsigned int _start;			/** Starting position of the partition */
	unsigned int _end; 				/** Ending position of the partition */
	unsigned int _partitionLength;	/** Number of sites (length) of the partition */
	short * _catToSite;				/** Assignment of categories to sites */

	std::vector<int> _missingSequences;	/** Vector of taxa with missing sequences */
	Model * _currentModel;				/** Model for computing the P matrix */
};

} /* namespace seqpred */

#endif /* PREDICTOR_H_ */
