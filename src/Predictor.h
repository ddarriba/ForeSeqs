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

typedef struct {
	unsigned int branchNumber;
	double scaler;
	double weight;
} branchInfo;

class Predictor {
public:
	/**
	* @brief Construct a new Predictor for a single partition
	*/
	Predictor(pllInstance * tree, partitionList * partitions, pllAlignmentData * phylip,
			unsigned int partitionNumber);
	Predictor(const Predictor&);
	virtual ~Predictor( void );

	/**
	 * @brief Predict the missing sequences for all taxa
	 */
	void predictMissingSequences( const pllAlignmentData * originalSequence = 0 );

	/**
	 * @brief Predict all sequences one by one, for testing purposes
	 */
	void predictAllSequences( void );

	/**
	 * @brief Get the number of taxa with missing sequences
	 */
	unsigned int getNumberOfMissingSequences(void) const {
		return (unsigned int) _missingSequences.size();
	}

	/**
	 * @brief Get the sequence similarity with the original data
	 */
	double getSequenceSimilarity(void) const {
		return _seqSimilarity;
	}

	/**
	 * @brief Get the number of predicted partitions
	 */
	unsigned int getMissingPartsCount(void) const {
		return _missingPartsCount;
	}


	Predictor& operator=(const Predictor&);

private:
	/**
	* @brief Check whether all taxa data in a subtree is missing
	*/
	boolean subtreeIsMissing(const nodeptr node) const;

	/**
	* @brief Find the farthest common ancestor having all missing data
	*
	* @param[in] startingNode Starting node with missing data
	* @return The root of the subtree with all missing data
	*/
	nodeptr findMissingDataAncestor( nodeptr startingNode ) const;

	/**
	* @brief Find all the taxa with missing sequence.
	* This method computes the data retrieved by getMissingSequences
	*/
	std::vector<int> findMissingSequences( void ) const;

	/**
	* @brief Mutates a sequence following the current model starting from the ancestor sequence
	*
	* @param[out] currentSequence The sequence to mutate
	* @param ancestralSequence The ancestral
	* @param branchLength The branch length
	*/
	void mutateSequence ( char * currentSequence, const char * ancestralSequence, double branchLength );

	/**
	* @brief Predict the sequences for a whole subtree
	*
	* @param node The node to evolve
	* @param ancestralSequence The ancestral
	*/
	void evolveNode(const nodeptr node, const char * ancestralSequence);

	/**
	 * @brief Compute the branch length for a node
	 *
	 * @param node The node to compute the branch length
	 * @return The branch length
	 */
	double computeBranchLength(const nodeptr node) const;

	/**
	* @brief Recursive algorithm for computing the weighted sum of branch lengths
	*
	* @param node The node for computing
	* @param depth The current depth in the algorithm
	* @param[out] weight The cummulative weight of subtree at the current node
	*/
	double getSumBranches(nodeptr node, int depth, double * weight) const;

	void getBranches(nodeptr node, int depth, double * cumWeight, std::vector<branchInfo> & branches) const;
	double drawBranchLengthScaler( void ) const;

	pllInstance * _pllTree;				/** PLL instance */
	partitionList * _pllPartitions;		/** PLL list of partitions */
	pllAlignmentData * _pllAlignment;	/** PLL alignment data */

	unsigned int _partitionNumber;	/** Partition for predicting the sequences */
	unsigned int _numberOfStates;	/** Number of different states (4 for NT, 20 for AA) */
	unsigned int _start;			/** Starting position of the partition */
	unsigned int _end; 				/** Ending position of the partition */
	unsigned int _partitionLength;	/** Number of sites (length) of the partition */
	short * _catToSite;				/** Assignment of categories to sites */

	std::vector<nodeptr> _missingSubtreesAncestors; /** Subtrees ancestors */
	std::vector<int> _missingSequences;	/** Vector of taxa with missing sequences */
	unsigned int _missingPartsCount;	/** Number of predicted partitions */
	Model * _currentModel;				/** Model for computing the P matrix */

	double _seqSimilarity;	/** Sequence similarity to original alignment (testing) **/
	double _branchLengthScaler; /** Scaler for branch length modes 's' and 'd' */
	std::vector<branchInfo> _scalers;	/** Vector of branch scalers information */
};

} /* namespace seqpred */

#endif /* PREDICTOR_H_ */
