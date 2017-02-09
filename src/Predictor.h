/*
 * Predictor.h
 *
 *  Created on: Oct 1, 2014
 *      Author: Diego Darriba
 *      E-mail: diego.darriba@h-its.org
 *
 *  This file is part of ForeSeqs.
 *
 *  ForeSeqs is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ForeSeqs is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ForeSeqs.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PREDICTOR_H_
#define PREDICTOR_H_

#include "Alignment.h"
#include "Model.h"
#include "PllDefs.h"
#include "Phylo.h"

#include <map>
#include <vector>

#define PLL_UNDEFINED_SITE 15
#define USE_FIXED_ANCESTRAL 1

#define MAX_LOCAL_SCALER 1000

namespace foreseqs {

typedef struct {
	size_t branchNumber;
	double scaler;
	double weight;
} branchInfo;

class Predictor {
public:
	/**
	 * @brief Construct a new Predictor for a single partition
	 *
	 * @param tree The PLL instance
	 * @param partitions The PLL partitions list
	 * @param phylip The PLL alignment
	 * @param partitionNumber The current partition
	 * @param missingSequences The missing sequences in the partition
	 * @param missingBranches The missing branches in the partition
	 */
	Predictor(pll_utree_t * tree,
		        Phylo & phylo,
						Alignment & alignment,
			      unsigned int partitionNumber);

	/**
	 * @brief Clone constructor
	 */
	Predictor(const Predictor&);

	/**
	 * @brief Destructor
	 */
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
	 *
	 * @return The number of taxa with missing sequences
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
	 *
	 * @return The number of predicted partitions
	 */
	unsigned int getMissingPartsCount(void) const {
		return _missingPartsCount;
	}

	/**
	 * Assign operator
	 */
	Predictor& operator=(const Predictor&);

private:

	/**
	 * @brief Check whether the missing subtree ancestors vector contains the node
	 *
	 * @param[in] node The node under consideration
	 *
	 * @return true, if the node is an ancestor in the missing subtree
	 */
	bool isAncestor(pll_utree_t * node) const;

	/**
	 * @brief Check if the branch is missing in a partition
	 *
	 * @param[in] node The node pointing to the branch
	 * @param[in] partition The partition
	 *
	 * @return true, if the branch has missing data in the partition
	 */
	bool isMissingBranch(const pll_utree_t * node, size_t partition) const;

	/**
	 * @brief Mutates a sequence following the current model starting from the ancestor sequence
	 *
	 * @param[in,out] currentSequence The sequence to mutate
	 * @param[in] ancestralSequence The ancestral
	 * @param[in] branchLength The branch length
	 */
	void mutateSequence ( char * currentSequence, const char * ancestralSequence, double branchLength );

	/**
	 * @brief Mutates the P matrix starting from the ancestor sequence
	 *        Multiplies the current P matrix by the ancestral P matrix
	 * @param[in,out] currentPMatrix The P matrix to mutate, and also the result
	 * @param[in] ancestralPMatrix The P matrix at the ancenstral node
	 * @param branchLength The branch length
	 */
	void mutatePMatrix(double * currentPMatrix,
			double * ancestralPMatrix, double branchLength);

	/**
	 * @brief Steals branch lengths recursively
	 *
	 * @param node The root node of the missing subtree
	 */
	void stealBranchRecursive(pll_utree_t * node);

	/**
	 * @brief Predict the sequences for a whole subtree starting from an ancestral sequence
	 *
	 * @param node The node to evolve
	 * @param ancestralSequence The ancestral
	 */
	void evolveNode(pll_utree_t * node, const char * ancestralSequence);

	/**
	 * @brief Predict the sequences for a whole subtree starting from Marginal Ancestral Probabilities
	 *
	 * @param node The node to evolve
	 * @param ancestralProbabilites The marginal ancestral probabilites
	 * @param ancestralPMatrix The P matrix for the parent node
	 */
	void evolveNode(pll_utree_t * node, const double * ancestralProbabilites, double * ancestralPMatrix = 0);

	/**
	 * @brief Compute the branch length for a node
	 *
	 * @param node The node to compute the branch length
	 *
	 * @return The branch length
	 */
	double computeBranchLength(pll_utree_t * node) const;

	/**
	* @brief Recursive algorithm for computing the weighted sum of branch lengths
	*
	* @param node The node for computing
	* @param depth The current depth in the algorithm
	* @param[out] weight The cummulative weight of subtree at the current node
	*/
	double getSumBranches(pll_utree_t * node, int depth, double * weight) const;

	/**
	 * @brief Get the information about the same branch in partitions where it exists
	 *
	 * @param[in] node The node pointing to the input branch
	 * @param[in] depth The depth of the branch in order to weight the scalers
	 * @param[out] cumWeight The cumulative weight of the existing branches
	 * @param[out] branches The list of the branches
	 */
	void getBranches(pll_utree_t * node, int depth, double * cumWeight, std::vector<branchInfo> & branches) const;

	/**
	 * @brief Draws a branch length scaler from a distribution
	 *
	 * @return The branch length scaler
	 */
	double drawBranchLengthScaler( void ) const;

	/**
	 * @brief Finds all rooting nodes for the current partition
	 */
	void getRootingNodes();

	pll_utree_t * _tree;	            /** PLL tree */
	pll_partition_t * _partition;	  /** PLL partition */
	Alignment & _alignment;             /** Alignment data */

	unsigned int _partitionNumber;  /** Partition for predicting the sequences */
	unsigned int _numberOfStates;   /** Number of different states (4 for NT, 20 for AA) */
	unsigned int _start;            /** Starting position of the partition */
	unsigned int _end;              /** Ending position of the partition */
	unsigned int _partitionLength;  /** Number of sites (length) of the partition */
	short * _catToSite;             /** Assignment of categories to sites */

	std::vector<pll_utree_t *> _missingSubtreesAncestors;                 /** Subtrees ancestors */
	std::vector<unsigned int> _missingSequences;                    /** Vector of taxa with missing sequences */
	const std::vector<pll_utree_t *> _missingBranches;   /** Vector of sorted branches with missing sequences */
	unsigned int _missingPartsCount;                                /** Number of predicted partitions */
	Model * _currentModel;                                          /** Model for computing the P matrix */

	double _seqSimilarity;             /** Sequence similarity to original alignment (testing) **/
	double _branchLengthScaler;        /** Scaler for branch length modes 's' and 'd' */
	std::vector<branchInfo> _scalers;  /** Vector of branch scalers information */
};

} /* namespace foreseqs */

#endif /* PREDICTOR_H_ */
