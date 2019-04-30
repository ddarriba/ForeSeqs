/*
 * Utils.h
 *
 *  Created on: Oct 2, 2014
 *      Author: Diego Darriba
 *      E-mail: diego.darriba@udc.es
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

#ifndef UTILS_H_
#define UTILS_H_

#define PRINT_TRACE 0
#define PRINT_SEQUENCES 0

#include "pll/pll.h"

#include <vector>
#include <iostream>
#include <cstdio>
#include <map>

#define TEST_SIM 0        /** enable testing features */
#define CFD 0		  	  /** correction for distance */

#define EPSILON 1e-6		  /** epsilon for comparing floating point values */

#define MIN_SCALER 1e-4		  /** minimum branch length scaler */
#define MAX_SCALER 1e+4		  /** maximum branch length scaler */

#define EX_OK EXIT_SUCCESS	  /** exit correctly */
#define EX_UNIMPLEMENTED 1	  /** exit due to an unimplemented feature */
#define EX_IOERR 2		  /** exit due to an input/output error */
#define EX_MEMORY 3		  /** exit due to a memory allocation error */

namespace foreseqs {

/** Data Type definition (nucleic or proteic) */
enum DataType { DT_NUCLEIC, DT_PROTEIC};

/** Mode for selecting per-site rate categories */
enum CatMode {
	CAT_RANDOM,	/** Random per-site category */
	CAT_ESTIMATE,	/** Estimated from existing data */
	CAT_AVERAGE,	/** Average of all categories */
	CAT_NONE /* Do not use per-site categories */
	};

/** Mode for stealing the branch lengths */
enum BLMode {
	BL_AVERAGE,	/** Average among partitions with existing data */
	BL_DRAW,	/** Draw a scaler from an inferred distribution */
	BL_SCALE	/** Find an average scaler */
	};

/** Mode for predicting the sequences */
enum PredMode {
	PRED_ANCSEQ,    /** Predict from ancestral sequences */
	PRED_MAP,	    /** Predict from Marginal Ancestral Probabilities */
	PRED_NONE       /** Skip prediction */
};

class Utils {
public:
	/**
	 * @brief Check the existence of a file
	 *
	 * @param filename The file to check
	 *
	 * @return true, if the file exists
	 */
	static inline bool existsFile(const std::string& filename) {
		if (FILE *file = fopen(filename.c_str(), "r")) {
			fclose(file);
			return true;
		} else {
			return false;
		}
	}

	/**
	 * @brief Generate an uniform random value between 0 and 1
	 *
	 * @return random value between 0 and 1
	 */
	static double genRand(void);

	/**
	 * @brief Print a raw sequence
	 *
	 * @param seq the sequence in raw format (e.j., "ACGCCGTA")
	 */
	static void printSequence(const char * seq);

	/**
	 * @brief Print a vector
	 *
	 * @param vec the vector of doubles
	 * @param len the length
	 */
	static void printVector(const double * vec, int len);

	/**
	 * @brief Compare equality of 2 float values
	 *
	 * @param v1 first value to compare
	 * @param v2 second value to compare
	 *
	 * @return true, if both differ less than an epsilon value
	 */
	static bool floatEquals(double v1, double v2);

	/**
	 * @brief Multiply A x B matrices
	 */
	static void matrixMultiply(size_t columns, size_t rows, const double * A, double * B, double * result);

	/**
	 * @brief Gets the data type (NT or AA) according to the number of states
	 *
	 * @param[in] pllPartitions The PLL partition list
	 * @param numberOfPartition The number of partition for checking the data type
	 *
	 * @return the data type (DT_NUCLEIC or DT_PROTEIC)
	 */
	static DataType getDataType(const partitionList * pllPartitions, size_t numberOfPartition);

	/**
	 * @brief Compare 2 nucleotide states
	 *
	 * @param state0 First state. It should be 'a', 'c', 'g' or 't'
	 * @param state1 Second state. It can be either a nucleotide, an ambiguity or a gap
	 * @param[out] validForComp 1, if state1 is 'a', 'c', 'g' or 't'
	 */
	static double compareNucStates(unsigned char state0, unsigned char state1, bool * validForComp);

	/**
	 * @brief Get the list of missing sequences for every partition
	 *
	 * @param[in] pllTree The tree instance
	 * @param[in] pllPartitions The partitions description
	 *
	 * @return The list of missing sequences in the tree
	 */
	static std::vector< std::vector<unsigned int> > findMissingSequences ( pllInstance * pllTree, partitionList * pllPartitions );

	/**
	 * @brief Return the rooting node for a missing branch
	 *
	 * @param[in] pllTree The tree instance
	 * @param[in,out] missingSequences The list of missing sequences. Those in the missing subtree are removed
	 * @param[in] startingNode The missing branch for start searching
	 * @param[out] missingBranches The list of missing branches in the missing subtree
	 *
	 * @return The rooting node
	 */
	static nodeptr findRootingNode( pllInstance * pllTree, std::vector<unsigned int> * missingSequences,
		nodeptr startingNode, std::vector<nodeptr> * missingBranches );

	/**
	 * @brief Get the list of missing branches for every partition
	 *
	 * @param[in] pllTree The tree instance
	 * @param[in] pllPartitions The partitions description
	 *
	 * @return The list of missing branches in the tree
	 */
	static std::vector< std::vector<nodeptr> > findMissingBranches ( pllInstance * pllTree,
		partitionList * pllPartitions, std::vector< std::vector<unsigned int> > missingSequences );

	/**
	 * @brief Optimizes the model parameters
	 *
	 * @param[in,out] pllTree The tree
	 * @param[in,out] pllPartitions The partitions
	 */
	static void optimizeModelParameters(pllInstance * pllTree, partitionList * pllPartitions);

	/**
	 * @brief Allocate memory
	 * @param n number of elements to allocate
	 * @param el_size per-element size
	 *
	 * @return Pointer to the allocated memory. 0, if error
	 */
	static void * allocate(size_t n, size_t el_size);

private:

	/**
	 * @brief Check if the subtree is a missing subtree
	 *
	 * @param[in] pllTree The tree instance
	 * @param[in,out] missingSequences The list of missing sequences. Those in the missing subtree are removed
	 * @param[in] node The missing branch for start searching
	 * @param[out] missingBranches The list of missing branches in the missing subtree, sorted by node number
	 *
	 * @return The rooting node
	 */
	static boolean subtreeIsMissing( pllInstance * pllTree, std::vector<unsigned int> * missingSequences,
		const nodeptr node, std::vector<nodeptr> * missingBranches );

	/**
	 * @brief Transposes a matrix
	 *
	 * @param[in,out] m The matrix
	 * @param size The size of the matrix
	 */
	static void transpose(double * m, size_t size);

};

extern BLMode branchLengthsMode;		/** Mode for stealing the branch lengths */
extern CatMode categoriesMode;			/** Mode for selecting per-site rate categories */
extern PredMode predictionMode;			/** Prior for predicting the sequences */
extern unsigned int numberOfRateCategories;	/** Number of gamma rate categories */
extern unsigned int numberOfTaxa;		/* Number of taxa */
extern unsigned int numberOfThreads;	/* Number of threads */
extern unsigned int numberOfPartitions; /* Number of partitions */
extern double threshold;				/* Threshold for considering sequences as missing data */
extern unsigned int sequenceLength;		/* Number of sites (aa/bp) */
extern char ** taxaNames;			/* Name list of taxa */
extern unsigned int * seqIndexTranslate;	/* Array for translating indexes in PLL instance to PLL alignment */
extern bool predictSequences;           /* Flag for trigger the sequence prediction */

} /* namespace foreseqs */

#endif /* UTILS_H_ */
