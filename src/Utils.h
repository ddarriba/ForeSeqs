/*
 * Utils.h
 *
 *  Created on: Oct 2, 2014
 *      Author: diego
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "pll.h"

#include <vector>
#include <iostream>
#include <cstdio>
#include <map>

#define TEST_SIM 1

#define EPSILON 1e-6		  /** epsilon for comparing floating point values */

#define EX_OK EXIT_SUCCESS	  /** exit correctly */
#define EX_UNIMPLEMENTED 1	  /** exit due to an unimplemented feature */
#define EX_IOERR 2			  /** exit due to an input/output error */
#define EX_MEMORY 3			  /** exit due to a memory allocation error */

namespace seqpred {

/** Data Type definition (nucleic or proteic) */
enum DataType { DT_NUCLEIC, DT_PROTEIC};
/** Mode for selecting per-site rate categories */
enum CatMode {
	CAT_RANDOM,		/** Random per-site category */
	CAT_ESTIMATE,	/** Estimated from other partitions */
	CAT_AVERAGE		/** Average of all categories */
	};
/** Mode for stealing the branch lengths */
enum BLMode {
	BL_AVERAGE,		/** Average among all other partitions */
	BL_DRAW,		/** Draw a scaler from an inferred distribution */
	BL_SCALE		/** Find an average scaler */
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
	 * @brief Gets the data type (NT or AA) according to the number of states
	 *
	 * @param[in] pllPartitions The PLL partition list
	 * @param numberOfPartition The number of partition for checking the data type
	 *
	 * @return the data type (DT_NUCLEIC or DT_PROTEIC)
	 */
	static DataType getDataType(const partitionList * pllPartitions, int numberOfPartition);

	/**
	 * @brief Compare 2 nucleotide states
	 *
	 * @param state0 First state. It should be 'a', 'c', 'g' or 't'
	 * @param state1 Second state. It can be either a nucleotide, an ambiguity or a gap
	 * @param[out] validForComp 1, if state1 is 'a', 'c', 'g' or 't'
	 */
	static double compareNucStates(unsigned char state0, unsigned char state1, bool * validForComp);
};

extern BLMode branchLengthsMode;			/** Mode for stealing the branch lengths */
extern CatMode categoriesMode;				/** Mode for selecting per-site rate categories */
extern unsigned int numberOfRateCategories;	/** Number of gamma rate categories */
extern unsigned int numberOfTaxa;			/* Number of taxa */
extern unsigned int sequenceLength;			/* Number of sites (aa/bp) */
extern char ** taxaNames;					/* Name list of taxa */
extern unsigned int * seqIndexTranslate;		/* Array for translating indexes in PLL instance to PLL alignment */

} /* namespace seqpred */

#endif /* UTILS_H_ */
