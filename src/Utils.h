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

#define EPSILON 1e-6		/** epsilon for comparing floating point values */

#define EX_OK EXIT_SUCCESS	/** exit correctly */
#define EX_UNIMPLEMENTED 1	/** exit due to an unimplemented feature */
#define EX_IOERR 2			/** exit due to an input/output error */

namespace seqpred {

/** Data Type definition (nucleic or proteic) */
enum DataType { DT_NUCLEIC, DT_PROTEIC};
/** Mode for selecting per-site rate categories */
enum CatMode {
	CAT_RANDOM,		/** Random per-site category */
	CAT_ESTIMATE,	/** Estimated from other partitions */
	CAT_AVERAGE		/** Average of all categories */
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
	 * @param v1,v2 values to compare
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
};

extern CatMode categoriesMode;				/** Mode for selecting per-site rate categories */
extern int numberOfRateCategories;			/** Number of gamma rate categories */
extern int numberOfTaxa;					/* Number of taxa */
extern int sequenceLength;					/* Number of sites (aa/bp) */
extern char ** taxaNames;					/* Name list of taxa */

} /* namespace seqpred */

#endif /* UTILS_H_ */
