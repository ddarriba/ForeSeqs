/*
 * Utils.h
 *
 *  Created on: Oct 2, 2014
 *      Author: diego
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <iostream>
#include <cstdio>
#include <map>

#define EX_OK EXIT_SUCCESS
#define EX_UNIMPLEMENTED 1
#define EX_IOERR 2

namespace seqpred {

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

	static double genRand(void);
	static void printSequence(const char * seq);
	static void printSequence(const unsigned char * seq, int length);
	static void printVector(const double * vec, int len);
	static bool floatEquals(double v1, double v2);
	static int abyx(double a, double x[], int n);
	static int xtoy(double x[], double y[], int n);
	static int matinv(double x[], int n, int m, double space[]);
	static int eigen(int job, double A[], int n, double rr[], double ri[],
			double vr[], double vi[], double work[]);
	static void init( void );
};

enum DataType { DT_NUCLEIC, DT_PROTEIC};
extern DataType dataType;
extern int numberOfStates;
extern std::vector<char> states;			/** Vector of the different states */
extern std::map<char, int> statesMap;		/** Map of the states index according to char */

extern bool initialized;

} /* namespace seqpred */

#endif /* UTILS_H_ */
