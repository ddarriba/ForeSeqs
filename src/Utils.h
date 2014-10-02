/*
 * Utils.h
 *
 *  Created on: Oct 2, 2014
 *      Author: diego
 */

#ifndef UTILS_H_
#define UTILS_H_

namespace seqpred {

class Utils {
public:
	static double genRand(void);
	static void printSequence(const char * seq);
	static void printVector(const double * vec, int len);
	static bool floatEquals(double v1, double v2);
	static int abyx(double a, double x[], int n);
	static int xtoy(double x[], double y[], int n);
	static int matinv(double x[], int n, int m, double space[]);
	static int eigen(int job, double A[], int n, double rr[], double ri[],
			double vr[], double vi[], double work[]);
};

} /* namespace seqpred */

#endif /* UTILS_H_ */
