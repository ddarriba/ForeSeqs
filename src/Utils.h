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
	static void printSequence(char * seq);
	static void printVector(double * vec, int len);
};

} /* namespace seqpred */

#endif /* UTILS_H_ */
