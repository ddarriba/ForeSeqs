/*
 * Utils.cpp
 *
 *  Created on: Oct 2, 2014
 *      Author: diego
 */

#include "Utils.h"

#include <climits>
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace std;

namespace seqpred {

/* Number of categories hardcoded to 4 */
CatMode categoriesMode = CAT_ESTIMATE;
int numberOfRateCategories = 4;
int numberOfTaxa, sequenceLength;
char ** taxaNames;

double Utils::genRand(void) {
	return rand() * (1.0 / INT_MAX);
}

void Utils::printSequence(const char * sequence) {
	for (unsigned int i = 0; i < strlen(sequence); i++) {
		cout << sequence[i];
	}
	cout << endl;
}

void Utils::printVector(const double * vec, int len) {
	cout << "(";
	for (int i = 0; i < len; i++) {
		cout << vec[i] << ",";
	}
	cout << ")" << endl;
}

bool Utils::floatEquals(double v1, double v2) {
	return (std::abs(v1 - v2) <= EPSILON);
}

DataType Utils::getDataType(const partitionList * pllPartitions, int numberOfPartition) {
	switch (pllPartitions->partitionData[numberOfPartition]->states) {
	case 4:
		return DT_NUCLEIC;
	case 20:
		return DT_PROTEIC;
	default:
		cerr << "[ERROR] Invalid number of states (" <<
		pllPartitions->partitionData[numberOfPartition]->states
		<< ") for partition " << numberOfPartition << endl;
		exit(EX_IOERR);
	}
}

} /* namespace seqpred */
