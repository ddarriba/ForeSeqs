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

BLMode branchLengthsMode = BL_AVERAGE;
CatMode categoriesMode = CAT_ESTIMATE;
PredMode predictionMode = PRED_ANCSEQ;

/* Number of categories hardcoded to 4 */
unsigned int numberOfRateCategories = 4;
unsigned int numberOfTaxa, sequenceLength;
char ** taxaNames;
unsigned int * seqIndexTranslate;

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

double Utils::compareNucStates(unsigned char state0, unsigned char state1, bool * validForComp) {

	*validForComp = state1=='A'||state1=='C'||state1=='G'||state1=='T';
	if(state0 == state1) {
		return 1.0;
	} else {
		/* TODO: Check partial similarity with ambiguities */
		return 0.0;
	}
}

} /* namespace seqpred */
