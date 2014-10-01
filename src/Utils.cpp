/*
 * Utils.cpp
 *
 *  Created on: Oct 2, 2014
 *      Author: diego
 */

#include "Utils.h"

#include <climits>
#include <cstdlib>
#include <iostream>
#include <cstring>

using namespace std;

namespace seqpred {

double Utils::genRand(void) {
		return rand()*(1.0/INT_MAX);
	}

void Utils::printSequence(char * sequence) {
	for (int i=0; i < strlen(sequence); i++) {
		cout << (int) sequence[i];
	}
	cout << endl;
}

void Utils::printVector(double * vec, int len) {
	for (int i=0; i<len; i++) {
		cout << vec[i] << " ";
	}
	cout << endl;
}

} /* namespace seqpred */
