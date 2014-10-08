/*
 * ProteinModel.cpp
 *
 *  Created on: Oct 8, 2014
 *      Author: diego
 */

#include "ProteinModel.h"
#include "Utils.h"

#include <cassert>
#include <iostream>
#include <cstring>
#include <cmath>

/* for exit() */
#include <cstdlib>

using namespace std;

namespace seqpred {

ProteinModel::ProteinModel(partitionList * pllPartitions, int partitionIndex) :
		Model(pllPartitions, partitionIndex) {

	assert(partitionInfo->states == NUM_AA);
	assert(numberOfStates == NUM_AA);

	int numFreqs = NUM_AA;
	int numRates = (NUM_AA - 1) * NUM_AA / 2;
	frequencies.resize(numFreqs);
	substRates.resize(numRates);
	memcpy(&(frequencies[0]), partitionInfo->frequencies,
			numFreqs * sizeof(double));
	memcpy(&(substRates[0]), partitionInfo->substRates,
			numRates * sizeof(double));

	char statesChar[23] = { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
			'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', '-' };
	states.assign(&statesChar[0], &statesChar[0] + 23);
	statesMap['a'] = statesMap['A'] = 0;
	statesMap['r'] = statesMap['R'] = 1;
	statesMap['n'] = statesMap['N'] = 2;
	statesMap['d'] = statesMap['D'] = 3;
	statesMap['c'] = statesMap['C'] = 4;
	statesMap['q'] = statesMap['Q'] = 5;
	statesMap['e'] = statesMap['E'] = 6;
	statesMap['g'] = statesMap['G'] = 7;
	statesMap['h'] = statesMap['H'] = 8;
	statesMap['i'] = statesMap['I'] = 9;
	statesMap['l'] = statesMap['L'] = 10;
	statesMap['k'] = statesMap['K'] = 11;
	statesMap['m'] = statesMap['M'] = 12;
	statesMap['f'] = statesMap['F'] = 13;
	statesMap['p'] = statesMap['P'] = 14;
	statesMap['s'] = statesMap['S'] = 15;
	statesMap['t'] = statesMap['T'] = 16;
	statesMap['w'] = statesMap['W'] = 17;
	statesMap['y'] = statesMap['Y'] = 18;
	statesMap['v'] = statesMap['V'] = 19;
	statesMap['b'] = statesMap['B'] = 20;
	statesMap['z'] = statesMap['Z'] = 21;
	statesMap['-'] = statesMap['?'] = 22;
	SetupGTR();
}

ProteinModel::~ProteinModel() {
	// TODO Auto-generated destructor stub
}

void ProteinModel::SetupGTR() {

	double fracchange = computeFracchange();
	for (int i = 0; i < NUM_AA; i++) {
		Root[i] = -partitionInfo->EIGN[i] / fracchange;
	}

	for (int i = 0; i < NUM_AA; i++) {
		for (int j = 0; j < NUM_AA; j++) {
			for (int k = 0; k < NUM_AA; k++) {
				Cijk[i * SQNUM_AA + j * NUM_AA + k] = partitionInfo->EI[i * NUM_AA + k]
										* partitionInfo->EV[j * NUM_AA + k];
			}
		}
	}

}

char ProteinModel::getState(double * P) const {
	int j;
	double r = Utils::genRand();
	for (j=0; r>(*P) && j<NUM_AA-1; j++) P++;
	return (states[j]);
}

void ProteinModel::setMatrix(double * matrix, double branchLength) const {

	int i, j, k;
	double expt[NUM_AA];
	double *P;

	P = matrix;
	if (branchLength < 1e-6) {
		for (i = 0; i < NUM_AA; i++) {
			for (j = 0; j < NUM_AA; j++) {
				if (i == j)
					*P = 1.0;
				else
					*P = 0.0;
				P++;
			}
		}
		return;
	}

	for (k = 1; k < NUM_AA; k++) {
		expt[k] = exp(branchLength * Root[k]);
	}
	for (i = 0; i < NUM_AA; i++) {
		for (j = 0; j < NUM_AA; j++) {
			(*P) = Cijk[i * SQNUM_AA + j * NUM_AA + 0];
			for (k = 1; k < NUM_AA; k++) {
				(*P) += Cijk[i * SQNUM_AA + j * NUM_AA + k] * expt[k];
			}
			P++;
		}
	}

	/* the rows are cumulative to help with picking one using
	 a random number */
	for (int i = 0; i < NUM_AA; i++) {
		for (int j = 1; j < NUM_AA; j++) {
			int nextIndex = NUM_AA * i + j;
			matrix[nextIndex] += matrix[nextIndex - 1];
		}
		assert(Utils::floatEquals(matrix[NUM_AA * (i + 1) - 1], 1.0));
	}

}

} /* namespace seqpred */
