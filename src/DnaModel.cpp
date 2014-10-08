/*
 * DnaModel.cpp
 *
 *  Created on: Oct 7, 2014
 *      Author: diego
 */

#include "DnaModel.h"
#include "Utils.h"

#include <iostream>
#include <cstring>
#include <cassert>
#include <cmath>

using namespace std;

namespace seqpred {

DnaModel::DnaModel(partitionList * pllPartitions, int partitionIndex) :
		Model(pllPartitions, partitionIndex) {

	assert (_pllPartitionInfo->states == NUM_NUC);
	assert (numberOfStates == NUM_NUC);

	int numFreqs = NUM_NUC;
	int numRates = (NUM_NUC - 1) * NUM_NUC / 2;
	_frequencies.resize(numFreqs);
	_substRates.resize(numRates);
	memcpy(&(_frequencies[0]), _pllPartitionInfo->frequencies, numFreqs * sizeof(double));
	memcpy(&(_substRates[0]), _pllPartitionInfo->substRates, numRates * sizeof(double));

	char statesChar[16] = {'_', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', '-'};
	_charStates.assign(&statesChar[0], &statesChar[0]+16);
	_statesToIntMap['a'] = _statesToIntMap['A'] = 0;
	_statesToIntMap['c'] = _statesToIntMap['C'] = 1;
	_statesToIntMap['g'] = _statesToIntMap['G'] = 2;
	_statesToIntMap['t'] = _statesToIntMap['T'] = 3;

	SetupGTR();
}

DnaModel::~DnaModel() {
	/* do nothing */
}

void DnaModel::SetupGTR( void ) {

	double fracchange = computeFracchange();
	for (int i = 0; i < NUM_NUC; i++) {
		_eigenValues[i] = -_pllPartitionInfo->EIGN[i] / fracchange;
	}

	for (int i = 0; i < NUM_NUC; i++) {
		for (int j = 0; j < NUM_NUC; j++) {
			for (int k = 0; k < NUM_NUC; k++) {
				_Cijk[i * SQNUM_NUC + j * NUM_NUC + k] = _pllPartitionInfo->EI[i * NUM_NUC + k]
										* _pllPartitionInfo->EV[j * NUM_NUC + k];
			}
		}
	}

}

char DnaModel::getState(const double * pMatrix) const {
	int j;
	const double *P = pMatrix;
	double r = Utils::genRand();
	for (j=0; r>(*P) && j<NUM_NUC-1; j++) P++;
	return (_charStates[pow(2,j)]);
}

void DnaModel::setMatrix(double * matrix, double branchLength) const {
	int i, j, k;
	double expt[NUM_NUC];
	double *P;

	P = matrix;
	if (branchLength < 1e-6) {
		for (i = 0; i < NUM_NUC; i++) {
			for (j = 0; j < NUM_NUC; j++) {
				if (i == j)
					*P = 1.0;
				else
					*P = 0.0;
				P++;
			}
		}
		return;
	}

	for (k = 1; k < NUM_NUC; k++) {
		expt[k] = exp(branchLength * _eigenValues[k]);
	}
	for (i = 0; i < NUM_NUC; i++) {
		for (j = 0; j < NUM_NUC; j++) {
			(*P) = _Cijk[i * SQNUM_NUC + j * NUM_NUC + 0];
			for (k = 1; k < NUM_NUC; k++) {
				(*P) += _Cijk[i * SQNUM_NUC + j * NUM_NUC + k] * expt[k];
			}
			P++;
		}
	}

	/* the rows are cumulative to help with picking one using
	 a random number */
	for (int i = 0; i < NUM_NUC; i++) {
		for (int j = 1; j < NUM_NUC; j++) {
			int nextIndex = NUM_NUC * i + j;
			matrix[nextIndex] += matrix[nextIndex - 1];
		}
		assert(Utils::floatEquals(matrix[NUM_NUC * (i + 1) - 1], 1.0));
	}

}

} /* namespace seqpred */
