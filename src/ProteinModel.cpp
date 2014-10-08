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

	assert(_pllPartitionInfo->states == NUM_AA);
	assert(numberOfStates == NUM_AA);

	int numFreqs = NUM_AA;
	int numRates = (NUM_AA - 1) * NUM_AA / 2;
	_frequencies.resize(numFreqs);
	_substRates.resize(numRates);
	memcpy(&(_frequencies[0]), _pllPartitionInfo->frequencies,
			numFreqs * sizeof(double));
	memcpy(&(_substRates[0]), _pllPartitionInfo->substRates,
			numRates * sizeof(double));

	char statesChar[23] = { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
			'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', '-' };
	_charStates.assign(&statesChar[0], &statesChar[0] + 23);
	_statesToIntMap['a'] = _statesToIntMap['A'] = 0;
	_statesToIntMap['r'] = _statesToIntMap['R'] = 1;
	_statesToIntMap['n'] = _statesToIntMap['N'] = 2;
	_statesToIntMap['d'] = _statesToIntMap['D'] = 3;
	_statesToIntMap['c'] = _statesToIntMap['C'] = 4;
	_statesToIntMap['q'] = _statesToIntMap['Q'] = 5;
	_statesToIntMap['e'] = _statesToIntMap['E'] = 6;
	_statesToIntMap['g'] = _statesToIntMap['G'] = 7;
	_statesToIntMap['h'] = _statesToIntMap['H'] = 8;
	_statesToIntMap['i'] = _statesToIntMap['I'] = 9;
	_statesToIntMap['l'] = _statesToIntMap['L'] = 10;
	_statesToIntMap['k'] = _statesToIntMap['K'] = 11;
	_statesToIntMap['m'] = _statesToIntMap['M'] = 12;
	_statesToIntMap['f'] = _statesToIntMap['F'] = 13;
	_statesToIntMap['p'] = _statesToIntMap['P'] = 14;
	_statesToIntMap['s'] = _statesToIntMap['S'] = 15;
	_statesToIntMap['t'] = _statesToIntMap['T'] = 16;
	_statesToIntMap['w'] = _statesToIntMap['W'] = 17;
	_statesToIntMap['y'] = _statesToIntMap['Y'] = 18;
	_statesToIntMap['v'] = _statesToIntMap['V'] = 19;
	_statesToIntMap['b'] = _statesToIntMap['B'] = 20;
	_statesToIntMap['z'] = _statesToIntMap['Z'] = 21;
	_statesToIntMap['-'] = _statesToIntMap['?'] = 22;
	SetupGTR();
}

ProteinModel::~ProteinModel() {
	/* do nothing */
}

void ProteinModel::SetupGTR( void ) {

	double fracchange = computeFracchange();
	for (int i = 0; i < NUM_AA; i++) {
		_eigenValues[i] = -_pllPartitionInfo->EIGN[i] / fracchange;
	}

	for (int i = 0; i < NUM_AA; i++) {
		for (int j = 0; j < NUM_AA; j++) {
			for (int k = 0; k < NUM_AA; k++) {
				_Cijk[i * SQNUM_AA + j * NUM_AA + k] = _pllPartitionInfo->EI[i * NUM_AA + k]
										* _pllPartitionInfo->EV[j * NUM_AA + k];
			}
		}
	}

}

char ProteinModel::getState(const double * pMatrix) const {
	int j;
	const double *P = pMatrix;
	double r = Utils::genRand();
	for (j=0; r>(*P) && j<NUM_AA-1; j++) P++;
	return (_charStates[j]);
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
		expt[k] = exp(branchLength * _eigenValues[k]);
	}
	for (i = 0; i < NUM_AA; i++) {
		for (j = 0; j < NUM_AA; j++) {
			(*P) = _Cijk[i * SQNUM_AA + j * NUM_AA + 0];
			for (k = 1; k < NUM_AA; k++) {
				(*P) += _Cijk[i * SQNUM_AA + j * NUM_AA + k] * expt[k];
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
