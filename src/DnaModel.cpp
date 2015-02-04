/*
 * DnaModel.cpp
 *
 *  Created on: Oct 7, 2014
 *      Author: Diego Darriba
 *      E-mail: diego.darriba@h-its.org
 *
 *  This file is part of SeqPred.
 *
 *  SeqPred is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SeqPred is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with SeqPred.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "DnaModel.h"
#include "Utils.h"

#include <iostream>
#include <cstring>
#include <cassert>
#include <cmath>

using namespace std;

namespace seqpred {

DnaModel::DnaModel(partitionList * pllPartitions, size_t partitionIndex) :
		Model(pllPartitions, partitionIndex) {

	assert (_pllPartitionInfo->states == NUM_NUC);

	size_t numFreqs = NUM_NUC;
	size_t numRates = (NUM_NUC - 1) * NUM_NUC / 2;
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
	for (size_t i = 0; i < NUM_NUC; i++) {
		_eigenValues[i] = -_pllPartitionInfo->EIGN[i] / fracchange;
	}

	for (size_t i = 0; i < NUM_NUC; i++) {
		for (size_t j = 0; j < NUM_NUC; j++) {
			for (size_t k = 0; k < NUM_NUC; k++) {
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
	return (_charStates[(size_t)pow(2,j)]);
}

char DnaModel::getMostProbableState(const double * probArray) const {
	int maxPos = -1;
	double max = 0.0;
	const double *P = probArray;
	for (int j = 0; j < NUM_NUC; j++) {
		if (*P > max) {
			max = *P;
			maxPos = j;
		}
		P++;
	}
	return (_charStates[(size_t)pow(2,maxPos)]);
}

void DnaModel::setMatrix(double * matrix, double branchLength, bool cummulative) const {
	double expt[NUM_NUC];
	double *P;

	P = matrix;
	if (branchLength < 1e-6) {
		for (size_t i = 0; i < NUM_NUC; i++) {
			for (size_t j = 0; j < NUM_NUC; j++) {
				if (i == j)
					*P = 1.0;
				else
					*P = 0.0;
				P++;
			}
		}
		return;
	}

	for (size_t k = 1; k < NUM_NUC; k++) {
		expt[k] = exp(branchLength * _eigenValues[k]);
	}
	for (size_t i = 0; i < NUM_NUC; i++) {
		for (size_t j = 0; j < NUM_NUC; j++) {
			(*P) = _Cijk[i * SQNUM_NUC + j * NUM_NUC + 0];
			for (size_t k = 1; k < NUM_NUC; k++) {
				(*P) += _Cijk[i * SQNUM_NUC + j * NUM_NUC + k] * expt[k];
			}
			P++;
		}
	}

	if (cummulative) {
		/* the rows are cumulative to help with picking one using
		 a random number */
		for (size_t i = 0; i < NUM_NUC; i++) {
			for (size_t j = 1; j < NUM_NUC; j++) {
				size_t nextIndex = NUM_NUC * i + j;
				matrix[nextIndex] += matrix[nextIndex - 1];
			}
			assert(Utils::floatEquals(matrix[NUM_NUC * (i + 1) - 1], 1.0));
		}
	} else {
		/* the matrix rows sum to 1.0 */
		for (size_t i = 0; i < NUM_NUC; i++) {
			double sum = 0.0;
			for (size_t j = 0; j < NUM_NUC; j++) {
				sum += matrix[NUM_NUC * i + j];
			}
			assert(Utils::floatEquals(sum, 1.0));
		}
	}

}

} /* namespace seqpred */
