/*
 * DnaModel.cpp
 *
 *  Created on: Oct 7, 2014
 *      Author: Diego Darriba
 *      E-mail: diego.darriba@h-its.org
 *
 *  This file is part of ForeSeqs.
 *
 *  ForeSeqs is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ForeSeqs is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ForeSeqs.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "DnaModel.h"
#include "Utils.h"

#include <iostream>
#include <cstring>
#include <cassert>
#include <cmath>

using namespace std;

namespace foreseqs {

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

void DnaModel::setMatrix(double * matrix, double branchLength, bool cummulative) {
	constructPMatrix(matrix, branchLength, cummulative, NUM_NUC, _eigenValues, _Cijk);
}

} /* namespace foreseqs */
