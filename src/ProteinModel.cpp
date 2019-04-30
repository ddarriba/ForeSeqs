/*
 * ProteinModel.cpp
 *
 *  Created on: Oct 8, 2014
 *      Author: Diego Darriba
 *      E-mail: diego.darriba@udc.es
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

#include "ProteinModel.h"
#include "Utils.h"

#include <cassert>
#include <iostream>
#include <cstring>
#include <cmath>

/* for exit() */
#include <cstdlib>

using namespace std;

namespace foreseqs {

ProteinModel::ProteinModel(partitionList * pllPartitions, size_t partitionIndex) :
		Model(pllPartitions, partitionIndex) {

	assert(_pllPartitionInfo->states == NUM_AA);

	size_t numFreqs = NUM_AA;
	size_t numRates = (NUM_AA - 1) * NUM_AA / 2;
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
	for (size_t i = 0; i < NUM_AA; i++) {
		_eigenValues[i] = -_pllPartitionInfo->EIGN[i] / fracchange;
	}

	for (size_t i = 0; i < NUM_AA; i++) {
		for (size_t j = 0; j < NUM_AA; j++) {
			for (size_t k = 0; k < NUM_AA; k++) {
				_Cijk[i * SQNUM_AA + j * NUM_AA + k] = _pllPartitionInfo->EI[i * NUM_AA + k]
										* _pllPartitionInfo->EV[j * NUM_AA + k];
			}
		}
	}

}

char ProteinModel::getState(const double * pMatrix) const {
	size_t j;
	const double *P = pMatrix;
	double r = Utils::genRand();
	for (j=0; r>(*P) && j<NUM_AA-1; j++) P++;
	return (_charStates[j]);
}

char ProteinModel::getMostProbableState(const double * probArray) const {
	int maxPos = -1;
	double max = 0.0;
	const double *P = probArray;
	for (int j = 0; j < NUM_AA; j++) {
		if (*P > max) {
			max = *P;
			maxPos = j;
		}
		P++;
	}
	return (_charStates[(size_t)pow(2,maxPos)]);
}

void ProteinModel::setMatrix(double * matrix, double branchLength, bool cummulative) {
	constructPMatrix(matrix, branchLength, cummulative, NUM_AA, _eigenValues, _Cijk);
}

} /* namespace foreseqs */
