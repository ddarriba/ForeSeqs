/*
 * Model.cpp
 *
 *  Created on: Oct 1, 2014
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

#include "Model.h"
#include "Utils.h"

#include <cmath>
#include <cassert>
#include <alloca.h>

using namespace std;

namespace foreseqs {

Model::Model(partitionList * pllPartitions, size_t partitionIndex) :
		_pllPartitionInfo(pllPartitions->partitionData[partitionIndex]),
		_frequencies(), _substRates(), _charStates(), _statesToIntMap() {
}

Model::Model(const Model& other) :
		_pllPartitionInfo(other._pllPartitionInfo),
		_frequencies(other._frequencies), _substRates(other._substRates),
		_charStates(other._charStates), _statesToIntMap(other._statesToIntMap) {

}

Model::~Model() {
	/* nothing to do */
}

Model& Model::operator=(const Model& other) {
	_pllPartitionInfo = other._pllPartitionInfo;
	_frequencies = other._frequencies;
	_substRates = other._substRates;
	_charStates = other._charStates;
	_statesToIntMap = other._statesToIntMap;
	return *this;
}

double Model::computeFracchange( void ) const {

	unsigned int numberOfStates = (unsigned int)_pllPartitionInfo->states;
	assert (numberOfStates == _frequencies.size());
	assert ((numberOfStates * (numberOfStates-1))/2 == _substRates.size());

	/* convert rates into matrix */
	double * r;
	r = (double *) alloca(numberOfStates * numberOfStates * sizeof(double));
	unsigned int i = 0;
	for (unsigned int j = 0; j < (numberOfStates-1); j++)
		for (unsigned int k = j + 1; k < numberOfStates; k++)
			r[j*numberOfStates+k] = _substRates[i++];
	for (unsigned int j = 0; j < numberOfStates; j++) {
		r[j*numberOfStates+j] = 0.0;
		for (unsigned int k = 0; k < j; k++)
			r[j*numberOfStates+k] = r[k*numberOfStates+j];
	}

	/* evaluate fracchange */
	double fracchange = 0.0;
	for (unsigned int j = 0; j < numberOfStates; j++)
		for (unsigned int k = 0; k < numberOfStates; k++)
			fracchange += _frequencies[j] * r[j*numberOfStates+k] * _frequencies[k];
	return fracchange;

}

void Model::constructPMatrix(double * matrix, double branchLength, bool cummulative, size_t numStates, const double * eigenValues, const double * Cijk) {
	double expt[numStates];
	double *P;
	size_t squaredNumStates = numStates * numStates;

	P = matrix;

	if (branchLength < MIN_BRANCH_LEN) {
		/* account for near zero branch lengths */
		for (size_t i = 0; i < numStates; i++) {
			for (size_t j = 0; j < numStates; j++) {
				if (cummulative)
					*P = (i<=j)?1.0:0.0;
				else
					*P = (i==j)?1.0:0.0;
				P++;
			}
		}
		return;
	}
	else
	{
		for (size_t k = 1; k < numStates; k++) {
			expt[k] = exp(branchLength * eigenValues[k]);
		}
		for (size_t i = 0; i < numStates; i++) {
			for (size_t j = 0; j < numStates; j++) {
				(*P) = Cijk[i * squaredNumStates + j * numStates + 0];
				for (size_t k = 1; k < numStates; k++) {
					(*P) += Cijk[i * squaredNumStates + j * numStates + k]
							* expt[k];
				}
				P++;
			}
		}
	}
	if (cummulative) {
		/* the rows are cumulative to help with picking one using
		 a random number */
		for (size_t i = 0; i < numStates; i++) {
			for (size_t j = 1; j < numStates; j++) {
				size_t nextIndex = numStates * i + j;
				matrix[nextIndex] += matrix[nextIndex - 1];
			}
			assert(Utils::floatEquals(matrix[numStates * (i + 1) - 1], 1.0));
		}
	} else {
		/* the matrix rows sum to 1.0 */
		for (size_t i = 0; i < numStates; i++) {
			double sum = 0.0;
			for (size_t j = 0; j < numStates; j++) {
				sum += matrix[numStates * i + j];
			}
			assert(Utils::floatEquals(sum, 1.0));
		}
	}
}

} /* namespace foreseqs */
