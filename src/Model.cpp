/*
 * Model.cpp
 *
 *  Created on: Oct 1, 2014
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

#include "Model.h"
#include "Utils.h"

#include <cassert>
#include <alloca.h>

using namespace std;

namespace seqpred {

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

} /* namespace seqpred */
