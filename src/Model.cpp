/*
 * Model.cpp
 *
 *  Created on: Oct 1, 2014
 *      Author: diego
 */

#include "Model.h"
#include "Utils.h"

#include <cassert>

using namespace std;

namespace seqpred {

Model::Model(partitionList * pllPartitions, int partitionIndex) :
		_pllPartitionInfo(pllPartitions->partitionData[partitionIndex]) {
}

Model::~Model() {
	/* nothing to do */
}

double Model::computeFracchange( void ) const {

	int numberOfStates = _pllPartitionInfo->states;
	assert ((unsigned int)numberOfStates == _frequencies.size());
	assert ((unsigned int) (numberOfStates * (numberOfStates-1))/2 == _substRates.size());

	/* convert rates into matrix */
	double r[numberOfStates][numberOfStates];
	int i = 0;
	for (int j = 0; j < (numberOfStates-1); j++)
		for (int k = j + 1; k < numberOfStates; k++)
			r[j][k] = _substRates[i++];
	for (int j = 0; j < numberOfStates; j++) {
		r[j][j] = 0.0;
		for (int k = 0; k < j; k++)
			r[j][k] = r[k][j];
	}

	/* evaluate fracchange */
	double fracchange = 0.0;
	for (int j = 0; j < numberOfStates; j++)
		for (int k = 0; k < numberOfStates; k++)
			fracchange += _frequencies[j] * r[j][k] * _frequencies[k];
	return fracchange;

}

} /* namespace seqpred */
