/*
 * Model.cpp
 *
 *  Created on: Oct 1, 2014
 *      Author: diego
 */

#include "Model.h"

using namespace std;

namespace seqpred {



Model::Model(partitionList * pllPartitions, int partitionIndex) :
		partitionInfo(pllPartitions->partitionData[partitionIndex]) {
}

Model::~Model() {
	/* nothing to do */
}

} /* namespace seqpred */
