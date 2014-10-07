/*
 * Model.h
 *
 *  Created on: Oct 1, 2014
 *      Author: diego
 */

#ifndef MODEL_H_
#define MODEL_H_

#include <vector>
#include "pll.h"

namespace seqpred {

class Model {
public:
	Model(partitionList * pllPartitions, int partitionIndex);
	virtual void setMatrix(double * matrix, double branchLength) const = 0;
	virtual ~Model();
protected:

	pInfo * partitionInfo;
	std::vector<double> frequencies;
	std::vector<double> substRates;
};

} /* namespace seqpred */

#endif /* MODEL_H_ */
