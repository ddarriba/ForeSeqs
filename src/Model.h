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

namespace seqgen {

class Model {
public:
	Model(partitionList * pllPartitions, int partitionIndex);
	void setMatrix(double * matrix, double branchLength);
	virtual ~Model();
private:
	pInfo * partitionInfo;
	std::vector<double> frequencies;
	std::vector<double> substRates;
};

} /* namespace seqgen */

#endif /* MODEL_H_ */
