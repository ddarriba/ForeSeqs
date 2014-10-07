/*
 * DnaModel.h
 *
 *  Created on: Oct 7, 2014
 *      Author: diego
 */

#ifndef DNAMODEL_H_
#define DNAMODEL_H_

#include "Model.h"

namespace seqpred {

class DnaModel: public Model {
public:
	DnaModel(partitionList * pllPartitions, int partitionIndex);
	virtual void setMatrix(double * matrix, double branchLength) const;
	virtual ~DnaModel();
private:
	void SetupGTR();
};

} /* namespace seqpred */

#endif /* DNAMODEL_H_ */
