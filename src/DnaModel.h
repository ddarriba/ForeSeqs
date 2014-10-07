/*
 * DnaModel.h
 *
 *  Created on: Oct 7, 2014
 *      Author: diego
 */

#ifndef DNAMODEL_H_
#define DNAMODEL_H_

#include "Model.h"

#define SQNUM_NUC 4*4
#define NUM_NUC 4
#define CUNUM_NUC 256

namespace seqpred {

class DnaModel: public Model {
public:
	DnaModel(partitionList * pllPartitions, int partitionIndex);
	virtual void setMatrix(double * matrix, double branchLength) const;
	virtual ~DnaModel();
private:
	void SetupGTR();
	double Cijk[CUNUM_NUC];
	double Root[NUM_NUC];
};

} /* namespace seqpred */

#endif /* DNAMODEL_H_ */
