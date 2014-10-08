/*
 * ProteinModel.h
 *
 *  Created on: Oct 8, 2014
 *      Author: diego
 */

#ifndef PROTEINMODEL_H_
#define PROTEINMODEL_H_

#include "Model.h"

#include <vector>

#define NUM_AA 20
#define SQNUM_AA NUM_AA*NUM_AA
#define CUNUM_AA NUM_AA*NUM_AA*NUM_AA

namespace seqpred {

class ProteinModel: public Model {
public:
	ProteinModel(partitionList * pllPartitions, int partitionIndex);
	virtual void setMatrix(double * matrix, double branchLength) const;
	virtual char getState(double * Pmatrix) const;
	virtual ~ProteinModel();
private:

	/**
	 * @brief Setup the Q matrix
	 */
	void SetupGTR();

	double Cijk[CUNUM_AA];
	double Root[NUM_AA];
};

} /* namespace seqpred */

#endif /* PROTEINMODEL_H_ */
