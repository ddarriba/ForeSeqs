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
	virtual void setMatrix(double * pMatrix, double branchLength) const;
	virtual char getState(const double * pMatrix) const;
	virtual ~ProteinModel();
private:

	/**
	 * @brief Setup the Q matrix
	 */
	void SetupGTR( void );

	double _Cijk[CUNUM_AA];
	double _eigenValues[NUM_AA];
};

} /* namespace seqpred */

#endif /* PROTEINMODEL_H_ */
