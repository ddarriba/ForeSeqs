/*
 * DnaModel.h
 *
 *  Created on: Oct 7, 2014
 *      Author: diego
 */

#ifndef DNAMODEL_H_
#define DNAMODEL_H_

#include "Model.h"

#define NUM_NUC 4
#define SQNUM_NUC 4*4
#define CUNUM_NUC 4*4*4

namespace seqpred {

class DnaModel: public Model {
public:
	DnaModel(partitionList * pllPartitions, int partitionIndex);
	virtual void setMatrix(double * pMatrix, double branchLength) const;
	virtual char getState(const double * pMatrix) const;
	virtual ~DnaModel();
private:

	/**
	 * @brief Setup the Q matrix
	 */
	void SetupGTR( void );

	double _Cijk[CUNUM_NUC];
	double _eigenValues[NUM_NUC];
};

} /* namespace seqpred */

#endif /* DNAMODEL_H_ */
