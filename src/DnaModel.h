/*
 * DnaModel.h
 *
 *  Created on: Oct 7, 2014
 *      Author: Diego Darriba
 *      E-mail: diego.darriba@h-its.org
 *
 *  This file is part of ForeSeqs.
 *
 *  ForeSeqs is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ForeSeqs is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ForeSeqs.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DNAMODEL_H_
#define DNAMODEL_H_

#include "Model.h"

#define NUM_NUC 4
#define SQNUM_NUC 4*4
#define CUNUM_NUC 4*4*4

namespace foreseqs {

class DnaModel: public Model {
public:

	/**
	 * @brief Constructs a new DNA model
	 *
	 * @param pllPartitions The PLL partition list
	 * @param partitionIndex The current partition
	 */
	DnaModel(pll_partition_t * partition);

	virtual void setMatrix(double * pMatrix, double branchLength, bool cummulative = true);
	virtual char getState(const double * pMatrix) const;
	virtual char getMostProbableState(const double * probArray) const;
	virtual ~DnaModel();
private:

	/**
	 * @brief Setup the Q matrix
	 */
	void SetupGTR( void );

	double _Cijk[CUNUM_NUC];        /** The cumulative transition probability matrix */
};

} /* namespace foreseqs */

#endif /* DNAMODEL_H_ */
