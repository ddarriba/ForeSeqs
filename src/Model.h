/*
 * Model.h
 *
 *  Created on: Oct 1, 2014
 *      Author: diego
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "pll.h"

#include <map>
#include <vector>

namespace seqpred {

class Model {
public:
	Model(partitionList * pllPartitions, int partitionIndex);

	/**
	 * @brief Compute the P-Matrix for selecting the character to insert
	 */
	virtual void setMatrix(double * matrix, double branchLength) const = 0;
	virtual char getState(double * Pmatrix) const = 0;

	int getStateIndex(char state) {
		return statesMap[state];
	}

	virtual ~Model();
protected:
	double computeFracchange( void ) const;
	pInfo * partitionInfo;				/** PLL Partition info */
	std::vector<double> frequencies;	/** Frequencies */
	std::vector<double> substRates;		/** Substitution rates */
	std::vector<char> states;			/** Vector of the different states */
	std::map<char, int> statesMap;		/** Map of the states index according to char */
};

} /* namespace seqpred */

#endif /* MODEL_H_ */
