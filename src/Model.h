/*
 * Model.h
 *
 *  Created on: Oct 1, 2014
 *      Author: diego
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "pll/pll.h"

#include <map>
#include <vector>

namespace seqpred {

class Model {
public:

	Model(partitionList * pllPartitions, int partitionIndex);
	Model(const Model&);

	virtual ~Model();

	/**
	 * @brief Compute the P-Matrix for selecting the character to insert
	 *
	 * @param[out] pMatrix the (already allocated) P matrix
	 * @param[in] branchLength the branch length
	 * @param[in] cummulative if true, the matrix columns are cummulative
	 */
	virtual void setMatrix(double * pMatrix, double branchLength, bool cummulative = true) const = 0;

	/**
	 * @brief Get the randomly selected state according to a probability matrix
	 *
	 * @param[in] pMatrix the commulative P matrix
	 * @return the randomly selected char
	 */
	virtual char getState(const double * pMatrix) const = 0;

	/**
	 * @brief Get the most probable state according to a probability array
	 *
	 * @param[in] probArray the probability array
	 * @return the randomly selected char
	 */
	virtual char getMostProbableState(const double * probArray) const = 0;

	/**
	 * @brief Get the index for a state character
	 *
	 * @param[in] state the state char
	 * @return the index of the state
	 */
	int getStateIndex(char state) {
		return _statesToIntMap[state];
	}

	Model& operator=(const Model&);

protected:

	/**
	 * @brief Compute the fracchange value according to frequencies and substRates
	 *
	 * @return the fracchange
	 */
	double computeFracchange( void ) const;

	pInfo * _pllPartitionInfo;				/** PLL Partition info */
	std::vector<double> _frequencies;		/** Frequencies */
	std::vector<double> _substRates;		/** Substitution rates */
	std::vector<char> _charStates;			/** Vector of the different states */
	std::map<char, int> _statesToIntMap;	/** Map of the states index according to char */
};

} /* namespace seqpred */

#endif /* MODEL_H_ */
