/*
 * Model.h
 *
 *  Created on: Oct 1, 2014
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

#ifndef MODEL_H_
#define MODEL_H_

#include "PllDefs.h"

#include <map>
#include <vector>

#define MIN_BRANCH_LEN 1e-6

namespace foreseqs {

class Model {
public:

	/**
	 * @brief Constructs a new model
	 *
	 * @param partition The partition data
	 */
	Model(pll_partition_t * partition);

	/**
	 * @brief Copy constructor
	 */
	Model(const Model&);

	/**
	 * @brief Destructor
	 */
	virtual ~Model();

	/**
	 * @brief Compute the P-Matrix for selecting the character to insert
	 *
	 * @param[out] pMatrix the (already allocated) P matrix
	 * @param[in] branchLength the branch length
	 * @param[in] cummulative if true, the matrix columns are cummulative
	 */
	virtual void setMatrix(double * pMatrix, double branchLength, bool cummulative = true) = 0;

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
	unsigned int getStateIndex(char state) {
		return _statesToIntMap[state];
	}

	/**
	 * @brief Clone operator
	 */
	Model& operator=(const Model&);

protected:

	/**
	 * @brief Compute the P-Matrix for selecting the character to insert
	 *
	 * @param[out] pMatrix the (already allocated) P matrix
	 * @param[in] branchLength the branch length
	 * @param[in] cummulative if true, the matrix columns are cummulative
	 * @param[in] numStates the number of states
	 * @param[in] eigenValues the EIGEN values
	 * @param[in] Cijk the probabilities
	 */
	void constructPMatrix(double * matrix,
                        double branchLength,
 											  bool cummulative,
 											  size_t numStates);

  pll_partition_t * _partition;
	size_t _numberOfStates;               /** Number of states */
	std::vector<double> _frequencies;     /** Frequencies */
	std::vector<double> _substRates;      /** Substitution rates */
	std::vector<char> _charStates;        /** Vector of the different states */
	std::map<char, unsigned int> _statesToIntMap;  /** Map of the states index according to char */

	double * _eigenVals;
	double * _eigenVecs;
	double * _invEigenVecs;
};

} /* namespace foreseqs */

#endif /* MODEL_H_ */
