/*
 * Model.cpp
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

#include "Model.h"
#include "Utils.h"

#include <cmath>
#include <cassert>
#include <alloca.h>

using namespace std;

namespace foreseqs {

Model::Model(pll_partition_t * partition) :
	  _partition(partition),
		_numberOfStates(partition->states),
		_frequencies(), _substRates(), _charStates(), _statesToIntMap() {
}

Model::Model(const Model& other) :
	  _partition(other._partition),
		_numberOfStates(other._numberOfStates),
		_frequencies(other._frequencies), _substRates(other._substRates),
		_charStates(other._charStates), _statesToIntMap(other._statesToIntMap) {

}

Model::~Model() {
	/* nothing to do */
}

Model& Model::operator=(const Model& other) {
	_partition = other._partition;
	_numberOfStates = other._numberOfStates;

	assert(_partition->states = _numberOfStates);

	_frequencies = other._frequencies;
	_substRates = other._substRates;
	_charStates = other._charStates;
	_statesToIntMap = other._statesToIntMap;
	return *this;
}

void Model::constructPMatrix(double * matrix,
	                           double branchLength,
														 bool cummulative,
														 size_t numStates)
{

	unsigned int n,j,k,m;
	double * expd;
	double * temp;

  double * evecs;
  double * inv_evecs;
  double * evals;
  double * pmat;

	expd = (double *)malloc(numStates * sizeof(double));
	temp = (double *)malloc(numStates*numStates*sizeof(double));

	assert(branchLength >= 0);

	/* compute effective pmatrix location */
	for (n = 0; n < _partition->rate_cats; ++n)
	{
		double * evecs = _eigenVecs;
		double * inv_evecs = _invEigenVecs;
		double * evals = _eigenVals;

		/* if branch length is zero then set the p-matrix to identity matrix */
		if (!branchLength)
		{
			for (j = 0; j < numStates; ++j)
				for (k = 0; k < numStates; ++k)
					matrix[j*numStates + k] = (j == k) ? 1 : 0;
		}
		else
		{
			/* exponentiate eigenvalues */
			for (j = 0; j < numStates; ++j)
				expd[j] = exp(evals[j] * _partition->rates[n] * branchLength);

			for (j = 0; j < numStates; ++j)
				for (k = 0; k < numStates; ++k)
					temp[j*numStates+k] = inv_evecs[j*numStates+k] * expd[k];

			for (j = 0; j < numStates; ++j)
			{
				for (k = 0; k < numStates; ++k)
				{
					matrix[j*numStates+k] = 0;
					for (m = 0; m < numStates; ++m)
					{
						matrix[j*numStates+k] +=
								temp[j*numStates+m] * evecs[m*numStates+k];
					}
				}
			}
		}
	}

	if (cummulative) {
			/* the rows are cumulative to help with picking one using
			 a random number */
			for (size_t i = 0; i < numStates; i++) {
				for (size_t j = 1; j < numStates; j++) {
					size_t nextIndex = numStates * i + j;
					matrix[nextIndex] += matrix[nextIndex - 1];
				}
				assert(Utils::floatEquals(matrix[numStates * (i + 1) - 1], 1.0));
			}
		} else {
			/* the matrix rows sum to 1.0 */
			for (size_t i = 0; i < numStates; i++) {
				double sum = 0.0;
				for (size_t j = 0; j < numStates; j++) {
					sum += matrix[numStates * i + j];
				}
				assert(Utils::floatEquals(sum, 1.0));
			}
		}

	free(expd);
	free(temp);
}

//
// 	double expt[numStates];
// 	double *P;
// 	size_t squaredNumStates = numStates * numStates;
//
// 	P = matrix;
//
// 	if (branchLength < MIN_BRANCH_LEN) {
// 		/* account for near zero branch lengths */
// 		for (size_t i = 0; i < numStates; i++) {
// 			for (size_t j = 0; j < numStates; j++) {
// 				if (cummulative)
// 					*P = (i<=j)?1.0:0.0;
// 				else
// 					*P = (i==j)?1.0:0.0;
// 				P++;
// 			}
// 		}
// 		return;
// 	}
// 	else
// 	{
// 		for (size_t k = 0; k < numStates; k++) {
// 			expt[k] = exp(branchLength * eigenValues[k]);
// 			cout << " [EI " << k << "] " << eigenValues[k] << " " << expt[k] << endl;
// 		}
// 		for (size_t i = 0; i < numStates; i++) {
// 			for (size_t j = 0; j < numStates; j++) {
// 				(*P) = Cijk[i * squaredNumStates + j * numStates + 0];
// 				for (size_t k = 1; k < numStates; k++) {
// 					(*P) += Cijk[i * squaredNumStates + j * numStates + k]
// 							* expt[k];
// 				}
// 				P++;
// 			}
// 		}
// 	}
// 	if (cummulative) {
// 		/* the rows are cumulative to help with picking one using
// 		 a random number */
// 		for (size_t i = 0; i < numStates; i++) {
// 			for (size_t j = 1; j < numStates; j++) {
// 				size_t nextIndex = numStates * i + j;
// 				matrix[nextIndex] += matrix[nextIndex - 1];
// 				cout << " [" << nextIndex << "] " << matrix[nextIndex];
// 			}
// 			cout << endl;
// 			assert(Utils::floatEquals(matrix[numStates * (i + 1) - 1], 1.0));
// 		}
// 	} else {
// 		/* the matrix rows sum to 1.0 */
// 		for (size_t i = 0; i < numStates; i++) {
// 			double sum = 0.0;
// 			for (size_t j = 0; j < numStates; j++) {
// 				sum += matrix[numStates * i + j];
// 			}
// 			assert(Utils::floatEquals(sum, 1.0));
// 		}
// 	}
// }

} /* namespace foreseqs */
