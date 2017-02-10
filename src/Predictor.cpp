/*
 * Predictor.cpp
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

#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <alloca.h>

#include "Predictor.h"
#include "DnaModel.h"
#include "ProteinModel.h"
#include "Utils.h"

#ifdef _USE_OPENBLAS_
#include <cblas.h>
#endif

#define PRINT_TRACE

using namespace std;

#define MIN_BR_LEN 0.00001

namespace foreseqs {

	static char * getAncestral(pll_partition_t * partition,
		                         pll_utree_t * node,
														 Model * model)
	{
		char * ancestral = (char *) malloc (partition->sites + 1);
		char * anc_ptr;
		double * clvp;

		anc_ptr = ancestral;
		clvp = partition->clv[node->clv_index];
		for (int i=0; i<partition->sites; i++)
		{
			*(anc_ptr++) = model->getMostProbableState(clvp);
			clvp += partition->states * partition->rate_cats;
		}
		*anc_ptr = '\0';

		return ancestral;
	}

Predictor::Predictor(pll_utree_t * tree,
					Phylo & phylo,
					Alignment & alignment,
					unsigned int partitionNumber) :
		    _tree(tree),
				_partition(phylo.getPartition(partitionNumber)),
		    _alignment(alignment),
				_phylo(phylo),
				_partitionNumber(partitionNumber), _numberOfStates(0),
				_start(alignment.getStartPosition(partitionNumber)),
				_end(alignment.getEndPosition(partitionNumber)),
				_partitionLength(_end - _start + 1), _catToSite(0),
				_missingSubtreesAncestors(),
				_missingSequences(phylo.getMissingSequences(partitionNumber)),
				_missingBranches(phylo.getMissingBranches(partitionNumber)),
				_missingPartsCount(0),
				_currentModel(0),
				_seqSimilarity(0),
				_branchLengthScaler(1.0),
				_scalers() {

	if (alignment.getDataType(partitionNumber) == DT_NUCLEIC) {
		_currentModel = new DnaModel(_partition);
		_numberOfStates = 4;
	} else {
		_currentModel = new ProteinModel(_partition);
		_numberOfStates = 20;
	}

	/* get taxa with missing data in the partition */
	_missingPartsCount = (unsigned int)_missingSequences.size();

	if (_missingSequences.size()) {
		if (categoriesMode != CAT_AVERAGE) {
			_catToSite = (short *) calloc((size_t) _partitionLength, sizeof(short));
			unsigned int * catToSiteCount;
			catToSiteCount = (unsigned int *) alloca (numberOfRateCategories*sizeof (unsigned int));
			for (size_t k = 0; k < numberOfRateCategories; k++) {
				catToSiteCount[k] = 0;
			}
			if (categoriesMode == CAT_RANDOM) {
				/* random assignment of sites to categories */
				for (unsigned int i = 0; i < _partitionLength; i++) {
					_catToSite[i] = (short) (numberOfRateCategories * Utils::genRand());
					catToSiteCount[_catToSite[i]]++;
				}
			} else {
				double * gammaRates;
				gammaRates = (double *) alloca (numberOfRateCategories * sizeof(double));
				double * perSiteLikelihoods = (double *) malloc(
						(size_t) _partitionLength * sizeof(double));
				catToSiteCount[0] = _partitionLength;
				if (!perSiteLikelihoods) {
					cerr << "ERROR: There was a problem allocating memory."
							<< endl;
					exit(EX_MEMORY);
				}
				memcpy(gammaRates,
						_partition->rates,
						_partition->rate_cats * sizeof(double));

				/* initialize the per-site likelihoods to a lower bound */
				for (size_t i = 0; i < _partitionLength; i++) {
					perSiteLikelihoods[i] = -10000;
				}

				double * iclv = (double *) malloc(_partition->states * _partitionLength * sizeof(double));
				for (short k = 0; k < (short)numberOfRateCategories; k++) {
					double * clvp = _partition->clv[_tree->clv_index] + k * _partition->states;
					for (int s = 0; s<_partitionLength; s++) {
						iclv[s] = *clvp;
						clvp += _partition->states;
					}
					unsigned int freqsIndex = 0;
					pll_core_edge_loglikelihood_ti(_partition->states,
					                              _partition->sites,
					                                      1,
					                                      iclv,
																								_partition->scale_buffer[_tree->scaler_index],
																							  _partition->tipchars[_tree->back->node_index],
																								_partition->tipmap,
																								_partition->maxstates,
			                                          _partition->pmatrix[_tree->pmatrix_index],
			                                          _partition->frequencies,
			                                          _partition->rate_weights,
			                                          _partition->pattern_weights,
			                                          _partition->prop_invar,
			                                          _partition->invariant,
			                                          &freqsIndex,
			                                          perSiteLikelihoods,
			                                          _partition->attributes);

					//TODO: USE ONE SINGLE CATEGORY

					// /* set all gamma rates to the same value */
					// for (unsigned int i = 0; i < numberOfRateCategories; i++) {
					// 	_partition->rates[i] =
					// 			gammaRates[k];
					// }
					//
					// /* get per-site likelihood */
					// pllEvaluateLikelihood(_tree, _partition,
					// 		_tree->start, true, true);
					//
					// double * X =
					// 		_partition->partitionData[partitionNumber]->perSiteLikelihoods;
					// double * Y = perSiteLikelihoods;
					// for (size_t i = 0; i < _partitionLength; i++) {
					// 	/* if likelihood improves, set this category */
					// 	if (*X > *Y) {
					// 		*Y = *X;
					// 		catToSiteCount[_catToSite[i]]--;
					// 		catToSiteCount[k]++;
					// 		_catToSite[i] = k;
					// 	}
					// 	X++;
					// 	Y++;
					// }
				}
				free(iclv);
				free(perSiteLikelihoods);

				// /* reset rates */
				// memcpy(
				// 		_partition->rates,
				// 		gammaRates, numberOfRateCategories * sizeof(double));
				// pllEvaluateLikelihood(_tree, _partition, _tree->start,
				// 		true, true);

			}

			/* print rates assignment summary */
			cout.setf(ios_base::fixed, ios_base::floatfield);
			cout << "Per site Gamma rate categories assignment:" << endl;
			for (unsigned int k = 0; k < numberOfRateCategories; k++) {
				cout << "  " << k << ": " << setprecision(5)
						<< _partition->rates[k]
						<< " (" << setprecision(2)
						<< (double) 100.0 * catToSiteCount[k] / _partitionLength
						<< "%)" << endl;
			}
			cout << endl;

		}
	}
}

Predictor::Predictor(const Predictor& other) :
				_tree(other._tree), _partition(other._partition),
						_alignment(other._alignment),
						_phylo(other._phylo),
						_partitionNumber(other._partitionNumber),
						_numberOfStates(other._numberOfStates),
						_start(other._start),_end(other._end),
						_partitionLength(other._partitionLength), _catToSite(other._catToSite),
						_missingSubtreesAncestors(other._missingSubtreesAncestors),
						_missingSequences(other._missingSequences),
						_missingBranches(other._missingBranches),
						_missingPartsCount(other._missingPartsCount),
						_currentModel(other._currentModel),	_seqSimilarity(other._seqSimilarity),
						_branchLengthScaler(other._branchLengthScaler),
						_scalers(other._scalers) {
	assert(_end > _start);
	assert(_partitionLength == _end - _start);
}

Predictor::~Predictor( void ) {
	delete _currentModel;
	if (_catToSite) {
		free(_catToSite);
	}
}

Predictor& Predictor::operator=(const Predictor& other) {
	_tree = other._tree;
	_partition = other._partition;
	_alignment = other._alignment;
	_partitionNumber = other._partitionNumber;
	_numberOfStates = other._numberOfStates;
	_start = other._start;
	_end = other._end;
	_partitionLength = other._partitionLength;
	_catToSite = other._catToSite;
	_missingSequences = other._missingSequences;
	_missingPartsCount = other._missingPartsCount;
	_currentModel = other._currentModel;
	_seqSimilarity = other._seqSimilarity;

	assert(_partitionLength == _end - _start);
	assert(_end > _start);

	return *this;
}

void Predictor::mutateSequence(char * currentSequence,
		const char * ancestralSequence, double branchLength) {

	char * seqPtr;
	int substitutionsCount = 0;

	if ( _partitionLength != strlen(ancestralSequence) ) {
		cerr << "ERROR: Length of ancestral sequence (" << strlen(ancestralSequence)
				<< ") differ from the expected length (" << _partitionLength << ")" << endl;
		assert(0);
	}

#ifdef PRINT_TRACE
	cout << "TRACE: Simulating sequence..." << endl;
#endif

	/* start mutating the ancestral sequence */
	strcpy(currentSequence, ancestralSequence);

	seqPtr = currentSequence;

	/* construct and validate P matrix */
	double ** matrix = new double*[numberOfRateCategories];
	for (size_t i=0; i<numberOfRateCategories; i++) {
		matrix[i] = new double[_numberOfStates*_numberOfStates];
	}

	for (unsigned int i = 0; i < numberOfRateCategories; i++) {
		/* compute cummulative P matrix */
		_currentModel->setMatrix(matrix[i], _partition->rates[i] * branchLength, true);

		/* verify that last per-row value is 1.0 */
		for (unsigned int j = 0; j < _numberOfStates; j++) {
			if (!Utils::floatEquals(
					matrix[i][j * _numberOfStates + _numberOfStates - 1],
					1.0)) {
				cerr << "ERROR: Pmatrix " << i << " sums to "
						<< matrix[i][j * _numberOfStates + _numberOfStates - 1]
						<< " instead of 1.0" << endl;
				assert(0);
			}
		}
	}

	switch (categoriesMode) {
	case CAT_AVERAGE:
	{
		double * averageMatrix = (double *) alloca (_numberOfStates * _numberOfStates * sizeof(double));
		for (unsigned int i = 0; i < _numberOfStates * _numberOfStates; i++) {
			averageMatrix[i] = 0.0;
			for (unsigned int j = 0; j < numberOfRateCategories; j++) {
				averageMatrix[i] += matrix[j][i];
			}
			averageMatrix[i] /= numberOfRateCategories;
			if (!((i+1) % _numberOfStates)) {
				assert(Utils::floatEquals(averageMatrix[i], 1.0));
			} else {
				assert(averageMatrix[i] < 1.0);
			}
		}

		for (unsigned int i = 0; i < _partitionLength; i++) {
			char newState = _currentModel->getState(
					averageMatrix
					+ (_currentModel->getStateIndex(*seqPtr)
									* _numberOfStates));
			if (*seqPtr != newState) substitutionsCount++;
			*seqPtr = newState;
			seqPtr++;
		}
		break;
	}
	case CAT_ESTIMATE:
	case CAT_RANDOM:
	{
		short * siteCatPtr = _catToSite;
		for (unsigned int i = 0; i < _partitionLength; i++) {
			char newState = _currentModel->getState(
					matrix[*siteCatPtr]
							+ (_currentModel->getStateIndex(*seqPtr)
									* _numberOfStates));
			if (*seqPtr != newState) substitutionsCount++;
				*seqPtr = newState;
			seqPtr++;
			siteCatPtr++;
		}
		break;
	}
	}

	for (size_t i=0; i<numberOfRateCategories; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;

	cout << "    - Substitutions from ancestral = " << setprecision(2) << (double)100.0*substitutionsCount/_partitionLength << "%" << endl;
}

#ifdef DEBUG
void printNodes(pllInstance * _tree, partitionList * _partition) {
	/* compute the branch length as the average over all other partitions */
	//TODO
	// for (int node = 0; node < 2 * _alignment.getSequenceCount() - 3; node++) {
	// 	cout << node << " " << _tree->nodep[node]->next->back->clv_index << " "
	// 			<< _tree->nodep[node]->next->next->back->clv_index << endl;
	// }
}

void printBranchLengths(pllInstance * _tree, partitionList * _partition) {
	/* compute the branch length as the average over all other partitions */
	//TODO
	// for (int i=0; i<_alignment.getNumberOfPartitions(); i++) {
	// 	cout << i << " -> ";
	// 	for (int node=0; node < 2*_alignment.getSequenceCount()a()-3; node++) {
	// 		cout << " " << getBranchLength(_tree->nodep[node], i);
	// 	}
	// 	cout << endl;
	// }
}
#endif

double Predictor::drawBranchLengthScaler( void ) const {
	double r = Utils::genRand();
	size_t j;
	branchInfo bInfo = _scalers[0];
	for (j=0; r>bInfo.weight; j++) bInfo=_scalers[j+1];
	return (_scalers[j].scaler);
}

static double getBranchLength(const pll_utree_t * node, unsigned int partitionId)
{
	return ((double *)node->data)[partitionId];
}

static void setBranchLength(pll_utree_t * node, unsigned int partitionId, double branchLength)
{
	((double *)node->data)[partitionId] =
		node->length = node->back->length = branchLength;
}

double Predictor::computeBranchLength(pll_utree_t * node) const {

	double branchLength = 0.0;

	if (isMissingBranch(node, _partitionNumber) && _alignment.getNumberOfPartitions() > 1) {
		int sumwgt = 0;
cout << " Branch " << _partitionNumber << " " << node->node_index << endl;

		/* compute the total weight */
		for (unsigned int i = 0;
						i < (unsigned int) _alignment.getNumberOfPartitions(); i++) {
			if (!isMissingBranch(node, i)) {
				cout << " Sum: " << i << " " << _alignment.getSequenceLength(i) << endl;
				sumwgt += _alignment.getSequenceLength(i);
			}
		}

		if (sumwgt == 0) {
			cerr << "Error: There is no information available for stealing branch length ("
					<< node->clv_index << "," << node->back->clv_index << ")" << endl;
		}

		assert(sumwgt > 0);

		double sumFactors = 0.0;
		/* compute the branch length as the average over all other partitions */
		for (unsigned int i = 0;
				i < (unsigned int) _alignment.getNumberOfPartitions(); i++) {
			if (!isMissingBranch(node, i)) {
				double factor = (double) _alignment.getSequenceLength(i) / (double) sumwgt;
				sumFactors += factor;
				double currentBL = getBranchLength(node, i) * factor;
				branchLength += currentBL;
			}
		}
		assert(Utils::floatEquals(sumFactors, 1.0));
	} else {
		/* return the current and only branch length */
		branchLength = getBranchLength(node, _partitionNumber);
	}

	switch (branchLengthsMode) {
	case BL_AVERAGE:
	case BL_SCALE: {
		branchLength *= _branchLengthScaler;
		break;
	}
	case BL_DRAW: {
		double scaler = drawBranchLengthScaler();
		branchLength *= scaler;
		cout << "  - Branch length scaler drawn: " << scaler << endl;
		break;
	}
	}

	/* branch length correction */
	branchLength = max(branchLength, MIN_BR_LEN);
	setBranchLength(node, _partitionNumber, branchLength);

	return branchLength;
}

/* WMOD=0 -> No correction-for-distance */
#if(CFD)
#define WMOD 0.5
#else
#define WMOD 0
#endif
#define WMEAN 0
#define WSIGMA 1

static double computeWeight(double x) {
	double w = 1/(WSIGMA*sqrt(2*3.1415926))*exp(-((x-WMEAN)*(x-WMEAN))/(2*WSIGMA*WSIGMA));
	return w;
}

bool Predictor::isAncestor(pll_utree_t * node) const {
	if (find(_missingSubtreesAncestors.begin(),
				_missingSubtreesAncestors.end(), node)
				== _missingSubtreesAncestors.end())
		return false;
	return true;
}

bool Predictor::isMissingBranch(pll_utree_t * node, size_t partition) const {
	return _phylo.isBranchMissing(node, partition);
}

void Predictor::getBranches(pll_utree_t * node, int depth, double * cumWeight, vector<branchInfo> & branches) const {

	double localSum = 0.0, ratio = 0.0, curWeight = 0.0;
	int sumWgt = 0;

	if (!(isAncestor(node->next) || isAncestor(node->next->next))) {

		/* add to scaler only if the branch exist */
		if (!isMissingBranch(node, _partitionNumber)) {

			for (unsigned int i = 0;
					i < (unsigned int) _alignment.getNumberOfPartitions();
					i++) {
				if (_partitionNumber != i && !isMissingBranch(node, i)) {
					sumWgt += _alignment.getSequenceLength(i);
				}
			}

			/* add to scaler only if there is info in other partitions */
			if (sumWgt > 0) {

				double branchLength = getBranchLength(node, _partitionNumber);

				/* compute the current weighted branch length ratio */
				double sumFactors = 0.0;
				for (unsigned int i = 0;
						i < (unsigned int) _alignment.getNumberOfPartitions();
						i++) {
					if (_partitionNumber != i && !isMissingBranch(node, i)) {
						double curBranchLength = getBranchLength(node, i);
						double factor =
								(double) _alignment.getSequenceLength(i)
										/ (double) sumWgt;
						sumFactors += factor;
						double r = branchLength / curBranchLength;
						if (r > MAX_LOCAL_SCALER || r < (1 / MAX_LOCAL_SCALER)) {
							/* skip */
							localSum += branchLength * factor;
						} else {
							localSum += curBranchLength * factor;
						}
					}
				}
				assert(Utils::floatEquals(sumFactors, 1.0));

				curWeight = computeWeight(depth * WMOD);
				ratio = branchLength / localSum;

				stringstream ss;
				ss << "(" << node->clv_index << "," << node->back->clv_index << ")";
				cout << setw(10) << right << ss.str() << " - weight: "
						<< curWeight << " ratio: " << ratio << endl;

				/* add the new sample */
				branchInfo bInfo;
				bInfo.branchNumber = (size_t) node->pmatrix_index;
				bInfo.scaler = ratio;
				bInfo.weight = curWeight;
				branches.push_back(bInfo);

			}
		}

		if ((unsigned int) node->clv_index > numberOfTaxa) {
			/* inspect children */
			double child1W = 0.0, child2W = 0.0;
			getBranches(node->next->back, depth+1, &child1W, branches);
			getBranches(node->next->next->back, depth+1, &child2W, branches);
			*cumWeight += child1W + child2W;
		}
		*cumWeight += curWeight;
	}
}

double Predictor::getSumBranches(pll_utree_t * node, int depth, double * weight) const {

	double childrenSum = 0.0, localSum = 0.0, ratio = 0.0;

	/* we define 2 weights: one according to the distance to the rooting branch
	 * and another one accorting to the number of sites in the partition.
	 */
	double distanceWeight = 0.0;
	int sumWgt = 0;

	if (!(isAncestor(node->next) || isAncestor(node->next->next))) {

		/* add to scaler only if the branch exist */
		if (!isMissingBranch(node, _partitionNumber)) {

			for (unsigned int i = 0;
					i < (unsigned int) _alignment.getNumberOfPartitions();
					i++) {
				if (_partitionNumber != i && !isMissingBranch(node, i)) {
					sumWgt += _alignment.getSequenceLength(i);
				}
			}

			/* add to scaler only if there is info in other partitions */
			if (sumWgt > 0) {

				double branchLength = getBranchLength(node, _partitionNumber);

				/* compute the current weighted branch length ratio */
				double sumFactors = 0.0;
				for (unsigned int i = 0;
						i < (unsigned int) _alignment.getNumberOfPartitions();
						i++) {

					if (_partitionNumber != i && !isMissingBranch(node, i)) {

						double curBranchLength = getBranchLength(node, i);
						double factor =
								(double) _alignment.getSequenceLength(i)
										/ (double) sumWgt;
						sumFactors += factor;
//						double r = branchLength / curBranchLength;
//						if (r > MAX_LOCAL_SCALER
//								|| r < (1.0 / MAX_LOCAL_SCALER)) {
//							/* skip */
//							localSum += branchLength * factor;
//						} else {
							localSum += curBranchLength * factor;
//						}

					}

				}
				assert(Utils::floatEquals(sumFactors, 1.0));

				distanceWeight = computeWeight(depth * WMOD);
				ratio = branchLength / localSum;

				stringstream ss;
				ss << "(" << node->clv_index << "," << node->back->clv_index << ")";
				cout << setw(10) << right << ss.str() << setprecision(6)
						<< " - weight: " << distanceWeight << " ratio: "
						<< ratio << endl;
				ratio *= distanceWeight;
			}
		}

		if ((unsigned int) node->clv_index > numberOfTaxa) {
			/* inspect children */
			double child1W = 0.0, child2W = 0.0;
			childrenSum += getSumBranches(node->next->back, depth+1, &child1W);
			childrenSum += getSumBranches(node->next->next->back, depth+1, &child2W);
			*weight += child1W + child2W;
		} else {
			assert(*weight == 0.0);
		}
		*weight += distanceWeight;

	}
	return (ratio + childrenSum);
}

//double Predictor::computeBranchLengthScaler( void ) const {
//
//	double scaler = 0.0;
//
//	if (_missingSubtreesAncestors.size()) {
//
//		/* traverse the tree starting at one arbitrary ancestor node */
//		pll_utree_t * startingNode = _missingSubtreesAncestors[0];
//
//		cout << "Branch length scaler:" << endl;
//		double cumWeight = 0;
//		scaler = getSumBranches(startingNode->next->back, 0, &cumWeight);
//		scaler += getSumBranches(startingNode->next->next->back, 0, &cumWeight);
//		scaler = scaler/cumWeight;
//	}
//
//	return scaler;
//}

void Predictor::stealBranchRecursive(pll_utree_t * node) {
	computeBranchLength(node);

	if (node->clv_index > _alignment.getSequenceCount()) {
		stealBranchRecursive(node->next->back);
		stealBranchRecursive(node->next->next->back);
	} else {
		/* remove visited taxon */
		_missingSequences.erase(
				remove(_missingSequences.begin(), _missingSequences.end(),
						node->clv_index), _missingSequences.end());
	}
}

void Predictor::evolveNode(pll_utree_t * node, const double * ancestralProbabilities, double * ancestralPMatrix) {

#ifdef PRINT_TRACE
	cout << "TRACE: Mutating node " << node->clv_index << endl;
#endif
#ifdef PRINT_ANCESTRAL
	if (ancestralSequence) {
		cout << "INFO: Ancestral: ";
	    Utils::printSequence(ancestralSequence);
	}
#endif
	unsigned int n = _partitionLength;
	double branchLength = computeBranchLength(node);
	cout << "  - Estimated branch length for node " << setprecision(5)
			<< node->clv_index << " (partition "
			<< _partitionNumber
			<< ") = " << branchLength << endl;

	if (node->clv_index > _alignment.getSequenceCount()) {

		double * currentPMatrix = (double *) Utils::allocate(
				(size_t) (numberOfRateCategories * _numberOfStates
						* _numberOfStates), sizeof(double));

		mutatePMatrix(currentPMatrix, ancestralPMatrix, branchLength);

		evolveNode(node->next->back, ancestralProbabilities, currentPMatrix);
		evolveNode(node->next->next->back, ancestralProbabilities, currentPMatrix);

		free(currentPMatrix);
	} else {
		/* compute the new probability matrix */
		/* construct the new marginal probabilities matrix */

		// Multiply P' x M
		double * currentProbabilities = (double *) Utils::allocate(
				(size_t) numberOfRateCategories * _numberOfStates * n,
				sizeof(double));

		for (size_t cat = 0; cat < numberOfRateCategories; cat++) {
			if (ancestralPMatrix) {
#ifdef _USE_OPENBLAS_
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					n, _numberOfStates, _numberOfStates,
					1.0, ancestralProbabilities, _numberOfStates,
					&ancestralPMatrix[cat * _numberOfStates* _numberOfStates], _numberOfStates,
					0.0, &currentProbabilities[cat * _numberOfStates * n], _numberOfStates);
#else
				Utils::matrixMultiply(_numberOfStates, n, ancestralProbabilities,
						&ancestralPMatrix[cat * _numberOfStates
								* _numberOfStates], &currentProbabilities[cat * _numberOfStates * n]);
#endif
			} else {
				/* copy the ancestral probabilities for each category */
				memcpy(&currentProbabilities[cat * _numberOfStates * n], ancestralProbabilities,
						_numberOfStates * n * sizeof(double));
			}

#ifdef USE_FIXED_ANCESTRAL
			double * M = &currentProbabilities[cat * _numberOfStates * n];
			for (unsigned int i = 0; i < n; i++) {
				double sum = 0.0;
				for (unsigned int j = 0; j < _numberOfStates; j++) {
					size_t nextIndex = _numberOfStates * i + j;
					sum += M[nextIndex];
				}
				assert(Utils::floatEquals(sum, 1.0));
			}
#else
			/* Make M matrix cummulative */
			double * M = &currentProbabilities[cat * _numberOfStates * n];
			for (unsigned int i = 0; i < n; i++) {
				for (unsigned int j = 1; j < _numberOfStates; j++) {
					size_t nextIndex = _numberOfStates * i + j;
					M[nextIndex] += M[nextIndex - 1];
				}
				assert(
						Utils::floatEquals(M[_numberOfStates * (i + 1) - 1],
								1.0));
			}
#endif
		}

		/* draw random ancestral states */
		char * ancestralSequence = (char *) malloc(n + 1);
		switch (categoriesMode) {
		case CAT_AVERAGE: {
			cerr << "Not implemented yet" << endl;
			exit(EX_UNIMPLEMENTED);
		}
		case CAT_RANDOM:
		case CAT_ESTIMATE: {

			short * siteCatPtr = _catToSite;
			char * seqPtr = ancestralSequence;
#ifdef USE_FIXED_ANCESTRAL
			for (unsigned int i = 0; i < _partitionLength; i++) {
				char newState = _currentModel->getMostProbableState(
						&currentProbabilities[(size_t)(*siteCatPtr)
								* _numberOfStates * _partitionLength
								+ i * _numberOfStates]);
				*seqPtr = newState;
				seqPtr++;
				siteCatPtr++;
			}
#else
			for (unsigned int i = 0; i < _partitionLength; i++) {
				char newState = _currentModel->getState(
						&currentProbabilities[(*siteCatPtr)
								* _numberOfStates * _partitionLength
								+ i * _numberOfStates]);
				*seqPtr = newState;
				seqPtr++;
				siteCatPtr++;
			}
#endif
			break;
		}
		}
		ancestralSequence[n] = '\0';

cout << "ANCESTRAL " << ancestralSequence << endl;
		/* compute the new sequence from the probabilities */
		char * currentSequence = (char *) malloc (n + 1);
		mutateSequence(currentSequence, ancestralSequence, branchLength);
		free(ancestralSequence);

		/* set the new sequence */
		//TODO: Use setter!
		//memcpy(&(_pllAlignment->sequenceData[seqIndexTranslate[node->clv_index]][_start]), currentSequence, _partitionLength);

		/* remove visited taxon */
		_missingSequences.erase(remove(_missingSequences.begin(), _missingSequences.end(), node->clv_index), _missingSequences.end());

		free(currentSequence);
		free(currentProbabilities);
	}
}

void Predictor::mutatePMatrix(double * currentPMatrix,
		double * ancestralPMatrix, double branchLength) {

	/* construct and validate P matrix */
	double * matrix = (double *) Utils::allocate(
			(size_t) numberOfRateCategories * _numberOfStates * _numberOfStates,
			sizeof(double));

	for (unsigned int i = 0; i < numberOfRateCategories; i++) {
		_currentModel->setMatrix(
				&(matrix[i * _numberOfStates * _numberOfStates]),
				_partition->rates[i] * branchLength, false);

		/* validation */
		for (unsigned int j = 0; j < _numberOfStates; j++) {
			double sum = 0.0;
			for (unsigned int k = 0; k < _numberOfStates; k++) {
				sum += matrix[i * _numberOfStates * _numberOfStates + j*_numberOfStates + k];
			}
			assert(Utils::floatEquals(sum, 1.0));
		}
	}

	if (ancestralPMatrix) {
		/* multiply P matrices */
		size_t n = _numberOfStates;
		for (size_t cat = 0; cat < numberOfRateCategories; cat++) {
			double * catMatrix = &(matrix[cat * _numberOfStates
					* _numberOfStates]);
#ifdef _USE_OPENBLAS_
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0,
					catMatrix, n, &ancestralPMatrix[cat * n * n], n, 0.0,
					&currentPMatrix[cat * n * n], n);
#else
			Utils::matrixMultiply(n, n, catMatrix, &ancestralPMatrix[cat * n * n],
					&currentPMatrix[cat * n * n]);
#endif
		}
	} else {
		/* directly set the P matrix */
		for (unsigned int i = 0; i < numberOfRateCategories; i++) {
			memcpy(currentPMatrix, matrix,
					numberOfRateCategories * _numberOfStates * _numberOfStates
							* sizeof(double));
		}
	}

	free(matrix);
}

void Predictor::evolveNode(pll_utree_t * node, const char * ancestralSequence) {

#ifdef PRINT_TRACE
	cout << "TRACE: Mutating node " << node->clv_index << endl;
	cout << "TRACE: From sequence " << ancestralSequence << endl;
#endif
#ifdef PRINT_ANCESTRAL
	if (ancestralSequence) {
		cout << "INFO: Ancestral: ";
	    Utils::printSequence(ancestralSequence);
	}
#endif

	double branchLength = computeBranchLength(node);
	cout << "  - Estimated branch length for node " << setprecision(5) << node->clv_index << " (partition " << _partitionNumber << ") = " << branchLength << endl;

	char * currentSequence = (char *) malloc(strlen(ancestralSequence) + 1);
	mutateSequence(currentSequence, ancestralSequence, branchLength);

	if (node->clv_index > _alignment.getSequenceCount()) {
#ifdef PRINT_TRACE
		cout << "TRACE: Recurse for Mutating node " << node->next->back->clv_index
				<< " and " << node->next->next->back->clv_index << endl;
#endif
		evolveNode(node->next->back, currentSequence);
		evolveNode(node->next->next->back, currentSequence);
	} else {
		/* set the new sequence */
		//TODO: Use setter!
		//memcpy(&(_pllAlignment->sequenceData[seqIndexTranslate[node->clv_index]][_start]), currentSequence, _partitionLength);
		//
		/* remove visited taxon */
		_missingSequences.erase(remove(_missingSequences.begin(), _missingSequences.end(), node->clv_index), _missingSequences.end());
	}
	free(currentSequence);
}

void Predictor::getRootingNodes() {
	assert ( !_missingSubtreesAncestors.size() );


#ifdef DEBUG_BRANCHES
		cout << pll_utree_export_newick(_tree) << endl;
		cout << "  Branches: " << _missingBranches.size() << endl;
		for (int j = 0; j < _missingBranches.size(); j++) {
			cout << "    * " << _missingBranches[j]->clv_index << endl;
		}
#endif

	if ( _missingBranches.size()) {
		unsigned int i;
		for (i = 0;
				i < _missingBranches.size()
						&& _missingBranches[i]->clv_index
								< _alignment.getSequenceCount(); i++)
			;

		pll_utree_t * currentNode = _missingBranches[i++];
		int counter = 0;
		for ( ; i<_missingBranches.size(); i++) {
			int testNode = _missingBranches[i]->clv_index;
			assert ( testNode >=  currentNode->clv_index );
			if ( testNode == currentNode->clv_index ) {
				counter++;
			} else {
				/* each node should appear either 1 or 3 times */
				assert ( counter == 0 || counter == 2);
				if ( counter == 0 ) {
					_missingSubtreesAncestors.push_back(currentNode);
				} else {
					counter = 0;
				}
				currentNode = _missingBranches[i];
			}
		}
		/* account for the last one */
		if (counter == 0) {
			_missingSubtreesAncestors.push_back(currentNode);
		}
	} else {
		return;
	}
}

void Predictor::predictMissingSequences( void ) {

#ifdef PRINT_TRACE
	for (size_t i = 0; i < _missingSequences.size(); i++) {
		cout << "TRACE: Missing sequence " << _missingSequences[i] << endl;
	}
#endif

	vector<unsigned int> missingSequencesCopy;
	// if (originalSequence) {
	// 	missingSequencesCopy = _missingSequences;
	// }

	getRootingNodes();

	if (_missingSequences.size()) {
		if (branchLengthsMode == BL_SCALE) {
			/* compute branch length scaler */
			if (_missingSubtreesAncestors.size()) {
				/* traverse the tree starting at one arbitrary ancestor node */
				pll_utree_t * startingNode = _missingSubtreesAncestors[0];

				cout << "Branch length scaler:" << endl;
				double cumWeight = 0.0, weight = 0.0;
				_branchLengthScaler = getSumBranches(startingNode->next->back, 0,
						&cumWeight);
				_branchLengthScaler += getSumBranches(
						startingNode->next->next->back, 0, &weight);
				cumWeight += weight;

				if (cumWeight > EPSILON) {
					_branchLengthScaler /= cumWeight;
				} else {
					_branchLengthScaler = 1.0;
				}
				_branchLengthScaler = min(_branchLengthScaler, MAX_SCALER);
				_branchLengthScaler = max(_branchLengthScaler, MIN_SCALER);
			}
			cout << "Applying branch length scaler " << _branchLengthScaler << endl << endl;
		} else if (branchLengthsMode == BL_DRAW) {
			/* traverse the tree starting at one arbitrary ancestor node */
			pll_utree_t * startingNode = _missingSubtreesAncestors[0];

			cout << "Branch length scaler:" << endl;
			double cumWeight = 0.0;
			getBranches(startingNode->next->back, 0, &cumWeight, _scalers);
			getBranches(startingNode->next->next->back, 0, &cumWeight, _scalers);
			cout << endl;

			/* normalize the weights */
			_scalers[0].weight /= cumWeight;
			for (size_t i = 1; i < _scalers.size(); i++) {
				_scalers[i].weight /= cumWeight;
				_scalers[i].weight += _scalers[i-1].weight;
			}
		}
	}

	/* loop over all possible subtrees with missing data */
	unsigned int nextAncestor = 0;
	while (_missingSequences.size()) {
		cout << "Predicting subtree" << endl;
		assert(nextAncestor < _missingSubtreesAncestors.size());
		pll_utree_t * ancestor = _missingSubtreesAncestors[nextAncestor++];

		pll_utree_t * startNode = ancestor;

#ifdef PRINT_TRACE
		cout << "TRACE: Generating ancestral for " << startNode->clv_index << endl;
#endif

		/*
		 * Compute probs size according to the maximum number of states.
		 * This is necessary in case there are mixed protein-dna partitions.
		 */

		// int probsSize = 0;
		// for (int i = 0; i < _alignment.getNumberOfPartitions(); i++) {
		// 	probsSize = max(_partition->states,
		// 			probsSize);
		// }
		// probsSize *= sequenceLength;

		// char * ancestral = (char *) malloc((size_t) sequenceLength + 1);
		// double * probs = (double *) malloc((size_t) probsSize * sizeof(double));
		char * ancestral = getAncestral(_partition, startNode, _currentModel);
		cout << "xANCESTRAL = " << ancestral << endl;
		// pllGetAncestralState(_tree, _partition, startNode, probs,
		// 		ancestral, false);

#ifdef DEBUG
		printNodes(_tree, _partition);
		printBranchLengths(_tree, _partition);
#endif

		switch (predictionMode) {
		case PRED_ANCSEQ: {
			/* discard MAP vector */
			// free(probs);

			evolveNode(ancestor->back, ancestral);

			free(ancestral);
			break;
		}
		case PRED_MAP: {
			/* discard ancestral sequence */
			free(ancestral);

			/* extract the ancestral for the partition */

			// size_t partProbsSize = ((size_t) probsSize / sequenceLength)
			// 		* _partitionLength;

			// double * probsAncestral = (double *) Utils::allocate(partProbsSize,
			// 		sizeof(double));
			//
			// memcpy(probsAncestral, &(probs[_start * _numberOfStates]),
			// 		partProbsSize * sizeof(double));
			// free(probs);

			evolveNode(ancestor->back, _partition->clv[ancestor->clv_index]);
			// free(probsAncestral);
		}
			break;
		case PRED_NONE:
			stealBranchRecursive(ancestor->back);
			break;
		}

		// TODO:
		// if (originalSequence) {
		// 	/* compute similarity */
		// 	_seqSimilarity = 0.0;
		// 	for (size_t i = 0; i < missingSequencesCopy.size(); i++) {
		// 		unsigned int seq = missingSequencesCopy[i];
		// 		double simCount = 0.0;
		// 		unsigned int seqLen = 0;
		// 		if (_alignment.getDataType(_partitionNumber)
		// 				== DT_NUCLEIC) {
		// 			for (size_t j = _start; j < _end; j++) {
		// 				bool validForComp;
		// 				simCount +=
		// 						Utils::compareNucStates(
		// 								_pllAlignment->sequenceData[seqIndexTranslate[seq]][j],
		// 								originalSequence->sequenceData[seqIndexTranslate[seq]][j],
		// 								&validForComp);
		// 				seqLen += validForComp;
		// 			}
		// 			cout << "Comparing sequence " << taxaNames[seq] << ": "
		// 					<< simCount << " valid sites." << endl;
		// 			_seqSimilarity += simCount / seqLen
		// 					/ missingSequencesCopy.size();
		// 		}
		// 		cout << "Similarity in sequence " << seq << "/"
		// 				<< _partition->partitionData[_partitionNumber]->partitionName
		// 				<< ": " << 100 * simCount / seqLen << "%" << endl;
		// 	}
		// }
	}

	/* ensure there is no subtree left */
	assert(nextAncestor == _missingSubtreesAncestors.size());

	cout << "Tree with stolen branches for partition " << _partitionNumber
			<< ":" << endl;
	cout << pll_utree_export_newick(_tree) << endl;
}

void Predictor::predictAllSequences( void ) {

	_missingSequences.resize(1);

	//TODO: Extract tip nodes
	pll_utree_t ** tipNodes;

	for (unsigned int taxon=0; taxon < numberOfTaxa; taxon++) {
		_missingSequences[0] = taxon;
		pll_utree_t * ancestor = tipNodes[taxon]->back;
		/*
		 * Compute probs size according to the maximum number of states.
		 * This is necessary in case there are mixed protein-dna partitions.
		 */
		// int probsSize = 0;
		// for (int i = 0; i < _alignment.getNumberOfPartitions(); i++) {
		// 	probsSize = max(_partition->states,
		// 			probsSize);
		// }
		// probsSize *= sequenceLength;

		// char * ancestral = (char *) malloc( (size_t) sequenceLength + 1);
		// double * probs = (double *) malloc(
		// 		(size_t) probsSize * sizeof(double));
		char * ancestral = getAncestral(_partition, ancestor, _currentModel);
		// free (probs);

		evolveNode(ancestor->back, ancestral);
		free(ancestral);
	}
}

} /* namespace foreseqs */
