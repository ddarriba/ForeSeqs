/*
 * Predictor.cpp
 *
 *  Created on: Oct 1, 2014
 *      Author: diego
 */

#include "Predictor.h"
#include "DnaModel.h"
#include "ProteinModel.h"
#include "Utils.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <cmath>
#include <alloca.h>

#include <cblas.h>

using namespace std;

#define MIN_BR_LEN 0.0001

namespace seqpred {

Predictor::Predictor(pllInstance * tree, partitionList * partitions,
		pllAlignmentData * phylip, unsigned int partitionNumber,
		std::vector<int> missingSequences,
		const std::vector< std::vector<nodeptr> > * missingBranches) :
		_pllTree(tree), _pllPartitions(partitions), _pllAlignment(phylip),
				_partitionNumber(partitionNumber), _numberOfStates(0),
				_start((unsigned int)partitions->partitionData[partitionNumber]->lower),
				_end((unsigned int)partitions->partitionData[partitionNumber]->upper),
				_partitionLength(_end - _start), _catToSite(0),
				_missingSubtreesAncestors(),
				_missingSequences(missingSequences), _missingBranches(missingBranches),
				_missingPartsCount(0),
				_currentModel(0), _seqSimilarity(0),
				_branchLengthScaler(1.0), _scalers() {

	if (Utils::getDataType(partitions, partitionNumber) == DT_NUCLEIC) {
		_currentModel = new DnaModel(partitions, partitionNumber);
		_numberOfStates = 4;
	} else {
		_currentModel = new ProteinModel(partitions, partitionNumber);
		_numberOfStates = 20;
	}

	/* get taxa with missing data in the partition */
	_missingPartsCount = _missingSequences.size();

	if (_missingSequences.size()) {
		if (categoriesMode != CAT_AVERAGE) {
			_catToSite = (short *) calloc((size_t) _partitionLength, sizeof(short));
			unsigned int * catToSiteCount;
			catToSiteCount = (unsigned int *) alloca (numberOfRateCategories*sizeof (unsigned int));
			for (unsigned int k = 0; k < numberOfRateCategories; k++) {
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
						_pllPartitions->partitionData[partitionNumber]->gammaRates,
						numberOfRateCategories * sizeof(double));

				/* initialize the per-site likelihoods to a lower bound, e.g., the global likelihood */
				for (unsigned int i = 0; i < _partitionLength; i++) {
					perSiteLikelihoods[i] = _pllTree->likelihood;
				}

				for (unsigned int k = 0; k < numberOfRateCategories; k++) {

					/* set all gamma rates to the same value */
					for (unsigned int i = 0; i < numberOfRateCategories; i++) {
						_pllPartitions->partitionData[partitionNumber]->gammaRates[i] =
								gammaRates[k];
					}

					/* get per-site likelihood */
					pllEvaluateLikelihood(_pllTree, _pllPartitions,
							_pllTree->start, true, true);

					double * X =
							_pllPartitions->partitionData[partitionNumber]->perSiteLikelihoods;
					double * Y = perSiteLikelihoods;
					for (unsigned int i = 0; i < _partitionLength; i++) {
						/* if likelihood improves, set this category */
						if (*X > *Y) {
							*Y = *X;
							catToSiteCount[_catToSite[i]]--;
							catToSiteCount[k]++;
							_catToSite[i] = k;
						}
						X++;
						Y++;
					}
				}

				free(perSiteLikelihoods);

				/* reset rates */
				memcpy(
						_pllPartitions->partitionData[partitionNumber]->gammaRates,
						gammaRates, numberOfRateCategories * sizeof(double));
				pllEvaluateLikelihood(_pllTree, _pllPartitions, _pllTree->start,
						true, true);

			}

			/* print rates assignment summary */
			cout.setf(ios_base::fixed, ios_base::floatfield);
			cout << "Per site Gamma rate categories assignment:" << endl;
			for (unsigned int k = 0; k < numberOfRateCategories; k++) {
				cout << "  " << k << ": " << setprecision(5)
						<< _pllPartitions->partitionData[partitionNumber]->gammaRates[k]
						<< " (" << setprecision(2)
						<< (double) 100.0 * catToSiteCount[k] / _partitionLength
						<< "%)" << endl;
			}
			cout << endl;

		}
	}
}

Predictor::Predictor(const Predictor& other) :
				_pllTree(other._pllTree), _pllPartitions(other._pllPartitions),
						_pllAlignment(other._pllAlignment),
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
	_pllTree = other._pllTree;
	_pllPartitions = other._pllPartitions;
	_pllAlignment = other._pllAlignment;
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
	double * gammaRates = (double *) alloca (numberOfRateCategories * sizeof(double));
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

	pllMakeGammaCats(_pllPartitions->partitionData[_partitionNumber]->alpha, gammaRates, numberOfRateCategories, false);

	seqPtr = currentSequence;

	/* construct and validate P matrix */
	double matrix[numberOfRateCategories][_numberOfStates*_numberOfStates];
	for (unsigned int i = 0; i < numberOfRateCategories; i++) {
		_currentModel->setMatrix(matrix[i], gammaRates[i] * branchLength);
		for (unsigned int j = 0; j < _numberOfStates; j++) {
			assert(Utils::floatEquals(matrix[i][j*_numberOfStates + _numberOfStates - 1], 1.0));
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
	case CAT_RANDOM:
	case CAT_ESTIMATE:
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
	cout << "    - Substitutions from ancestral = " << setprecision(2) << (double)100.0*substitutionsCount/_partitionLength << "%" << endl;
}

#ifdef DEBUG
void printNodes(pllInstance * _pllTree, partitionList * _pllPartitions) {
	/* compute the branch length as the average over all other partitions */
	for (int node = 1; node <= 2 * _pllTree->mxtips - 3; node++) {
		cout << node << " " << _pllTree->nodep[node]->next->back->number << " "
				<< _pllTree->nodep[node]->next->next->back->number << endl;
	}
}

void printBranchLengths(pllInstance * _pllTree, partitionList * _pllPartitions) {
	/* compute the branch length as the average over all other partitions */
	for (int i=0; i<_pllPartitions->numberOfPartitions; i++) {
		cout << i << " -> ";
		for (int node=1; node <= 2*_pllTree->mxtips-3; node++) {
			cout << " " << pllGetBranchLength(_pllTree, _pllTree->nodep[node], i);
		}
		cout << endl;
	}
}
#endif

double Predictor::drawBranchLengthScaler( void ) const {
	double r = Utils::genRand();
	int j=0;
	branchInfo bInfo = _scalers[0];
	for (j=0; r>bInfo.weight; j++) bInfo=_scalers[j+1];
	return (_scalers[j].scaler);
}

double Predictor::computeBranchLength(const nodeptr node) const {

	double branchLength = 0.0;

	if (_pllPartitions->numberOfPartitions > 1) {
		int sumwgt = 0;

		/* compute the total weight */
		for (unsigned int i = 0;
						i < (unsigned int) _pllPartitions->numberOfPartitions; i++) {
			if (!isMissingBranch(node, i)) {
				sumwgt = _pllPartitions->partitionData[i]->width;
			}
		}

		if (sumwgt == 0) {
			cerr << "Error: There is no information available for stealing branch length ("
					<< node->number << "," << node->back->number << ")" << endl;
		}

		assert(sumwgt > 0);

		/* compute the branch length as the average over all other partitions */
		for (unsigned int i = 0;
				i < (unsigned int) _pllPartitions->numberOfPartitions; i++) {
			if (_partitionNumber != i && !isMissingBranch(node, i)) {
				double factor = (double) _pllPartitions->partitionData[i]->width / (double) sumwgt;
				branchLength += pllGetBranchLength(_pllTree, node, i) * factor;
			}
		}
	} else {
		/* return the current and only branch length */
		branchLength = pllGetBranchLength(_pllTree, node, _partitionNumber);
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

	return branchLength;
}

#define WMOD 0.5
#define WMEAN 0
#define WSIGMA 1

static double computeWeight(double x) {
	double w = 1/(WSIGMA*sqrt(2*3.1415926))*exp(-((x-WMEAN)*(x-WMEAN))/(2*WSIGMA*WSIGMA));
	return w;
}

bool Predictor::isAncestor(nodeptr node) const {
	if (find(_missingSubtreesAncestors.begin(),
				_missingSubtreesAncestors.end(), node)
				== _missingSubtreesAncestors.end())
		return false;
	return true;
}

bool Predictor::isMissingBranch(const nodeptr node, int partition) const {
	if (find(_missingBranches->at(partition).begin(),
			_missingBranches->at(partition).end(), node)
				== _missingBranches->at(partition).end())
		return false;
	return true;
}

void Predictor::getBranches(nodeptr node, int depth, double * cumWeight, vector<branchInfo> & branches) const {

	double localSum = 0.0, ratio = 0.0, curWeight;
	int sumWgt = 0;

	if (!(isAncestor(node->next) || isAncestor(node->next->next))) {

		for (unsigned int i = 0;
				i < (unsigned int) _pllPartitions->numberOfPartitions; i++) {
			if (!isMissingBranch(node, i)) {
				sumWgt += _pllPartitions->partitionData[i]->width;
			}
		}

		/* compute the current weighted branch length ratio */
		for (unsigned int i = 0;
				i < (unsigned int) _pllPartitions->numberOfPartitions; i++) {
			if (!isMissingBranch(node, i)) {
				double factor = (double) _pllPartitions->partitionData[i]->width
						/ (double) sumWgt;
				localSum += pllGetBranchLength(_pllTree, node, i) * factor;
			}
		}
		curWeight = computeWeight(depth * WMOD);
		ratio = pllGetBranchLength(_pllTree, node, _partitionNumber)/localSum;

		stringstream ss;
						ss << "(" << node->number << "," << node->back->number << ")";
		cout << setw(10) << right << ss.str() << " - weight: " << curWeight << " ratio: "
				<< ratio << endl;

		/* add the new sample */
		branchInfo bInfo;
		bInfo.branchNumber = node->number;
		bInfo.scaler = ratio;
		bInfo.weight = curWeight;
		branches.push_back(bInfo);

		if ((unsigned int) node->number > numberOfTaxa) {
			/* inspect children */
			double child1W = 0.0, child2W = 0.0;
			getBranches(node->next->back, depth+1, &child1W, branches);
			getBranches(node->next->next->back, depth+1, &child2W, branches);
			*cumWeight += child1W + child2W;
		}
		*cumWeight += curWeight;
	}
}

double Predictor::getSumBranches(nodeptr node, int depth, double * weight) const {

	double childrenSum = 0.0, localSum = 0.0, ratio = 0.0;

	/* we define 2 weights: one according to the distance to the rooting branch
	 * and another one accorting to the number of sites in the partition.
	 */
	double distanceWeight;
	int sumWgt = 0;

	if (!(isAncestor(node->next) || isAncestor(node->next->next))) {

		for (unsigned int i = 0;
				i < (unsigned int) _pllPartitions->numberOfPartitions; i++) {
			if (!isMissingBranch(node, i)) {
				sumWgt += _pllPartitions->partitionData[i]->width;
			}
		}

		/* compute the current weighted branch length ratio */
		for (unsigned int i = 0;
				i < (unsigned int) _pllPartitions->numberOfPartitions; i++) {
			if (!isMissingBranch(node, i)) {
				double factor = (double) _pllPartitions->partitionData[i]->width / (double) sumWgt;
				localSum += pllGetBranchLength(_pllTree, node, i) * factor;
			}
		}

		distanceWeight = computeWeight(depth * WMOD);
		ratio = pllGetBranchLength(_pllTree, node, _partitionNumber)/localSum;

		stringstream ss;
				ss << "(" << node->number << "," << node->back->number << ")";
		cout << setw(10) << right << ss.str() << setprecision(6) << " - weight: " << distanceWeight << " ratio: "
				<< ratio << endl;
		ratio *= distanceWeight;

		if ((unsigned int) node->number > numberOfTaxa) {
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
//		nodeptr startingNode = _missingSubtreesAncestors[0];
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

void Predictor::evolveNode(const nodeptr node, const double * ancestralProbabilities, const double * ancestralPMatrix) {

#ifdef PRINT_TRACE
	cout << "TRACE: Mutating node " << node->number << endl;
#endif
#ifdef PRINT_ANCESTRAL
	if (ancestralSequence) {
		cout << "INFO: Ancestral: ";
	    Utils::printSequence(ancestralSequence);
	}
#endif
	unsigned int n = _partitionLength;
	double branchLength = computeBranchLength(node);
	cout << "  - Estimated branch length for node " << setprecision(5) << node->number << " (partition " << _pllPartitions->partitionData[_partitionNumber]->partitionName << ") = " << branchLength << endl;

	if (node->number > _pllTree->mxtips) {
		double * currentPMatrix = (double *) malloc(
				(size_t) (numberOfRateCategories * _numberOfStates
						* _numberOfStates) * sizeof(double *));
		mutatePMatrix(currentPMatrix, ancestralPMatrix, branchLength);

		evolveNode(node->next->back, ancestralProbabilities, currentPMatrix);
		evolveNode(node->next->next->back, ancestralProbabilities, currentPMatrix);

		free(currentPMatrix);
	} else {
		/* compute the new probability matrix */
		/* construct the new marginal probabilities matrix */
	//	void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
	//	                 const int M, const int N, const int K,
	//	                 const double alpha, const double *A, const int lda
	//	                 const double *B, const int ldb,
	//	                 const double beta, double *C, const int ldc);
		// Multiply P' x M

		double * currentProbabilities = (double *) malloc((size_t) numberOfRateCategories * _numberOfStates * n * sizeof(double));
		for (int cat = 0; cat < (int)numberOfRateCategories; cat++) {
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
				n, _numberOfStates, _numberOfStates,
				1.0, ancestralProbabilities, _numberOfStates,
				&ancestralPMatrix[cat * _numberOfStates* _numberOfStates], _numberOfStates,
				0.0, &currentProbabilities[cat * _numberOfStates * n], _numberOfStates);

			if (USE_FIXED_ANCESTRAL) {
				double * M = &currentProbabilities[cat * _numberOfStates * n];
				for (unsigned int i = 0; i < n; i++) {
					double sum = 0.0;
					for (unsigned int j = 0; j < _numberOfStates; j++) {
						int nextIndex = _numberOfStates * i + j;
						sum += M[nextIndex];
					}
					assert(Utils::floatEquals(sum, 1.0));
				}
			} else {
				/* Make M matrix cummulative */
				double * M = &currentProbabilities[cat * _numberOfStates * n];
				for (unsigned int i = 0; i < n; i++) {
					for (unsigned int j = 1; j < _numberOfStates; j++) {
						int nextIndex = _numberOfStates * i + j;
						M[nextIndex] += M[nextIndex - 1];
					}
					assert(
							Utils::floatEquals(M[_numberOfStates * (i + 1) - 1],
									1.0));
				}
			}
		}

		/* draw random ancestral states */
		char * ancestralSequence = (char *) malloc(n + 1);
		switch (categoriesMode) {
		case CAT_AVERAGE: {
			cerr << "Not implemented yet" << endl;
			exit(EX_UNIMPLEMENTED);
			break;
		}
		case CAT_RANDOM:
		case CAT_ESTIMATE: {

			short * siteCatPtr = _catToSite;
			char * seqPtr = ancestralSequence;
			if (USE_FIXED_ANCESTRAL) {
				for (unsigned int i = 0; i < _partitionLength; i++) {
					char newState = _currentModel->getMostProbableState(
							&currentProbabilities[(*siteCatPtr)
									* _numberOfStates * _partitionLength
									+ i * _numberOfStates]);
					*seqPtr = newState;
					seqPtr++;
					siteCatPtr++;
				}
			} else {
				for (unsigned int i = 0; i < _partitionLength; i++) {
					char newState = _currentModel->getState(
							&currentProbabilities[(*siteCatPtr)
									* _numberOfStates * _partitionLength
									+ i * _numberOfStates]);
					*seqPtr = newState;
					seqPtr++;
					siteCatPtr++;
				}
			}
			break;
		}
		}
		ancestralSequence[n] = '\0';

		/* compute the new sequence from the probabilities */
		char * currentSequence = (char *) malloc (n + 1);
		mutateSequence(currentSequence, ancestralSequence, branchLength);
		free(ancestralSequence);

		/* set the new sequence */
		memcpy(&(_pllAlignment->sequenceData[seqIndexTranslate[node->number]][_start]), currentSequence, _partitionLength);
		/* remove visited taxon */
		_missingSequences.erase(remove(_missingSequences.begin(), _missingSequences.end(), node->number), _missingSequences.end());

		free(currentSequence);
		free(currentProbabilities);
	}
}

void Predictor::mutatePMatrix(double * currentPMatrix,
		const double * ancestralPMatrix, double branchLength) {

	double * gammaRates = (double *) alloca (numberOfRateCategories * sizeof(double));

	pllMakeGammaCats(_pllPartitions->partitionData[_partitionNumber]->alpha, gammaRates, numberOfRateCategories, false);

	/* construct and validate P matrix */
	double matrix[numberOfRateCategories][_numberOfStates*_numberOfStates];
	for (unsigned int i = 0; i < numberOfRateCategories; i++) {
		_currentModel->setMatrix(matrix[i], gammaRates[i] * branchLength, false);

		/* validation */
		for (unsigned int j = 0; j < _numberOfStates; j++) {
			double sum = 0.0;
			for (unsigned int k = 0; k < _numberOfStates; k++) {
				sum += matrix[i][j*_numberOfStates + k];
			}
			assert(Utils::floatEquals(sum, 1.0));
		}
	}

	if (ancestralPMatrix) {
		/* multiply P matrices */
		int n = _numberOfStates;
		for (int cat = 0; cat < (int) numberOfRateCategories; cat++) {
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0,
					matrix[cat], n, &ancestralPMatrix[cat * n * n], n, 0.0,
					&currentPMatrix[cat * n * n], n);
		}
	} else {
		/* directly set the P matrix */
		for (unsigned int i = 0; i < numberOfRateCategories; i++) {
			memcpy(&currentPMatrix[i * _numberOfStates * _numberOfStates],
					matrix[i],
					_numberOfStates * _numberOfStates * sizeof(double));
		}
	}
}

void Predictor::evolveNode(const nodeptr node, const char * ancestralSequence) {

#ifdef PRINT_TRACE
	cout << "TRACE: Mutating node " << node->number << endl;
#endif
#ifdef PRINT_ANCESTRAL
	if (ancestralSequence) {
		cout << "INFO: Ancestral: ";
	    Utils::printSequence(ancestralSequence);
	}
#endif

	double branchLength = computeBranchLength(node);
	cout << "  - Estimated branch length for node " << setprecision(5) << node->number << " (partition " << _pllPartitions->partitionData[_partitionNumber]->partitionName << ") = " << branchLength << endl;

	char * currentSequence = (char *) malloc(strlen(ancestralSequence) + 1);
	mutateSequence(currentSequence, ancestralSequence, branchLength);

	if (node->number > _pllTree->mxtips) {
#ifdef PRINT_TRACE
		cout << "TRACE: Recurse for Mutating node " << node->next->back->number
				<< " and " << node->next->next->back->number << endl;
#endif
		evolveNode(node->next->back, currentSequence);
		evolveNode(node->next->next->back, currentSequence);
	} else {
		/* set the new sequence */
		memcpy(&(_pllAlignment->sequenceData[seqIndexTranslate[node->number]][_start]), currentSequence, _partitionLength);
		/* remove visited taxon */
		_missingSequences.erase(remove(_missingSequences.begin(), _missingSequences.end(), node->number), _missingSequences.end());
	}
	free(currentSequence);
}

void Predictor::getRootingNodes() {
	assert ( !_missingSubtreesAncestors.size() );
	if ( _missingBranches->at(_partitionNumber).size()) {
		unsigned int i;
		for (i = 0;
				i < _missingBranches->at(_partitionNumber).size()
						&& _missingBranches->at(_partitionNumber)[i]->number
								<= _pllTree->mxtips; i++)
			;

		nodeptr currentNode = _missingBranches->at(_partitionNumber)[i++];
		int counter = 0;
		for ( ; i<_missingBranches->at(_partitionNumber).size(); i++) {
			int testNode = _missingBranches->at(_partitionNumber)[i]->number;
			assert ( testNode >=  currentNode->number );
			if ( testNode == currentNode->number ) {
				counter++;
			} else {
				/* each inner node should appear either 1 or 3 times */
				assert ( counter == 0 || counter == 2);
				if ( counter == 0 ) {
					_missingSubtreesAncestors.push_back(currentNode);
				} else {
					counter = 0;
				}
				currentNode = _missingBranches->at(_partitionNumber)[i];
			}
		}
	} else {
		return;
	}
}

void Predictor::predictMissingSequences( const pllAlignmentData * originalSequence ) {

#ifdef PRINT_TRACE
	for (size_t i = 0; i < _missingSequences.size(); i++) {
		cout << "TRACE: Missing sequence " << _missingSequences[i] << endl;
	}
#endif

	vector<int> missingSequencesCopy;
	if (originalSequence) {
		missingSequencesCopy = _missingSequences;
	}

	getRootingNodes();

	if (_missingSequences.size()) {
		if (branchLengthsMode == BL_SCALE) {
			/* compute branch length scaler */
			if (_missingSubtreesAncestors.size()) {
				/* traverse the tree starting at one arbitrary ancestor node */
				nodeptr startingNode = _missingSubtreesAncestors[0];

				cout << "Branch length scaler:" << endl;
				double cumWeight = 0.0, weight = 0.0;
				_branchLengthScaler = getSumBranches(startingNode->next->back, 0,
						&cumWeight);
				_branchLengthScaler += getSumBranches(
						startingNode->next->next->back, 0, &weight);
				cumWeight += weight;
				_branchLengthScaler /= cumWeight;
			}
			cout << "Applying branch length scaler " << _branchLengthScaler << endl << endl;
		} else if (branchLengthsMode == BL_DRAW) {
			/* traverse the tree starting at one arbitrary ancestor node */
			nodeptr startingNode = _missingSubtreesAncestors[0];

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
		nodeptr ancestor = _missingSubtreesAncestors[nextAncestor++];

		nodeptr startNode = ancestor;

	#ifdef PRINT_TRACE
		cout << "TRACE: Updating partials for " << startNode->number << endl;
	#endif
		pllUpdatePartialsAncestral(_pllTree, _pllPartitions, startNode, false);
	#ifdef PRINT_TRACE
		cout << "TRACE: Generating ancestral for " << startNode->number << endl;
	#endif

		/*
		 * Compute probs size according to the maximum number of states.
		 * This is necessary in case there are mixed protein-dna partitions.
		 */
		int probsSize = 0;
		for (int i=0; i<_pllPartitions->numberOfPartitions; i++) {
			probsSize = max(_pllPartitions->partitionData[i]->states, probsSize);
		}
		probsSize *= sequenceLength;

		char * ancestral = (char *) malloc( (size_t) sequenceLength + 1);
		double * probs = (double *) malloc(
				(size_t) probsSize * sizeof(double));
		pllGetAncestralState(_pllTree, _pllPartitions, startNode, probs, ancestral, false);

	#ifdef DEBUG
		printNodes(_pllTree, _pllPartitions);
		printBranchLengths(_pllTree, _pllPartitions);
	#endif

		switch (predictionMode) {
		case PRED_ANCSEQ:
		{
			/* discard MAP vector */
			free (probs);

			/* extract the ancestral for the partition */
			char * partAncestral = (char *) malloc(_partitionLength + 1);
			memcpy(partAncestral, &(ancestral[_start]), _partitionLength);
			partAncestral[_partitionLength] = '\0';
			free(ancestral);

			evolveNode(ancestor->back, partAncestral);

			free(partAncestral);
			break;
		}
		case PRED_MAP:
			/* discard ancestral sequence */
			free(ancestral);

			/* extract the ancestral for the partition */
			size_t partProbsSize = (probsSize/sequenceLength) * _partitionLength;
			double * probsAncestral = (double *) malloc(partProbsSize * sizeof(double));
			memcpy(probsAncestral, &(probs[_start * _numberOfStates]), partProbsSize * sizeof(double));
			free (probs);

			/* validate */
//			for (int i=0; i<_partitionLength; i++) {
//				for (int j=0; j<_numberOfStates; j++) {
//					cout << probsAncestral[i*_numberOfStates + j] << " ";
//				}
//				cout << endl;
//			}

			evolveNode(ancestor->back, probsAncestral);
			free (probsAncestral);

			break;
		}

		if (originalSequence) {
			/* compute similarity */
			_seqSimilarity = 0.0;
			for (size_t i=0; i<missingSequencesCopy.size(); i++) {
				int seq = missingSequencesCopy[i];
				double simCount = 0.0;
				unsigned int seqLen = 0;
				if (Utils::getDataType(_pllPartitions, _partitionNumber)
						== DT_NUCLEIC) {
					for (size_t j = _start; j < _end; j++) {
						bool validForComp;
						simCount += Utils::compareNucStates(_pllAlignment->sequenceData[seqIndexTranslate[seq]][j],
								originalSequence->sequenceData[seqIndexTranslate[seq]][j], &validForComp);
						seqLen += validForComp;
					}
					cout << "Comparing sequence " << taxaNames[seq] << ": " << simCount << " valid sites." << endl;
					_seqSimilarity += simCount/seqLen/missingSequencesCopy.size();
				}
				cout << "Similarity in sequence " << seq <<"/" << _pllPartitions->partitionData[_partitionNumber]->partitionName << ": " << 100*simCount/seqLen << "%" << endl;
			}
		}
	}

	/* ensure there is no subtree left */
	assert(nextAncestor == _missingSubtreesAncestors.size());
}

void Predictor::predictAllSequences( void ) {

	_missingSequences.resize(1);
	for (unsigned int taxon=1; taxon <= numberOfTaxa; taxon++) {
		_missingSequences[0] = taxon;
		nodeptr ancestor = _pllTree->nodep[taxon]->back;
		/*
		 * Compute probs size according to the maximum number of states.
		 * This is necessary in case there are mixed protein-dna partitions.
		 */
		int probsSize = 0;
		for (int i = 0; i < _pllPartitions->numberOfPartitions; i++) {
			probsSize = max(_pllPartitions->partitionData[i]->states,
					probsSize);
		}
		probsSize *= sequenceLength;
		char * ancestral = (char *) malloc( (size_t) sequenceLength + 1);
		double * probs = (double *) malloc(
				(size_t) probsSize * sizeof(double));
		pllGetAncestralState(_pllTree, _pllPartitions, ancestor, probs, ancestral, false);
		free (probs);

		/* extract the ancestral for the partition */
		char * partAncestral = (char *) malloc(_partitionLength + 1);
		memcpy(partAncestral, &(ancestral[_start]), _partitionLength);
		partAncestral[_partitionLength] = '\0';
		free(ancestral);

		evolveNode(ancestor->back, partAncestral);
		free(partAncestral);
	}
}

} /* namespace seqpred */
