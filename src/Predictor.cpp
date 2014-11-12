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
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <cassert>

using namespace std;

namespace seqpred {

Predictor::Predictor(pllInstance * tree, partitionList * partitions,
		pllAlignmentData * phylip, int partitionNumber) :
		_pllTree(tree), _pllPartitions(partitions), _pllAlignment(phylip),
				_partitionNumber(partitionNumber), _numberOfStates(0),
				_start(partitions->partitionData[partitionNumber]->lower),
				_end(partitions->partitionData[partitionNumber]->upper),
				_partitionLength(_end - _start), _catToSite(0),
				_missingSequences(), _missingPartsCount(0),
				_currentModel(0), _seqSimilarity(0) {

	if (Utils::getDataType(partitions, partitionNumber) == DT_NUCLEIC) {
		_currentModel = new DnaModel(partitions, partitionNumber);
		_numberOfStates = 4;
	} else {
		_currentModel = new ProteinModel(partitions, partitionNumber);
		_numberOfStates = 20;
	}

	/* get taxa with missing data in the partition */
	_missingSequences = findMissingSequences();
	_missingPartsCount = _missingSequences.size();

	if (_missingSequences.size()) {
		if (categoriesMode != CAT_AVERAGE) {
			_catToSite = (short *) calloc((size_t) _partitionLength, sizeof(short));
			unsigned int catToSiteCount[numberOfRateCategories];
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
				double gammaRates[numberOfRateCategories];
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

//				for (unsigned int i=0; i<_partitionLength; i++) {
//					cout << _catToSite[i];
//				}
//				cout << endl;
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
						_missingSequences(other._missingSequences),
						_missingPartsCount(other._missingPartsCount),
						_currentModel(other._currentModel),	_seqSimilarity(other._seqSimilarity) {
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

vector<int> Predictor::findMissingSequences( void ) const {
	vector<int> missingSeqs;

	unsigned char undefinedSite = (_numberOfStates==4)?15:22;
	int missing;
	for (int i = 1; i <= _pllTree->ntips; i++) {
		missing = 1;
		for (unsigned int j = _start; j < _end; j++) {
			missing &= (_pllTree->yVector[i][j] == undefinedSite);
		}
		if (missing) {
			missingSeqs.push_back(i);
		}
	}
	return missingSeqs;
}

/*
 * Search if all the taxa in a subtree have missing sequences.
 * TODO: This algorithm can be improved using a hashtable for storing the already evaluated nodes.
 * However, the time taken for this part is negligible compared to the rest.
 */
boolean Predictor::subtreeIsMissing(const nodeptr node) const {
	if (node->number > _pllTree->mxtips) {
		return (subtreeIsMissing(node->next->back)
				& subtreeIsMissing(node->next->next->back));
	} else {
		return (find(_missingSequences.begin(), _missingSequences.end(),
				node->number) != _missingSequences.end());
	}
}

nodeptr Predictor::findMissingDataAncestor( void ) const {

	nodeptr currentNode;

	if (!_missingSequences.size()) {
		cerr << "ERROR: No missing sequences" << endl;
		return(0);
	}

	/* start searching in a random missing node */
	currentNode = _pllTree->nodep[_missingSequences[0]]->back;

	while (true) {
		bool missingRight = subtreeIsMissing(currentNode->next->back);
		bool missingLeft = subtreeIsMissing(currentNode->next->next->back);

		if (missingRight && missingLeft) {
			cerr << "ERROR: Everything is missing!!" << endl;
			exit(EX_IOERR);
		} else if (!(missingRight || missingLeft)) {
#ifdef PRINT_TRACE
			cout << "TRACE: Found ancestor in " << currentNode->number << endl;
#endif
			return currentNode;
		} else {
			/* move to next position */
			currentNode =
					(!missingRight) ?
							currentNode->next->back :
							currentNode->next->next->back;
#ifdef PRINT_TRACE
			cout << "TRACE: Moving node to " << currentNode->number << endl;
#endif
		}
	}
	return 0;
}

void Predictor::mutateSequence(char * currentSequence,
		const char * ancestralSequence, double branchLength) {

	char * seqPtr;
	double gammaRates[numberOfRateCategories];
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

	pllMakeGammaCats(_pllPartitions->partitionData[_partitionNumber]->alpha,gammaRates, numberOfRateCategories, false);

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
		double averageMatrix[_numberOfStates * _numberOfStates];
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

double Predictor::computeBranchLength(const nodeptr node) const {

	double branchLength = 0.0;

	switch (branchLengthsMode) {
	case BL_AVERAGE: {
		if (_pllPartitions->numberOfPartitions > 1) {
			/* compute the branch length as the average over all other partitions */
			for (unsigned int i = 0;
					i < (unsigned int) _pllPartitions->numberOfPartitions;
					i++) {
				if (_partitionNumber != i) {
					branchLength += pllGetBranchLength(_pllTree, node, i);
				}
			}
			branchLength /= (_pllPartitions->numberOfPartitions - 1);
		} else {
			/* return the current and only branch length */
			branchLength = pllGetBranchLength(_pllTree, node, _partitionNumber);
		}
		break;
	}
	case BL_DRAW: {
		cerr << "I AM SORRY: Unimplemented branch length stealing mode" << endl;
		exit(EX_UNIMPLEMENTED);
		break;
	}
	case BL_SCALE: {
		cerr << "I AM SORRY: Unimplemented branch length stealing mode" << endl;
		exit(EX_UNIMPLEMENTED);
		break;
	}
	}
	return branchLength;
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
		memcpy(&(_pllAlignment->sequenceData[node->number][_start]), currentSequence, _partitionLength);
		/* remove visited taxon */
		_missingSequences.erase(remove(_missingSequences.begin(), _missingSequences.end(), node->number), _missingSequences.end());
	}
	free(currentSequence);
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

	/* loop over all possible subtrees with missing data */
	while (_missingSequences.size()) {
		cout << "Predicting subtree" << endl;
		nodeptr ancestor = findMissingDataAncestor();

		nodeptr startNode = ancestor;

	#ifdef PRINT_TRACE
		cout << "TRACE: Updating partials for " << startNode->number << endl;
	#endif
		pllUpdatePartialsAncestral(_pllTree, _pllPartitions, startNode);
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
		pllGetAncestralState(_pllTree, _pllPartitions, startNode, probs, ancestral);
		free (probs);

		/* extract the ancestral for the partition */
		char * partAncestral = (char *) malloc(_partitionLength + 1);
		memcpy(partAncestral, &(ancestral[_start]), _partitionLength);
		partAncestral[_partitionLength] = '\0';
		free (ancestral);

	#ifdef DEBUG
		printNodes(_pllTree, _pllPartitions);
		printBranchLengths(_pllTree, _pllPartitions);
	#endif

		evolveNode(ancestor->back, partAncestral);
		free(partAncestral);

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
						simCount += Utils::compareNucStates(_pllAlignment->sequenceData[seq][j],
								originalSequence->sequenceData[seq][j], &validForComp);
						seqLen += validForComp;
					}
					_seqSimilarity += simCount/seqLen/missingSequencesCopy.size();
				}
				cout << "Similarity in sequence " << seq <<"/" << _pllPartitions->partitionData[_partitionNumber]->partitionName << ": " << 100*simCount/seqLen << "%" << endl;
			}
		}
	}
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
		pllGetAncestralState(_pllTree, _pllPartitions, ancestor, probs, ancestral);
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
