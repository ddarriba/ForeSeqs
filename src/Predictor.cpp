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
#include <algorithm>
#include <cstring>
#include <cassert>

using namespace std;

namespace seqpred {

Predictor::Predictor(pllInstance * tree, partitionList * partitions,
		pllAlignmentData * phylip, int partitionNumber) :
		_pllTree(tree), _pllPartitions(partitions), _pllAlignment(phylip), _partitionNumber(
				partitionNumber) {

	if (Utils::getDataType(partitions, partitionNumber) == DT_NUCLEIC) {
		_currentModel = new DnaModel(partitions, partitionNumber);
		_numberOfStates = 4;
	} else {
		_currentModel = new ProteinModel(partitions, partitionNumber);
		_numberOfStates = 20;
	}

	/* get information from the partition */
	_start = partitions->partitionData[partitionNumber]->lower;
	_end = partitions->partitionData[partitionNumber]->upper;
	_partitionLength = _end - _start;

	/* get taxa with missing data in the partition */
	_missingSequences = findMissingSequences();

}

Predictor::~Predictor( void ) {
	delete _currentModel;
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

	double gammaRates[numberOfRateCategories];
	pllMakeGammaCats(_pllPartitions->partitionData[_partitionNumber]->alpha,gammaRates, numberOfRateCategories, false);

	/* random assignment of sites to categories */
	short categories[_partitionLength];
	for (unsigned int i=0; i<_partitionLength; i++) {
		categories[i]=(int)(numberOfRateCategories * Utils::genRand());
	}

	short * siteCatPtr = categories;
	char * seqPtr = currentSequence;

	/* construct P matrix */
	double matrix[numberOfRateCategories][_numberOfStates*_numberOfStates];
	for (int i = 0; i < numberOfRateCategories; i++) {
		_currentModel->setMatrix(matrix[i], gammaRates[i] * branchLength);
	}
	for (unsigned int i = 0; i < _partitionLength; i++) {
		*seqPtr = _currentModel->getState(matrix[*siteCatPtr]+(_currentModel->getStateIndex(*seqPtr) * _numberOfStates));
		seqPtr++;
		siteCatPtr++;
	}
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
	if (_pllPartitions->numberOfPartitions > 1) {
		/* compute the branch length as the average over all other partitions */
		for (int i = 0; i < _pllPartitions->numberOfPartitions; i++) {
			if (_partitionNumber != i) {
				branchLength += pllGetBranchLength(_pllTree, node, i);
			}
		}
		branchLength /= (_pllPartitions->numberOfPartitions - 1);
	} else {
		/* return the current and only branch length */
		branchLength = pllGetBranchLength(_pllTree, node,
				_partitionNumber);
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
	cout << "  - Estimated branch length for node " << node->number << " (partition " << _pllPartitions->partitionData[_partitionNumber]->partitionName << ") = " << branchLength << endl;

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

void Predictor::predictMissingSequences( void ) {

#ifdef PRINT_TRACE
	for (size_t i = 0; i < _missingSequences.size(); i++) {
		cout << "TRACE: Missing sequence " << _missingSequences[i] << endl;
	}
#endif

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
	}
}

} /* namespace seqpred */
