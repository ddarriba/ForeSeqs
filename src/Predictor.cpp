/*
 * Predictor.cpp
 *
 *  Created on: Oct 1, 2014
 *      Author: diego
 */

#include "Predictor.h"
#include "Utils.h"

#include <iostream>
#include <algorithm>
#include <cstring>
#include "DnaModel.h"

#include <cassert>
#include <cmath>

using namespace std;

namespace seqpred {

Predictor::Predictor(pllInstance * tree, partitionList * partitions,
		pllAlignmentData * phylip, int partitionNumber) :
		tree(tree), partitions(partitions), phylip(phylip), partitionNumber(
				partitionNumber) {

	if (dataType == DT_NUCLEIC) {
		curModel = new DnaModel(partitions, partitionNumber);
	} else {
		cerr << "Unimplemented Data Type" << endl;
		exit(EX_IOERR);
	}

	/* get information from the partition */
	start = partitions->partitionData[partitionNumber]->lower;
	end = partitions->partitionData[partitionNumber]->upper;
	partitionLength = end - start;
	numStates = partitions->partitionData[partitionNumber]->states;

	/* number of gamma rate categories fixed to 4 */
	numRateCategories = 4;

	/* get taxa with missing data in the partition */
	missingSequences = findMissingSequences();

	/* create the different states */
	/* TODO: So far we set it to DNA states */
	if (numStates != 4) {
		cerr << " Not implemented yet for " << numStates << " states" << endl;
		exit(1);
	}
}

Predictor::~Predictor( void ) {
	delete curModel;
}

vector<int> Predictor::findMissingSequences( void ) const {
	vector<int> missingSeqs;

	int missing;
	for (int i = 1; i <= tree->ntips; i++) {
		missing = 1;
		for (unsigned int j = start; j < end; j++) {
			missing &= (tree->yVector[i][j] == PLL_UNDEFINED_SITE);
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
 */
boolean Predictor::subtreeIsMissing(nodeptr node) const {
	if (node->number > tree->mxtips) {
		return (subtreeIsMissing(node->next->back)
				& subtreeIsMissing(node->next->next->back));
	} else {
		return (find(missingSequences.begin(), missingSequences.end(),
				node->number) != missingSequences.end());
	}
}

nodeptr Predictor::findMissingDataAncestor( void ) const {

	nodeptr currentNode;

	if (!missingSequences.size()) {
		cerr << "ERROR: No missing sequences" << endl;
		return(0);
	}

	/* start searching in a random missing node */
	currentNode = tree->nodep[missingSequences[0]]->back;

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

char Predictor::getState(double * P) const {
	int j;
	double r = Utils::genRand();
	for (j=0; r>(*P) && j<numStates-1; j++) P++;
	return (states[pow(2,j)]);
}

void Predictor::mutateSequence(char * currentSequence,
		char * ancestralSequence, double branchLength) {

	if ( partitionLength != strlen(ancestralSequence) ) {
		cerr << "ERROR: Length of ancestral sequence (" << strlen(ancestralSequence)
				<< ") differ from the expected length (" << partitionLength << ")" << endl;
		assert(0);
	}

#ifdef PRINT_TRACE
	cout << "TRACE: Simulating sequence..." << endl;
#endif

	/* start mutating the ancestral sequence */
	strcpy(currentSequence, ancestralSequence);

	double gammaRates[numRateCategories];
	pllMakeGammaCats(partitions->partitionData[partitionNumber]->alpha,gammaRates, numRateCategories, false);

	/* random assignment of sites to categories */
	short categories[partitionLength];
	for (unsigned int i=0; i<partitionLength; i++) {
		categories[i]=(int)(numRateCategories * Utils::genRand());
	}

	short * siteCatPtr = categories;
	char * seqPtr = currentSequence;

	/* construct P matrix */
	double matrix[4][16];
	for (int i = 0; i < numRateCategories; i++) {
		curModel->setMatrix(matrix[i], gammaRates[i] * branchLength);
	}
	for (unsigned int i = 0; i < partitionLength; i++) {
		*seqPtr = getState(matrix[*siteCatPtr]+((statesMap[*seqPtr]) * numStates));
		seqPtr++;
		siteCatPtr++;
	}
}

#ifdef DEBUG
void printNodes(pllInstance * tree, partitionList * partitions) {
	/* compute the branch length as the average over all other partitions */
	for (int node = 1; node <= 2 * tree->mxtips - 3; node++) {
		cout << node << " " << tree->nodep[node]->next->back->number << " "
				<< tree->nodep[node]->next->next->back->number << endl;
	}
}

void printBranchLengths(pllInstance * tree, partitionList * partitions) {
	/* compute the branch length as the average over all other partitions */
	for (int i=0; i<partitions->numberOfPartitions; i++) {
		cout << i << " -> ";
		for (int node=1; node <= 2*tree->mxtips-3; node++) {
			cout << " " << pllGetBranchLength(tree, tree->nodep[node], i);
		}
		cout << endl;
	}
}
#endif

double Predictor::computeBranchLength(nodeptr node) const {
	double branchLength = 0.0;
	if (partitions->numberOfPartitions > 1) {
		/* compute the branch length as the average over all other partitions */
		for (int i = 0; i < partitions->numberOfPartitions; i++) {
			if (partitionNumber != i) {
				branchLength += pllGetBranchLength(tree, node, i);
			}
		}
		branchLength /= (partitions->numberOfPartitions - 1);
	} else {
		/* return the current and only branch length */
		branchLength = pllGetBranchLength(tree, node,
				partitionNumber);
	}
	return branchLength;
}

void Predictor::evolveNode(nodeptr node, char * ancestralSequence) {

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
	cout << "  - Estimated branch length for node " << node->number << " (partition " << partitions->partitionData[partitionNumber]->partitionName << ") = " << branchLength << endl;

	char * currentSequence = (char *) malloc(strlen(ancestralSequence) + 1);
	mutateSequence(currentSequence, ancestralSequence, branchLength);

	if (node->number > tree->mxtips) {
#ifdef PRINT_TRACE
		cout << "TRACE: Recurse for Mutating node " << node->next->back->number
				<< " and " << node->next->next->back->number << endl;
#endif
		evolveNode(node->next->back, currentSequence);
		evolveNode(node->next->next->back, currentSequence);
	} else {
		/* set the new sequence */
		memcpy(&(phylip->sequenceData[node->number][start]), currentSequence, partitionLength);
		/* remove visited taxon */
		missingSequences.erase(remove(missingSequences.begin(), missingSequences.end(), node->number), missingSequences.end());
	}
	free(currentSequence);
}

void Predictor::predictMissingSequences( void ) {

#ifdef PRINT_TRACE
	for (size_t i = 0; i < missingSequences.size(); i++) {
		cout << "TRACE: Missing sequence " << missingSequences[i] << endl;
	}
#endif

	/* loop over all possible subtrees with missing data */
	while (missingSequences.size()) {
		cout << "Predicting subtree" << endl;
		nodeptr ancestor = findMissingDataAncestor();

		nodeptr startNode = ancestor;

	#ifdef PRINT_TRACE
		cout << "TRACE: Updating partials for " << startNode->number << endl;
	#endif
		pllUpdatePartialsAncestral(tree, partitions, startNode);
	#ifdef PRINT_TRACE
		cout << "TRACE: Generating ancestral for " << startNode->number << endl;
	#endif

		char * ancestral = (char *) malloc(sequenceLength + 1);
		double * probs = (double *) malloc(
				sequenceLength * numStates * sizeof(double));
		pllGetAncestralState(tree, partitions, startNode, probs, ancestral);
		free (probs);

		/* extract the ancestral for the partition */
		char * partAncestral = (char *) malloc(partitionLength + 1);
		memcpy(partAncestral, &(ancestral[start]), partitionLength);
		partAncestral[partitionLength] = '\0';
		free (ancestral);

	#ifdef DEBUG
		printNodes(tree, partitions);
		printBranchLengths(tree, partitions);
	#endif

		evolveNode(ancestor->back, partAncestral);
		free(partAncestral);
	}
}

} /* namespace seqpred */
