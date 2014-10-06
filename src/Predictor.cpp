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
#include <cassert>
#include <cmath>

using namespace std;

namespace seqpred {

Predictor::Predictor(pllInstance * tree, partitionList * partitions,
		int partitionNumber) :
		tree(tree), partitions(partitions), partitionNumber(partitionNumber),
		curModel(partitions, partitionNumber) {

	/* get information from the partition */
	start = partitions->partitionData[partitionNumber]->lower;
	end = partitions->partitionData[partitionNumber]->upper;
	length = end - start;
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

Predictor::~Predictor() {
	for (size_t i = 0; i < missingSequences.size(); i++) {
		free(predictedSequences[missingSequences[i]]);
	}
}

vector<int> Predictor::findMissingSequences() {
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
boolean Predictor::subtreeIsMissing(nodeptr node) {
	if (node->number > tree->mxtips) {
		return (subtreeIsMissing(node->next->back)
				& subtreeIsMissing(node->next->next->back));
	} else {
		return (find(missingSequences.begin(), missingSequences.end(),
				node->number) != missingSequences.end());
	}
}

nodeptr Predictor::findMissingDataAncestor() {

	nodeptr currentNode;

	if (!missingSequences.size()) {
		cerr << "ERROR: No missing sequences" << endl;
		return 0;
	}

	/* start searching in a random missing node */
	currentNode = tree->nodep[missingSequences[0]]->back;

	while (true) {
		bool missingRight = subtreeIsMissing(currentNode->next->back);
		bool missingLeft = subtreeIsMissing(currentNode->next->next->back);

		if (missingRight && missingLeft) {
			cerr << "ERROR: Everything is missing!!" << endl;
			return 0;
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

char Predictor::getState(double * P) {
	int j;
	double r = Utils::genRand();
	for (j=0; r>(*P) && j<numStates-1; j++) P++;
	return (states[pow(2,j)]);
}

void Predictor::mutateSequence(char * currentSequence,
		char * ancestralSequence, double branchLength) {

	if ( length != strlen(ancestralSequence) ) {
		cerr << "ERROR: Length of ancestral sequence (" << strlen(ancestralSequence)
				<< ") differ from the expected length (" << length << ")" << endl;
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
	short categories[length];
	for (unsigned int i=0; i<length; i++) {
		categories[i]=(int)(numRateCategories * Utils::genRand());
	}

	short * R = categories;
	char * P = currentSequence;

	double *matrix[4];
	for (int i=0; i<4; i++) {
		matrix[i] = (double*)malloc(16 * sizeof(double));
	}
	for (int i = 0; i < numRateCategories; i++) {
		curModel.setMatrix(matrix[i], gammaRates[i] * branchLength);
	}
	for (unsigned int i = 0; i < length; i++) {
		*P = getState(matrix[*R]+((statesMap[*P]) * numStates));
		P++;
		R++;
	}
	for (int i = 0; i < 4; i++) {
		free(matrix[i]);
	}

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

	double branchLength = pllGetBranchLength(tree, node, partitionNumber);
	char * currentSequence = (char *) malloc(strlen(ancestralSequence) + 1);
	mutateSequence(currentSequence, ancestralSequence, branchLength);

	if (node->number > tree->mxtips) {
#ifdef PRINT_TRACE
		cout << "TRACE: Recurse for Mutating node " << node->next->back->number
				<< " and " << node->next->next->back->number << endl;
#endif
		evolveNode(node->next->back, currentSequence);
		evolveNode(node->next->next->back, currentSequence);
		free(currentSequence);
	} else {
		predictedSequences[node->number] = currentSequence;
	}
}

void Predictor::predictMissingSequences() {

#ifdef PRINT_TRACE
	for (size_t i = 0; i < missingSequences.size(); i++) {
		cout << "TRACE: Missing sequence " << missingSequences[i] << endl;
	}
#endif

	nodeptr ancestor = findMissingDataAncestor();

	nodeptr startNode = ancestor; //tree->nodep[5]->back->next->next->back;

#ifdef PRINT_TRACE
	cout << "TRACE: Updating partials for " << startNode->number << endl;
#endif
	pllUpdatePartialsAncestral(tree, partitions, startNode);
#ifdef PRINT_TRACE
	cout << "TRACE: Generating ancestral for " << startNode->number << endl;
#endif

	char * ancestral = (char *) malloc(length + 1);
	double * probs = (double *) malloc(
			length * numStates * sizeof(double));
	pllGetAncestralState(tree, partitions, startNode, probs, ancestral);

	free (probs);

	evolveNode(ancestor->back, ancestral);
	free (ancestral);
}

} /* namespace seqpred */
