/*
 * Predictor.cpp
 *
 *  Created on: Oct 1, 2014
 *      Author: diego
 */

#include "Predictor.h"
#include "Utils.h"
#include "Model.h"

#include <iostream>
#include <algorithm>
#include <cstring>
#include <cassert>

#define PRINT_ANCESTRAL 1

using namespace std;

namespace seqpred {

Predictor::Predictor(pllInstance * tree, partitionList * partitions,
		int partitionNumber) :
		tree(tree), partitions(partitions), partitionNumber(partitionNumber) {

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
	states.resize(numStates);
	states[0] = 'A'; //1; //A
	states[1] = 'C'; //2; //C
	states[2] = 'G'; //4; //G
	states[3] = 'T'; //8; //T
	statesMap['a'] = 0;
	statesMap['A'] = 0;
	statesMap['c'] = 1;
	statesMap['C'] = 1;
	statesMap['g'] = 2;
	statesMap['G'] = 2;
	statesMap['t'] = 3;
	statesMap['T'] = 3;
}

Predictor::~Predictor() {
	/* nothing to do */
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
			cout << "INFO: Found ancestor in " << currentNode->number << endl;
			return currentNode;
		} else {
			/* move to next position */
			currentNode =
					(!missingRight) ?
							currentNode->next->back :
							currentNode->next->next->back;
			cout << "INFO: Moving node to " << currentNode->number << endl;
		}
	}
	return 0;
}

char Predictor::getState(double * P) {
	int j;
	double r = Utils::genRand();
	for (j=0; r>(*P) && j<numStates-1; j++) P++;
	return (states[j]);
}

void Predictor::mutateSequence(char * currentSequence,
		char * ancestralSequence, double branchLength) {

	if ( length != strlen(ancestralSequence) ) {
		cerr << "ERROR: Length of ancestral sequence (" << strlen(ancestralSequence)
				<< ") differ from the expected length (" << length << ")" << endl;
		assert(0);
	}

	cout << "   Simulating sequence..." << endl;

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
	Model curModel(partitions, partitionNumber);
	for (int i = 0; i < numRateCategories; i++) {
		curModel.setMatrix(matrix[i], gammaRates[i] * branchLength);
	}
	for (unsigned int i = 0; i < length; i++) {
		*P = getState(matrix[*R]+((statesMap[*P]) * numStates));
		P++;
		R++;
	}


}

void Predictor::evolveNode(nodeptr node, char * ancestralSequence) {

	cout << "INFO: Mutating node " << node->number << endl;
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
		cout << "   Recurse for Mutating node " << node->next->back->number
				<< " and " << node->next->next->back->number << endl;
		evolveNode(node->next->back, currentSequence);
		evolveNode(node->next->next->back, currentSequence);
	} else {
		cout << "Finished sequence for taxon " << node->number << endl;
		Utils::printSequence(currentSequence);
	}
	free(currentSequence);
}

void Predictor::predictMissingSequences() {

	for (size_t i = 0; i < missingSequences.size(); i++) {
		cout << "Missing sequence " << missingSequences[i] << endl;
	}

	nodeptr ancestor = findMissingDataAncestor();

	nodeptr startNode = ancestor; //tree->nodep[5]->back->next->next->back;

	cout << "Updating partials for " << startNode->number << endl;
	pllUpdatePartialsAncestral(tree, partitions, startNode);
	cout << "Generating ancestral for " << startNode->number << endl;

	char * ancestral = (char *) malloc(length + 1);
	double * probs = (double *) malloc(
			length * numStates * sizeof(double));
	pllGetAncestralState(tree, partitions, startNode, probs, ancestral);
//	for(unsigned int i=0; i<length; i++) {
//		switch (ancestral[i]) {
//		case 'A':
//		case 'a':
//			ancestral[i] = states[0];
//			break;
//		case 'C':
//		case 'c':
//			ancestral[i] = states[1];
//			break;
//		case 'G':
//		case 'g':
//			ancestral[i] = states[2];
//			break;
//		case 'T':
//		case 't':
//			ancestral[i] = states[3];
//			break;
//		default:
//			assert(0);
//		}
//	}
	cout << "STARTING SEQ: "; Utils::printSequence(ancestral);
	free (probs);

	evolveNode(ancestor->back, ancestral);
	free (ancestral);

}

} /* namespace seqpred */
