/*
 * Predictor.cpp
 *
 *  Created on: Oct 1, 2014
 *      Author: diego
 */

#include "Predictor.h"
#include <iostream>
#include <algorithm>
#include <cstring>

#define PRINT_ANCESTRAL 1

using namespace std;

namespace seqpred {

Predictor::Predictor(pllInstance * tree, partitionList * partitions,
		int partitionNumber) :
		tree(tree), partitions(partitions) {

	/* get information from the partition */
	start = partitions->partitionData[partitionNumber]->lower;
	end = partitions->partitionData[partitionNumber]->upper;
	length = end - start + 1;
	numStates = partitions->partitionData[partitionNumber]->states;

	/* number of gamma rate categories fixed to 4 */
	numRateCategories = 4;

	/* get taxa with missing data in the partition */
	missingSequences = findMissingSequences();
}

Predictor::~Predictor() {
	// TODO Auto-generated destructor stub
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

void Predictor::mutateSequence(char * currentSequence,
		char * ancestralSequence) {

	cout << "   Simulating sequence..." << endl;
	strcpy(currentSequence, ancestralSequence);
}

void Predictor::evolveNode(nodeptr node, char * ancestralSequence) {

	cout << "INFO: Mutating node " << node->number << endl;
#ifdef PRINT_ANCESTRAL
	if (ancestralSequence)
		cout << "INFO: Ancestral: " << ancestralSequence << endl;
#endif

	char * currentSequence = (char *) malloc(strlen(ancestralSequence) + 1);
	mutateSequence(currentSequence, ancestralSequence);

	if (node->number > tree->mxtips) {
		cout << "   Recurse for Mutating node " << node->next->back->number
				<< " and " << node->next->next->back->number << endl;
		evolveNode(node->next->back, currentSequence);
		evolveNode(node->next->next->back, currentSequence);
	} else {
		cout << "Finished sequence for taxon " << node->number << endl;
		cout << currentSequence << endl;
	}
	free(currentSequence);
}

void Predictor::predictMissingSequences() {

	for (size_t i = 0; i < missingSequences.size(); i++) {
		cout << "Missing sequence " << missingSequences[i] << endl;
	}

	nodeptr ancestor = findMissingDataAncestor();

	nodeptr startNode = ancestor; //tree->nodep[5]->back->next->next->back;
	pllEvaluateLikelihood(tree, partitions, tree->start, true, false);
	pllTreeToNewick(tree->tree_string, tree, partitions, tree->start->back,
			true, true, false, false, false, PLL_SUMMARIZE_LH, false, false);
	cout << "Tree: " << tree->tree_string << endl;
	cout << "Updating partials for " << startNode->number << endl;
	pllUpdatePartialsAncestral(tree, partitions, startNode);
	cout << "Generating ancestral for " << startNode->number << endl;
	char * ancestral = (char *) malloc(length + 1);
	double * probs = (double *) malloc(
			length * numStates * sizeof(double));
	pllGetAncestralState(tree, partitions, startNode, probs, ancestral);
	free (probs);

	evolveNode(ancestor->back, ancestral);
	free (ancestral);

}

} /* namespace seqpred */
