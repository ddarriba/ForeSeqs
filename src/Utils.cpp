/*
 * Utils.cpp
 *
 *  Created on: Oct 2, 2014
 *      Author: diego
 */

#include "Utils.h"

#include <algorithm>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace std;

namespace seqpred {

BLMode branchLengthsMode = BL_AVERAGE;
CatMode categoriesMode = CAT_ESTIMATE;
PredMode predictionMode = PRED_ANCSEQ;

/* Number of categories hardcoded to 4 */
unsigned int numberOfRateCategories = 4;
unsigned int numberOfTaxa, sequenceLength;
char ** taxaNames;
unsigned int * seqIndexTranslate;

double Utils::genRand(void) {
	return rand() * (1.0 / INT_MAX);
}

void Utils::printSequence(const char * sequence) {
	for (unsigned int i = 0; i < strlen(sequence); i++) {
		cout << sequence[i];
	}
	cout << endl;
}

void Utils::printVector(const double * vec, int len) {
	cout << "(";
	for (int i = 0; i < len; i++) {
		cout << vec[i] << ",";
	}
	cout << ")" << endl;
}

bool Utils::floatEquals(double v1, double v2) {
	return (std::abs(v1 - v2) <= EPSILON);
}

DataType Utils::getDataType(const partitionList * pllPartitions, int numberOfPartition) {
	switch (pllPartitions->partitionData[numberOfPartition]->states) {
	case 4:
		return DT_NUCLEIC;
	case 20:
		return DT_PROTEIC;
	default:
		cerr << "[ERROR] Invalid number of states (" <<
		pllPartitions->partitionData[numberOfPartition]->states
		<< ") for partition " << numberOfPartition << endl;
		exit(EX_IOERR);
	}
}

double Utils::compareNucStates(unsigned char state0, unsigned char state1, bool * validForComp) {

	*validForComp = state1=='A'||state1=='C'||state1=='G'||state1=='T';
	if(state0 == state1) {
		return 1.0;
	} else {
		/* TODO: Check partial similarity with ambiguities */
		return 0.0;
	}
}

vector< vector<int> > Utils::findMissingSequences( pllInstance * pllTree, partitionList * pllPartitions ) {

	vector< vector<int> > missingSeqs(pllPartitions->numberOfPartitions);

	for (int part = 0; part < pllPartitions->numberOfPartitions; part++) {

		unsigned char undefinedSite = (pllPartitions->partitionData[part]->states==4)?15:22;

		int start = pllPartitions->partitionData[part]->lower,
			  end = pllPartitions->partitionData[part]->upper;

		int missing;
		for (int i = 1; i <= pllTree->mxtips; i++) {
			missing = 1;
			for (int j = start; j < end; j++) {
				missing &= (pllTree->yVector[i][j] == undefinedSite);
			}
			if (missing) {
				missingSeqs.at(part).push_back(i);
			}
		}
	}

	return missingSeqs;
}

boolean Utils::subtreeIsMissing( pllInstance * pllTree, vector<int> * missingSequences, const nodeptr node, vector<nodeptr> * missingBranches ) {
	if (find(missingBranches->begin(), missingBranches->end(), node)
			!= missingBranches->end()) {
		return true;
	}
	if (node->number > pllTree->mxtips) {
		if (subtreeIsMissing(pllTree, missingSequences, node->next->back,
				missingBranches)
				& subtreeIsMissing(pllTree, missingSequences,
						node->next->next->back, missingBranches)) {
			missingBranches->push_back(node);
			missingBranches->push_back(node->back);
			return true;
		} else {
			return false;
		}
	} else {
		if (find(missingSequences->begin(), missingSequences->end(), node->number)
				!= missingSequences->end()) {

			missingBranches->push_back(node);
						missingBranches->push_back(node->back);

			/* remove visited taxon */
			missingSequences->erase(
					remove(missingSequences->begin(), missingSequences->end(),
							node->number), missingSequences->end());

			return true;
		} else {
			return false;
		}
	}
}

nodeptr Utils::findRootingNode( pllInstance * pllTree, vector<int> * missingSequences, nodeptr startingNode, vector<nodeptr> * missingBranches ) {

	nodeptr currentNode;
	bool moveRoot = false;

	if (!missingSequences->size()) {
		return(0);
	}

	missingBranches->push_back(startingNode);
	missingBranches->push_back(startingNode->back);
	if (startingNode->back->number <= pllTree->mxtips) {
		missingSequences->erase(
				remove(missingSequences->begin(), missingSequences->end(),
						startingNode->back->number), missingSequences->end());
	}

	/* start searching in a random missing node */
	currentNode = startingNode;

	while (true) {
		bool missingRight = subtreeIsMissing(pllTree, missingSequences, currentNode->next->back, missingBranches);
		bool missingLeft = subtreeIsMissing(pllTree, missingSequences, currentNode->next->next->back, missingBranches);

		if (missingRight && missingLeft) {
			cerr << "ERROR: Everything is missing!!" << endl;
			exit(EX_IOERR);
		} else if (!(missingRight || missingLeft)) {
#ifdef PRINT_TRACE
			cout << "TRACE: Found ancestor in " << currentNode->number << endl;
#endif
			if (moveRoot) {
				missingBranches->push_back(currentNode);
				missingBranches->push_back(currentNode->back);
			}
			return currentNode;
		} else {
			/* move to next position */
			moveRoot = true;
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

static bool compareNodes (nodeptr a, nodeptr b) { return (a->number < b->number); }

std::vector< std::vector<nodeptr> > Utils::findMissingBranches ( pllInstance * pllTree, partitionList * pllPartitions, vector< vector<int> > missingSequences ) {

	vector< vector<nodeptr> > missingBranches(pllPartitions->numberOfPartitions);

	for (int part = 0; part < pllPartitions->numberOfPartitions; part++) {
		while (missingSequences[part].size()) {
			nodeptr startingNode =
					pllTree->nodep[missingSequences[part][0]]->back;
			findRootingNode(pllTree, &missingSequences[part], startingNode,
					&missingBranches[part]);
		}
		sort(missingBranches[part].begin(), missingBranches[part].end(), compareNodes);
	}

	return missingBranches;
}

} /* namespace seqpred */
