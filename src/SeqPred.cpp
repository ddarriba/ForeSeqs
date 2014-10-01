#include "Predictor.h"
#include "pll.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

///*
// * Search if all the taxa in a subtree have missing sequences.
// * TODO: This algorithm can be improved using a hashtable for storing the already evaluated nodes.
// */
//boolean subtreeIsMissing(nodeptr node, vector<int> missing, pllInstance * tree) {
//	if (node->number > tree->mxtips) {
//		return (
//				subtreeIsMissing(node->next->back, missing, tree)
//				&
//				subtreeIsMissing(node->next->next->back, missing, tree)
//				);
//	} else {
//		return (find(missing.begin(), missing.end(), node->number) != missing.end());
//	}
//}

//nodeptr missingDataAncestor(vector<int> missing, pllInstance * tree) {
//
//	nodeptr currentNode;
//
//	if (!missing.size()) {
//		cerr << "ERROR: No missing sequences" << endl;
//		return 0;
//	}
//
//	/* start searching in a random missing node */
//	currentNode = tree->nodep[missing[0]]->back;
//
//	while(true) {
//		bool missingRight = subtreeIsMissing(currentNode->next->back, missing, tree);
//		bool missingLeft = subtreeIsMissing(currentNode->next->next->back, missing, tree);
//
//		if (missingRight && missingLeft) {
//			cerr << "ERROR: Everything is missing!!" << endl;
//			return 0;
//		} else if (!(missingRight || missingLeft)) {
//			cout << "INFO: Found ancestor in " << currentNode->number << endl;
//			return currentNode;
//		} else {
//			/* move to next position */
//			currentNode = (!missingRight)?currentNode->next->back:currentNode->next->next->back;
//			cout << "INFO: Moving node to " << currentNode->number << endl;
//		}
//	}
//	return 0;
//}

//vector<int> getMissingSequences(pllInstance * tree, partitionList * partitions,
//		int start, int end) {
//	vector<int> missingSeqs;
//
//	int missing;
//	for (int i = 1; i <= tree->ntips; i++) {
//		missing = 1;
//		for (int j = start; j < end; j++) {
//			missing &= (tree->yVector[i][j] == PLL_UNDEFINED_SITE);
//		}
//		if (missing) {
//			missingSeqs.push_back(i);
//		}
//	}
//	return missingSeqs;
//}

//void mutateSequence ( char * currentSequence, char * ancestralSequence ) {
//	/* TODO : THis */
//}
//
//void evolveNode(nodeptr node, pllInstance * tree, char * ancestralSequence) {
//
//	char * currentSequence;
//	mutateSequence(currentSequence, ancestralSequence);
//	cout << "INFO: Mutating node " << node->number << endl;
//
//	if (node->number > tree->mxtips) {
//		cout << "   Recurse for Mutating node " << node->next->back->number << " and " << node->next->next->back->number << endl;
//		evolveNode(node->next->back, tree, currentSequence);
//		evolveNode(node->next->next->back, tree, currentSequence);
//	}
//}
//
//void predictMissingSequences(pllInstance * tree, partitionList * partitions) {
//
//	int start = partitions->partitionData[0]->lower;
//	int end = partitions->partitionData[0]->upper;
//	int seqlen = end - start + 1;
//
//	vector<int> missingSeqs = getMissingSequences(tree, partitions, start, end);
//	for (int i = 0; i < missingSeqs.size(); i++) {
//		cout << "Missing sequence " << missingSeqs[i] << endl;
//	}
//
//	nodeptr ancestor = missingDataAncestor(missingSeqs, tree);
//	evolveNode(ancestor->back, tree, 0);
//
//	nodeptr startNode = ancestor; //tree->nodep[5]->back->next->next->back;
//	pllEvaluateLikelihood(tree, partitions, tree->start, true, false);
//	pllTreeToNewick(tree->tree_string, tree, partitions, tree->start->back, true, true, false, false, false, PLL_SUMMARIZE_LH, false, false);
//	cout << "Tree: " << tree->tree_string << endl;
//	cout << "Updating partials for " << startNode->number << endl;
//	pllUpdatePartialsAncestral(tree, partitions, startNode );
//	cout << "Generating ancestral for " << startNode->number << endl;
//	char * ancestral = (char *) malloc (seqlen + 1);
//	double * probs = (double *) malloc (seqlen * partitions->partitionData[0]->states * sizeof(double));
//	pllGetAncestralState(tree, partitions, startNode, probs, ancestral);
//	cout << "Got ancestral" << endl;
//	cout << "Ancestral Sequence " << ancestral << endl;
//
//
//}

int main(int argc, char * argv[]) {

	pllQueue * pllPartsQueue = 0;
	pllInstance * pllTree = 0;
	partitionList * pllPartitions = 0;
	pllAlignmentData * pllPhylip = 0;
	string filename("testdata/alignment");
	string treefile("testdata/tree");

	{
		pllInstanceAttr pllInstanceAttr;
		pllInstanceAttr.fastScaling = PLL_FALSE;
		pllInstanceAttr.randomNumberSeed = 12345;
		pllInstanceAttr.rateHetModel = PLL_GAMMA;
		pllInstanceAttr.saveMemory = PLL_FALSE;
		pllInstanceAttr.useRecom = PLL_FALSE;
		pllInstanceAttr.numberOfThreads = 1;
		pllTree = pllCreateInstance(&pllInstanceAttr);
	}

	pllPhylip = pllParseAlignmentFile(PLL_FORMAT_PHYLIP, filename.c_str());
	if (!pllPhylip) {
		cerr << "[ERROR] There was an error parsing input data." << endl;
		exit(1);
	}

	char partitionString[256];
	sprintf(partitionString, "DNA, P0 = 1-%d", pllPhylip->sequenceLength);

	pllPartsQueue = pllPartitionParseString(partitionString);
	pllPartitions = pllPartitionsCommit(pllPartsQueue, pllPhylip);
	if (!pllPartitions) {
		cerr << "[ERROR] There was an error parsing partitions data." << endl;
		exit(1);
	}

	cout << "TRACE: Make Tree " << endl;
	pllTreeInitTopologyRandom(pllTree, pllPhylip->sequenceCount,
			pllPhylip->sequenceLabels);

	/* NOTE: We need to initialize the model first. Otherwise fracchange (average subst rate) is 0 */
	cout << "TRACE: Load Alignment " << endl;
	pllLoadAlignment(pllTree, pllPhylip, pllPartitions);
	cout << "TRACE: Initialize Model " << endl;
	pllInitModel(pllTree, pllPartitions);

	pllNewickTree * nt;
	nt = pllNewickParseFile(treefile.c_str());
	if (!nt) {
		fprintf(stderr, "Error while parsing newick file %s\n", argv[2]);
		return (EXIT_FAILURE);
	}
	if (!pllValidateNewick(nt)) /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */
	{
		fprintf(stderr, "Invalid phylogenetic tree\n");
		printf("%d\n", errno);
		return (EXIT_FAILURE);
	}

	pllTreeInitTopologyNewick(pllTree, nt, PLL_FALSE);

	pllNewickParseDestroy(&nt);

	for (int i = 1; i <= (2 * pllTree->ntips - 3); i++) {
		double bl = pllGetBranchLength(pllTree, pllTree->nodep[i], 0);
		cout << pllTree->nodep[i]->z[0] << " " << bl << endl;
	}

	seqpred::Predictor sequencePredictor(pllTree, pllPartitions, 0);
	sequencePredictor.predictMissingSequences();
	//predictMissingSequences(pllTree, pllPartitions);
}

