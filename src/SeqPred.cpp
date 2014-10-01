#include "Predictor.h"
#include "pll.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

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
	pllQueuePartitionsDestroy(&pllPartsQueue);

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
	pllAlignmentDataDestroy(pllPhylip);

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

//	for (int i = 1; i <= (2 * pllTree->ntips - 3); i++) {
//		double bl = pllGetBranchLength(pllTree, pllTree->nodep[i], 0);
//		cout << pllTree->nodep[i]->z[0] << " " << bl << endl;
//	}
	pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true, false);
	cout << "TRACE: Initial log likelihood: " << pllTree->likelihood << endl;
	pllOptimizeModelParameters(pllTree, pllPartitions, 0.1);

	seqpred::Predictor sequencePredictor(pllTree, pllPartitions, 0);
	sequencePredictor.predictMissingSequences();
	//predictMissingSequences(pllTree, pllPartitions);

	pllPartitionsDestroy(pllTree, &pllPartitions);
	pllDestroyInstance(pllTree);
}

