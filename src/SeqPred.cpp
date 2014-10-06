#include "Predictor.h"
#include "Utils.h"
#include "pll.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;

void optimizeModelParameters(pllInstance * pllTree,
		partitionList * pllPartitions) {
	pllOptRatesGeneric(pllTree, pllPartitions, 1.0, pllPartitions->rateList);
	pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true, false);
	pllOptBaseFreqs(pllTree, pllPartitions, 1.0, pllPartitions->freqList);
	pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true, false);
	pllOptAlphasGeneric(pllTree, pllPartitions, 1.0, pllPartitions->alphaList);
	pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true, false);
	double lk = 0.0;
	do {
		lk = pllTree->likelihood;
		pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
				false);
		pllOptRatesGeneric(pllTree, pllPartitions, 0.1,
				pllPartitions->rateList);
		pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
				false);
		pllOptBaseFreqs(pllTree, pllPartitions, 0.1, pllPartitions->freqList);
		pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
				false);
		pllOptAlphasGeneric(pllTree, pllPartitions, 0.1,
				pllPartitions->alphaList);
		pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
				false);
	} while (fabs(lk - pllTree->likelihood) > 0.01);
}

int main(int argc, char * argv[]) {

	int numberOfTaxa, sequenceLength;

	pllQueue * pllPartsQueue = 0;
	pllInstance * pllTree = 0;
	partitionList * pllPartitions = 0;
	pllAlignmentData * pllPhylip = 0;
	string inputfile, treefile, partitionsfile;

	if (argc == 1) {
		inputfile = "testdata/alignment";
		treefile = "testdata/tree";
	} else {
		if (argc != 4) {
			cerr << "Usage: " << argv[0] << " seqFile treeFile partitionsFile";
			inputfile = argv[1];
			treefile = argv[2];
			partitionsfile = argv[3];
		}
	}

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

	pllPhylip = pllParseAlignmentFile(PLL_FORMAT_PHYLIP, inputfile.c_str());
	if (!pllPhylip) {
		cerr << "[ERROR] There was an error parsing input data." << endl;
		exit(1);
	}
	numberOfTaxa = pllPhylip->sequenceCount;
	sequenceLength = pllPhylip->sequenceLength;

	char partitionString[256];
	sprintf(partitionString, "DNA, P0 = 1-%d", pllPhylip->sequenceLength);

	pllPartsQueue = pllPartitionParseString(partitionString);
	pllPartitions = pllPartitionsCommit(pllPartsQueue, pllPhylip);
	pllQueuePartitionsDestroy(&pllPartsQueue);

	if (!pllPartitions) {
		cerr << "[ERROR] There was an error parsing partitions data." << endl;
		exit(1);
	}

#ifdef PRINT_TRACE
	cout << "TRACE: Make Tree " << endl;
#endif
	pllTreeInitTopologyRandom(pllTree, numberOfTaxa,
			pllPhylip->sequenceLabels);

	/* NOTE: We need to initialize the model first. Otherwise fracchange (average subst rate) is 0 */
#ifdef PRINT_TRACE
	cout << "TRACE: Load Alignment " << endl;
#endif
	pllLoadAlignment(pllTree, pllPhylip, pllPartitions);
	pllAlignmentDataDestroy(pllPhylip);

#ifdef PRINT_TRACE
	cout << "TRACE: Initialize Model " << endl;
#endif
	pllInitModel(pllTree, pllPartitions);

	pllNewickTree * nt;
	nt = pllNewickParseFile(treefile.c_str());
	if (!nt) {
		cerr << "Error while parsing newick file " << argv[2] << endl;
		return (EXIT_FAILURE);
	}
	if (!pllValidateNewick(nt)) /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */
	{
		cerr << "Invalid phylogenetic tree" << endl;
		printf("%d\n", errno);
		return (EXIT_FAILURE);
	}

	pllTreeInitTopologyNewick(pllTree, nt, PLL_FALSE);

	pllNewickParseDestroy(&nt);

	pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true, false);
	pllTreeToNewick(pllTree->tree_string, pllTree, pllPartitions,
			pllTree->start->back, true, true, false, false, false,
			PLL_SUMMARIZE_LH, false, false);
	cout << "Tree: " << pllTree->tree_string << endl;

	pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true, false);
#ifdef PRINT_TRACE
	cout << "TRACE: Initial log likelihood: " << pllTree->likelihood << endl;
#endif

	optimizeModelParameters(pllTree, pllPartitions);
//	pllOptimizeModelParameters(pllTree, pllPartitions, 0.1);

	pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true, false);
	pllTreeToNewick(pllTree->tree_string, pllTree, pllPartitions,
			pllTree->start->back, true, true, false, false, false,
			PLL_SUMMARIZE_LH, false, false);

#ifdef PRINT_TRACE
	cout << "TRACE: Tree=" << pllTree->tree_string << endl;
#endif

	seqpred::Utils::init();
	for (int currentPartition=0; currentPartition < pllPartitions->numberOfPartitions; currentPartition++) {
		seqpred::Predictor sequencePredictor(pllTree, pllPartitions, 0);
		sequencePredictor.predictMissingSequences();
		//predictMissingSequences(pllTree, pllPartitions);

		if (sequencePredictor.getNumberOfMissingSequences() > 0) {
			vector<int> missingSequences = sequencePredictor.getMissingSequences();
			map<int, char*> predictedSequences =
					sequencePredictor.getPredictedSequences();
			for (int i = 1; i <= numberOfTaxa; i++) {
				cout << i << " : ";
				if (predictedSequences.find(i) != predictedSequences.end()) {
					seqpred::Utils::printSequence(predictedSequences[i]);
				} else {
					seqpred::Utils::printSequence(pllPartitions->partitionData[currentPartition]->yVector[i], pllPartitions->partitionData[currentPartition]->width);
				}
			}
		}
	}

	pllPartitionsDestroy(pllTree, &pllPartitions);
	pllDestroyInstance(pllTree);
}

