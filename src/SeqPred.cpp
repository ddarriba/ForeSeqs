#include "Predictor.h"
#include "Utils.h"
#include "pll.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>
#include <cassert>
#include <getopt.h>

using namespace std;

void optimizeModelParameters(pllInstance * pllTree,
		partitionList * pllPartitions) {
	double lk = 0.0;
	double epsilon = 0.01;
	bool optimizeBranchLengths = pllPartitions->numberOfPartitions > 1;

	if (optimizeBranchLengths) {
		/*
		 * Optimize per-gene branch lengths.
		 */
		int smoothIterations = 64;
		do {
			lk = pllTree->likelihood;
			pllOptimizeBranchLengths(pllTree, pllPartitions, smoothIterations);
			pllOptimizeModelParameters(pllTree, pllPartitions, 0.1);
		} while (fabs(lk - pllTree->likelihood) > epsilon);
	} else {
		/*
		 * In case there is one single partition, we do not optimize the branch lengths.
		 * Otherwise we would have weird results in the branches with missing data.
		 */
		pllOptRatesGeneric(pllTree, pllPartitions, 1.0,
				pllPartitions->rateList);
		pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
				false);
		pllOptBaseFreqs(pllTree, pllPartitions, 1.0, pllPartitions->freqList);
		pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
				false);
		pllOptAlphasGeneric(pllTree, pllPartitions, 1.0,
				pllPartitions->alphaList);
		pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
				false);
		do {
			lk = pllTree->likelihood;
			pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
					false);
			pllOptRatesGeneric(pllTree, pllPartitions, 0.1,
					pllPartitions->rateList);
			pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
					false);
			pllOptBaseFreqs(pllTree, pllPartitions, 0.1,
					pllPartitions->freqList);
			pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
					false);
			pllOptAlphasGeneric(pllTree, pllPartitions, 0.1,
					pllPartitions->alphaList);
			pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
					false);
		} while (fabs(lk - pllTree->likelihood) > epsilon);
	}
}

void exit_with_usage(char * command) {
	printf("\n");
	printf("Usage:\n");
	printf(
			"  %s -i inputFileName -t treeFileName [-q partitionsFileName] [-s seed] [-o outputFileName]\n",
			command);
	printf("\n");
	printf("  -h, --help                           Shows this help message");
	printf("\n");
	printf(
			"  -i, --input      inputFileName       Set the input alignment (mandatory)");
	printf("\n");
	printf(
			"  -o, --output     outputFileName      Set the output filename (default: [inputFile].predicted)");
	printf("\n");
	printf(
			"  -q, --partitions partitionsFileName  Set the partitions definition");
	printf("\n");
	printf(
			"  -s, --seed       randomNumberSed     Set a custom seed (default: 12345)");
	printf("\n");
	printf(
			"  -t, --tree       treeFileName        Set the input tree (mandatory)");
	printf("\n\n");
	printf("Examples:\n");
	printf("  %s -i input_file -t tree_file", command);
	printf("\n\n");
	exit(EXIT_SUCCESS);
}

int main(int argc, char * argv[]) {

	bool partitionsFileDefined = false;

	pllQueue * pllPartsQueue = 0;
	pllInstance * pllTree = 0;
	partitionList * pllPartitions = 0;
	pllAlignmentData * pllAlignment = 0;
	string inputfile, treefile, partitionsfile, outputfile;
	int randomNumberSeed = 12345;

	static struct option long_options[] = { { "help", no_argument, 0, 'h' }, {
			"input", required_argument, 0, 'i' }, { "tree", required_argument,
			0, 't' }, { "partitions", required_argument, 0, 'q' }, { "seed",
			required_argument, 0, 's' },
			{ "output", required_argument, 0, 'o' }, { 0, 0, 0, 0 } };

	int opt = 0, long_index = 0;
	while ((opt = getopt_long(argc, argv, "hi:t:q:n:o:", long_options,
			&long_index)) != -1) {
		switch (opt) {
		case 'h':
			exit_with_usage(argv[0]);
			break;
		case 'i':
			inputfile = optarg;
			break;
		case 't':
			treefile = optarg;
			break;
		case 'q':
			partitionsfile = optarg;
			partitionsFileDefined = true;
			break;
		case 'o':
			outputfile = optarg;
			break;
		case 's':
			randomNumberSeed = atoi(optarg);
			break;
		default:
			exit(EX_IOERR);
		}
	}

	if (argc == 1) {
		inputfile = "testdata/alignment";
		treefile = "testdata/tree";
	} else {
		if (inputfile.length() == 0) {
			cerr << "Alignment file (-i) is required" << endl;
			exit(EX_IOERR);
		}
		if (treefile.length() == 0) {
			cerr << "Tree file (-t) is required" << endl;
			exit(EX_IOERR);
		}
	}

	{
		pllInstanceAttr pllInstanceAttr;
		pllInstanceAttr.fastScaling = PLL_FALSE;
		pllInstanceAttr.randomNumberSeed = randomNumberSeed;
		pllInstanceAttr.rateHetModel = PLL_GAMMA;
		pllInstanceAttr.saveMemory = PLL_FALSE;
		pllInstanceAttr.useRecom = PLL_FALSE;
		pllInstanceAttr.numberOfThreads = 1;
		pllTree = pllCreateInstance(&pllInstanceAttr);
		pllTree->perGeneBranchLengths = PLL_TRUE;
	}

	if (!seqpred::Utils::existsFile(inputfile)) {
		cerr << "[ERROR] Input alignment (" << inputfile << ") does not exist."
				<< endl;
		exit(EX_IOERR);
	}

	if (!seqpred::Utils::existsFile(treefile)) {
		cerr << "[ERROR] Tree file (" << treefile << ") does not exist."
				<< endl;
		exit(EX_IOERR);
	}

	if (partitionsFileDefined && !seqpred::Utils::existsFile(treefile)) {
		cerr << "[ERROR] Partitions file (" << partitionsfile
				<< ") does not exist." << endl;
		exit(EX_IOERR);
	}

	pllAlignment = pllParseAlignmentFile(PLL_FORMAT_PHYLIP, inputfile.c_str());
	if (!pllAlignment) {
		cerr << "[ERROR] There was an error parsing input data." << endl;
		exit(EX_IOERR);
	}
	seqpred::numberOfTaxa = pllAlignment->sequenceCount;
	seqpred::sequenceLength = pllAlignment->sequenceLength;

	if (partitionsFileDefined) {
		pllPartsQueue = pllPartitionParse(partitionsfile.c_str());
		if (!pllPartitionsValidate(pllPartsQueue, pllAlignment)) {
			cerr
					<< "[ERROR] There was an error with the partitions description."
					<< endl;
			exit(EX_IOERR);
		}
	} else {
		char partitionString[256];
		sprintf(partitionString, "DNA, P0 = 1-%d",
				pllAlignment->sequenceLength);
		pllPartsQueue = pllPartitionParseString(partitionString);
		assert(pllPartitionsValidate(pllPartsQueue, pllAlignment));
	}

	pllPartitions = pllPartitionsCommit(pllPartsQueue, pllAlignment);
	pllQueuePartitionsDestroy(&pllPartsQueue);

	if (!pllPartitions) {
		cerr << "[ERROR] There was an error parsing partitions data." << endl;
		exit(EX_IOERR);
	}

#ifdef PRINT_TRACE
	cout << " Make Tree " << endl;
#endif
	pllTreeInitTopologyRandom(pllTree, seqpred::numberOfTaxa,
			pllAlignment->sequenceLabels);

	/* NOTE: We need to initialize the model first. Otherwise fracchange (average subst rate) is 0 */
	cout << endl << "Loading alignment " << endl;
	pllLoadAlignment(pllTree, pllAlignment, pllPartitions);
	seqpred::taxaNames = pllTree->nameList;

	cout << "Initializing model " << endl;
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
		return (EX_IOERR);
	}

	pllTreeInitTopologyNewick(pllTree, nt, PLL_FALSE);

	pllNewickParseDestroy(&nt);

	pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true, false);
//	pllTreeToNewick(pllTree->tree_string, pllTree, pllPartitions,
//			pllTree->start->back, true, true, false, false, false,
//			PLL_SUMMARIZE_LH, false, false);
//	cout << "Tree: " << pllTree->tree_string << endl;

#ifdef PRINT_TRACE
	cout << "TRACE: Initial log likelihood: " << pllTree->likelihood << endl;
#endif

	optimizeModelParameters(pllTree, pllPartitions);
//	pllOptimizeModelParameters(pllTree, pllPartitions, 0.1);

//	pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true, false);
//	pllTreeToNewick(pllTree->tree_string, pllTree, pllPartitions,
//			pllTree->start->back, true, true, false, false, false,
//			PLL_SUMMARIZE_LH, false, false);
//	cout << "Tree: " << pllTree->tree_string << endl;

#ifdef PRINT_TRACE
	cout << "TRACE: Tree=" << pllTree->tree_string << endl;
#endif

	cout << "Predicting sequences..." << endl << endl;
	seqpred::Utils::init();
	for (int currentPartition = 0;
			currentPartition < pllPartitions->numberOfPartitions;
			currentPartition++) {
		seqpred::Predictor sequencePredictor(pllTree, pllPartitions,
				pllAlignment, currentPartition);
		sequencePredictor.predictMissingSequences();
	}
	cout << endl << "...Done!" << endl << endl;
#ifdef PRINT_TRACE
	pllAlignmentDataDumpConsole(pllAlignment);
#endif

	if (outputfile.length() == 0) {
		outputfile = (inputfile + ".prediction");
	}

	pllAlignmentDataDumpFile(pllAlignment, PLL_FORMAT_PHYLIP,
			outputfile.c_str());

	cout << "Alignment dumped to " << outputfile << endl << endl;

	pllAlignmentDataDestroy(pllAlignment);
	pllPartitionsDestroy(pllTree, &pllPartitions);
	pllDestroyInstance(pllTree);

	return (EX_OK);
}

