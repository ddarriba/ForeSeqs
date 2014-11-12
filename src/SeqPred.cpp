#include "Predictor.h"
#include "Utils.h"

#include "pll.h"
#include "config.h"

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <cmath>
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
		cout << "Optimizing per-gene branch lengths / model parameters " << endl;

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
		cout << "Optimizing model parameters " << endl;

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
			"  %s -i inputFileName -t treeFileName [-q partitionsFileName] [-s seed] [-o outputFileName]",
			command);
	printf("\n\n");
	printf("  -b, --branches    branchLengthsMode Set the mode for stealing the branch lengths (default: a)\n");
	printf("      --branches    a (average)       Average of branch lengths for other partitions\n");
	printf("      --branches    d (draw)          Draw a branch length from an inferred distribution\n");
	printf("      --branches    s (scale)         Compute an average branch-length scaler\n");
	printf("\n\n");
	printf("  -c, --categories  categoriesMode    Set the mode for selecting the per-site category (default: r)\n");
	printf("      --categories  r (random)        Random per-site category for each sequence\n");
	printf("      --categories  e (estimate)      Estimated from other partitions\n");
	printf("      --categories  a (average)       Average of all categories\n");
	printf("\n\n");
	printf("  -h, --help                           Shows this help message");
	printf("\n\n");
	printf(
			"  -i, --input      inputFileName       Set the input alignment (mandatory)");
	printf("\n\n");
#ifdef TEST_SIM
	printf(
			"  -I, --original   originalFileName    Set the original alignment (for testing purposes)");
	printf("\n\n");
#endif
	printf(
			"  -o, --output     outputFileName      Set the output filename (default: [inputFile].predicted)");
	printf("\n\n");
	printf(
			"  -q, --partitions partitionsFileName  Set the partitions definition (PLL like)");
	printf("\n\n");
	printf(
			"                   If no partitions file is set, one single partition and DNA data is assumed\n");
	printf(
				"                   and branch lengths are taken directly from the input tree.");
	printf("\n\n");
	printf(
			"  -r, --replicates numberOfReplicates  Produce several replicates");
	printf("\n\n");
	printf(
			"  -s, --seed       randomNumberSed     Set a custom seed (default: 12345)");
	printf("\n\n");
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
#ifdef TEST_SIM
	string originalfile;
	pllAlignmentData * originalAlignment = 0;
	bool originalFileDefined = false;
#endif
	unsigned int randomNumberSeed = 12345;
	unsigned int numberOfReplicates = 1;

	static struct option long_options[] = {
			{ "branches", required_argument, 0, 'b' },
			{ "categories", required_argument, 0, 'c' },
			{ "help", no_argument, 0, 'h' },
			{ "input", required_argument, 0, 'i' },
			{ "original", required_argument, 0, 'I' },
			{ "tree", required_argument, 0, 't' },
			{ "partitions", required_argument, 0, 'q' },
			{ "replicates", required_argument, 0, 'r' },
			{ "seed", required_argument, 0, 's' },
			{ "output", required_argument, 0, 'o' },
			{ 0, 0, 0, 0 } };

	int opt = 0, long_index = 0;
	while ((opt = getopt_long(argc, argv, "b:c:hi:I:t:q:r:s:o:", long_options,
			&long_index)) != -1) {
		switch (opt) {
		case 'b':
				{
					if (strlen(optarg) > 1) {
						cerr << "[ERROR] Invalid branch length stealing mode: " << optarg << endl;
						exit(EX_IOERR);
					}
					char catMode = optarg[0];
					switch(catMode) {
					case 'a':
						seqpred::branchLengthsMode = seqpred::BL_AVERAGE;
						break;
					case 'd':
						seqpred::branchLengthsMode = seqpred::BL_DRAW;
						break;
					case 's':
						seqpred::branchLengthsMode = seqpred::BL_SCALE;
						break;
					default:
						cerr << "[ERROR] Invalid branch length stealing mode: " << optarg << endl;
						exit(EX_IOERR);
					}
					break;
				}
		case 'c':
		{
			if (strlen(optarg) > 1) {
				cerr << "[ERROR] Invalid categories mode: " << optarg << endl;
				exit(EX_IOERR);
			}
			char catMode = optarg[0];
			switch(catMode) {
			case 'r':
				seqpred::categoriesMode = seqpred::CAT_RANDOM;
				break;
			case 'e':
				seqpred::categoriesMode = seqpred::CAT_ESTIMATE;
				break;
			case 'a':
				seqpred::categoriesMode = seqpred::CAT_AVERAGE;
				break;
			default:
				cerr << "[ERROR] Invalid categories mode: " << optarg << endl;
				exit(EX_IOERR);
			}
			break;
		}
		case 'h':
			exit_with_usage(argv[0]);
			break;
		case 'i':
			inputfile = optarg;
			break;
		case 'I':
#ifdef TEST_SIM
			originalfile = optarg;
			originalFileDefined = true;
#else
			cerr << "[ERROR] -I argument is not available." << endl;
			exit(EX_IOERR);
#endif
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
		case 'r':
			numberOfReplicates = atoi(optarg);
			break;
		case 's':
			randomNumberSeed = atoi(optarg);
			break;
		default:
			exit(EX_IOERR);
		}
	}

	srand(randomNumberSeed);
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

	if (partitionsFileDefined && !seqpred::Utils::existsFile(partitionsfile)) {
		cerr << "[ERROR] Partitions file (" << partitionsfile
				<< ") does not exist." << endl;
		exit(EX_IOERR);
	}

	if (outputfile.length() == 0) {
		outputfile = (inputfile + ".prediction");
	}

	pllAlignment = pllParseAlignmentFile(PLL_FORMAT_PHYLIP, inputfile.c_str());
	if (!pllAlignment) {
		cerr << "[ERROR] There was an error parsing input data." << endl;
		exit(EX_IOERR);
	}
	seqpred::numberOfTaxa = pllAlignment->sequenceCount;
	seqpred::sequenceLength = pllAlignment->sequenceLength;

#ifdef TEST_SIM
	if (originalFileDefined && !seqpred::Utils::existsFile(originalfile)) {
		cerr << "[ERROR] Alignment file (" << originalfile
				<< ") does not exist." << endl;
		exit(EX_IOERR);
	}

	if (originalFileDefined) {
		originalAlignment = pllParseAlignmentFile(PLL_FORMAT_PHYLIP,
				originalfile.c_str());
		if (!originalAlignment) {
			cerr
					<< "[ERROR] There was an error parsing original alignment data."
					<< endl;
			exit(EX_IOERR);
		}
		if (seqpred::numberOfTaxa
				!= (unsigned int) originalAlignment->sequenceCount) {
			cerr << "[ERROR] Number of taxa in the original alignment ("
					<< originalAlignment->sequenceCount
					<< ") does not match the input file ("
					<< seqpred::numberOfTaxa << ")." << endl;
			exit(EX_IOERR);
		}
		if (seqpred::sequenceLength
				!= (unsigned int) originalAlignment->sequenceLength) {
			cerr << "[ERROR] Sequence length in the original alignment ("
					<< originalAlignment->sequenceLength
					<< ") does not match the input file ("
					<< seqpred::sequenceLength << ")." << endl;
			exit(EX_IOERR);
		}
	}
#endif

	if (partitionsFileDefined) {
		/* Create partitions */
		pllPartsQueue = pllPartitionParse(partitionsfile.c_str());
		if (!pllPartitionsValidate(pllPartsQueue, pllAlignment)) {
			cerr
					<< "[ERROR] There was an error with the partitions description."
					<< endl;
			exit(EX_IOERR);
		}
	} else {
		/* Set one single partition */
		char partitionString[20];
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

	/* Initialization done */

	/* Print header */

	cout << setfill('-') << setw(60) << "" << setfill(' ') << endl;
	char header[60];
	sprintf(header, " Sequence predictor v%s", PACKAGE_VERSION);
	unsigned int padding = 30 + (strlen(header) / 2);
	cout << setw(padding) << header << endl;
	cout << setfill('-') << setw(60) << "" << setfill(' ') << endl;
	cout << setw(20) << left << "Input alignment:" << inputfile << endl;
	cout << setw(20) << left << "Input tree:" << treefile << endl;
	cout << setw(20) << left << "Partitions file:" << ((partitionsfile.length() > 0)?partitionsfile:"-") << endl;
	cout << setw(20) << left << "Output file:" << outputfile << endl;
	cout << setw(20) << left << "Gamma rates:";
	switch(seqpred::categoriesMode) {
	case seqpred::CAT_RANDOM:
		cout<< "Random" << endl;
		break;
	case seqpred::CAT_AVERAGE:
			cout<< "Average" << endl;
			break;
	case seqpred::CAT_ESTIMATE:
			cout<< "Estimated" << endl;
			break;
	}
	cout << setw(20) << left << "Random seed:" << randomNumberSeed << endl;
	cout << setfill('-') << setw(60) << "" << setfill(' ') << endl;

	/* Check partitions / Warn if needed */

	if (pllPartitions->numberOfPartitions == 1) {
		cout << endl << setfill('*') << setw(54) << "" << setfill(' ') << endl
				<< "WARNING: There is only one partition." << endl
				<< "         Branch lengths are taken from the input tree." << endl
				<< setfill('*') << setw(54) << "" << setfill(' ') << endl;
	} else {
		pllPartitions->perGeneBranchLengths = true; /* IMPORTANT!! */
	}

	/* Start! */

	pllTreeInitTopologyRandom(pllTree, seqpred::numberOfTaxa,
			pllAlignment->sequenceLabels);

	cout << endl << "Loading alignment " << endl;
	pllLoadAlignment(pllTree, pllAlignment, pllPartitions);
	seqpred::taxaNames = pllTree->nameList;

	cout << "Initializing model " << endl;
	pllInitModel(pllTree, pllPartitions);

	pllNewickTree * nt;
	nt = pllNewickParseFile(treefile.c_str());
	if (!nt) {
		cerr << "[ERROR] There was an error parsing newick file " << treefile << endl;
		return (EXIT_FAILURE);
	}
	if (!pllValidateNewick(nt))
	{
		cerr << "Invalid phylogenetic tree" << endl;
		printf("%d\n", errno);
		return (EX_IOERR);
	}

	cout << "Seting fixed topology " << endl;
	pllTreeInitTopologyNewick(pllTree, nt, PLL_FALSE);

	pllNewickParseDestroy(&nt);

	pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true, false);

#ifdef PRINT_TRACE
	cout << "TRACE: Initial log likelihood: " << pllTree->likelihood << endl;
#endif

	optimizeModelParameters(pllTree, pllPartitions);

#ifdef PRINT_TRACE
	cout << "TRACE: Tree=" << pllTree->tree_string << endl;
#endif


	cout << "Predicting sequences..." << endl << endl;
	for (unsigned int rep = 0; rep < numberOfReplicates; rep++) {

		if (numberOfReplicates > 1) {
			cout << "Replicate " << rep+1 << " of " << numberOfReplicates << endl;
		}
		for (int currentPartition = 0;
				currentPartition < pllPartitions->numberOfPartitions;
				currentPartition++) {
			seqpred::Predictor sequencePredictor(pllTree, pllPartitions,
					pllAlignment, currentPartition);
#ifdef TEST_SIM
			if (originalFileDefined) {
				sequencePredictor.predictMissingSequences(originalAlignment);
				if (sequencePredictor.getMissingPartsCount()) {
					cout << "Similarity in partition "
							<< pllPartitions->partitionData[currentPartition]->partitionName
							<< ": "
							<< 100 * sequencePredictor.getSequenceSimilarity()
							<< "%" << endl;
				}
			} else {
				sequencePredictor.predictMissingSequences();
			}
#else
			sequencePredictor.predictMissingSequences();
#endif
		}
		cout << endl << "...Done!" << endl << endl;
	#ifdef PRINT_TRACE
		pllAlignmentDataDumpConsole(pllAlignment);
	#endif

		stringstream rep_outputfile;
		rep_outputfile << outputfile;
		if (numberOfReplicates > 1) {
			rep_outputfile << ".R" << rep;
		}
		pllAlignmentDataDumpFile(pllAlignment, PLL_FORMAT_PHYLIP,
				rep_outputfile.str().c_str());

		cout << "Alignment dumped to " << rep_outputfile.str() << endl << endl;
	}

#ifdef TEST_SIM
			if (originalFileDefined) {
				pllAlignmentDataDestroy(originalAlignment);
			}
#endif
	pllAlignmentDataDestroy(pllAlignment);
	pllPartitionsDestroy(pllTree, &pllPartitions);
	pllDestroyInstance(pllTree);

	return (EX_OK);
}

