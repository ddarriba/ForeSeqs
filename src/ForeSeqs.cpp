/**
 * ForeSeqs.cpp
 *
 *  Created on: Oct 2, 2014
 *      Author: Diego Darriba
 *      E-mail: diego.darriba@h-its.org
 *
 *  This file is part of ForeSeqs.
 *
 *  ForeSeqs is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ForeSeqs is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ForeSeqs.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <cmath>
#include <cassert>
#include <getopt.h>
#include <time.h>

#include "PllDefs.h"
#include "Utils.h"
#include "Phylo.h"
#include "Alignment.h"
#include "Predictor.h"

using namespace std;

void exit_with_usage(char *) __attribute__ ((noreturn));


void exit_with_usage(char * command) {
	printf("\n");
	printf("Usage:\n");
	printf(
			"  %s -i inputFileName -t treeFileName [-q partitionsFileName] [-s seed] [-o outputFileName]",
			command);
	printf("\n\n");
	printf("  -b, --branches    branchLengthsMode Set the mode for stealing the branch lengths (default: a)\n");
	printf("      --branches a    (average)       Average of branch lengths in partitions with existing data\n");
	printf("      --branches d    (draw)          Draw a branch length from an inferred distribution\n");
	printf("      --branches s    (scale)         Compute an average branch-length scaler");
	printf("\n\n");
	printf("  -c, --categories  categoriesMode    Set the mode for selecting the per-site category (default: r)\n");
	printf("      --categories r  (random)        Random per-site category for each sequence\n");
	printf("      --categories e  (estimate)      Estimated from the existing data\n");
	printf("      --categories a  (average)       Average of all categories");
	printf("\n\n");
	printf("  -h, --help                           Shows this help message");
	printf("\n\n");
	printf("  -H, --threshold                      Set a threshold in [0.0, 1.0) for considering a missing sequence\n");
	printf("                                       A sequence when the proportion of existing sites is less than the threshold");
		printf("\n\n");
	printf(
			"  -i, --input      inputFileName       Set the input alignment (mandatory)");
	printf("\n\n");
#if(TEST_SIM)
	printf(
			"  -I, --original   originalFileName    Set the original alignment (for testing purposes)");
	printf("\n\n");
#endif
	printf(
			"  -o, --output     outputFileName      Set the output filename (default: [inputFile].predicted)");
	printf("\n\n");
	printf("  -p, --prior      predictionPrior     Set the prior used for predicting the sequences (default: s)\n");
	printf("      --prior s      (sequences)       Predict from the ancestral most likely sequence\n");
	printf("      --prior m      (MAPs)            Predict from the marginal ancestral probabilities (MAPs)");
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
	printf(
			"  -T, --threads    numberOfThreads     Set the number of threads (default: 1)");
	printf("\n\n");
	printf(
			"  -u, --skip-prediction                Skips sequence prediction (prints tree with stolen branches)");
		printf("\n\n");
	printf("Examples:\n");
	printf("  %s -i input_file -t tree_file", command);
	printf("\n\n");
	exit(EXIT_SUCCESS);
}

#define DEFAULT_RANDOM_SEED 12345

int main(int argc, char * argv[]) {

	bool partitionsFileDefined = false;

	time_t startTime, currentTime;

	pll_utree_t * tree = 0;
	foreseqs::Alignment * alignment = 0;
	foreseqs::Phylo * phylo = 0;

	string inputfile, treefile, partitionsfile, outputfile;
	foreseqs::numberOfThreads = 1;
#if(TEST_SIM)
	string originalfile;
	foreseqs::Alignment * originalAlignment = 0;
	bool originalFileDefined = false;
#endif
	unsigned int randomNumberSeed = DEFAULT_RANDOM_SEED;
	unsigned int numberOfReplicates = 1;

	static struct option long_options[] = {
			{ "branches", required_argument, 0, 'b' },
			{ "categories", required_argument, 0, 'c' },
			{ "help", no_argument, 0, 'h' },
			{ "input", required_argument, 0, 'i' },
			{ "original", required_argument, 0, 'I' },
			{ "prior", required_argument, 0, 'p' },
			{ "tree", required_argument, 0, 't' },
			{ "partitions", required_argument, 0, 'q' },
			{ "replicates", required_argument, 0, 'r' },
			{ "seed", required_argument, 0, 's' },
			{ "output", required_argument, 0, 'o' },
			{ "threshold", required_argument, 0, 'H' },
			{ "threads", required_argument, 0, 'T' },
			{ "skip-prediction", required_argument, 0, 'u' },
			{ 0, 0, 0, 0 } };

	int opt = 0, long_index = 0;
	while ((opt = getopt_long(argc, argv, "b:c:hi:I:t:q:r:s:o:p:H:T:u", long_options,
			&long_index)) != -1) {
		switch (opt) {
		case 'b':
				{
					if (strlen(optarg) > 1) {
						cerr << "[ERROR] Invalid branch length stealing mode: " << optarg << endl;
						exit(EX_IOERR);
					}
					char branchMode = optarg[0];
					switch(branchMode) {
					case 'a':
						foreseqs::branchLengthsMode = foreseqs::BL_AVERAGE;
						break;
					case 'd':
						foreseqs::branchLengthsMode = foreseqs::BL_DRAW;
						break;
					case 's':
						foreseqs::branchLengthsMode = foreseqs::BL_SCALE;
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
				foreseqs::categoriesMode = foreseqs::CAT_RANDOM;
				break;
			case 'e':
				foreseqs::categoriesMode = foreseqs::CAT_ESTIMATE;
				break;
			case 'a':
				foreseqs::categoriesMode = foreseqs::CAT_AVERAGE;
				break;
			default:
				cerr << "[ERROR] Invalid categories mode: " << optarg << endl;
				exit(EX_IOERR);
			}
			break;
		}
		case 'h':
			exit_with_usage(argv[0]);
		case 'i':
			inputfile = optarg;
			break;
		case 'I':
#if(TEST_SIM)
			originalfile = optarg;
			originalFileDefined = true;
			break;
#else
			cerr << "[ERROR] -I argument is not available." << endl;
			exit(EX_IOERR);
#endif
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
		case 'p':
		{
			if (strlen(optarg) > 1) {
				cerr << "[ERROR] Invalid prediction mode: " << optarg << endl;
				exit(EX_IOERR);
			}
			char predMode = optarg[0];
			switch (predMode) {
			case 's':
				foreseqs::predictionMode = foreseqs::PRED_ANCSEQ;
				break;
			case 'm':
				foreseqs::predictionMode = foreseqs::PRED_MAP;
				break;
			default:
				cerr << "[ERROR] Invalid prediction mode: " << optarg << endl;
				exit(EX_IOERR);
			}
			break;
		}
		case 'r':
			numberOfReplicates = (unsigned int) atoi(optarg);
			break;
		case 's':
			randomNumberSeed = (unsigned int) atoi(optarg);
			break;
		case 'H':
			foreseqs::threshold = (double) atof(optarg);
			break;
		case 'T':
			foreseqs::numberOfThreads = (unsigned int) atoi(optarg);
			break;
		case 'u':
			foreseqs::predictSequences = false;
			foreseqs::predictionMode = foreseqs::PRED_NONE;
			break;
		default:
			exit(EX_IOERR);
		}
	}

	srand(randomNumberSeed);
	if (argc == 1) {
		/* Test run */
		cout << endl << "****** TEST RUN ******" << endl << endl;
		inputfile      = "testdata/alignment";
		partitionsfile = "testdata/partitions";
		partitionsFileDefined = true;
		treefile       = "testdata/tree";
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

	if (!foreseqs::Utils::existsFile(inputfile)) {
		cerr << "[ERROR] Input alignment (" << inputfile << ") does not exist."
				<< endl;
		exit(EX_IOERR);
	}

	if (!foreseqs::Utils::existsFile(treefile)) {
		cerr << "[ERROR] Tree file (" << treefile << ") does not exist."
				<< endl;
		exit(EX_IOERR);
	}

	if (partitionsFileDefined && !foreseqs::Utils::existsFile(partitionsfile)) {
		cerr << "[ERROR] Partitions file (" << partitionsfile
				<< ") does not exist." << endl;
		exit(EX_IOERR);
	}

	if (outputfile.length() == 0) {
		outputfile = (inputfile + ".prediction");
	}

	if (foreseqs::threshold < 0 || foreseqs::threshold >= 1)
	{
		cerr << "[ERROR] Threshold must be in range [0,1)" << endl;
		exit(EX_IOERR);
	}

	//TODO: try-catch
	alignment = new foreseqs::Alignment(inputfile);
	// if (!alignment) {
	// 	cerr << "[ERROR] There was an error parsing input data." << endl;
	// 	exit(EX_IOERR);
	// }

	foreseqs::numberOfTaxa = (unsigned int) alignment->getSequenceCount();
	foreseqs::sequenceLength = (unsigned int) alignment->getSequenceLength();

	if (!foreseqs::predictSequences) {
		if (numberOfReplicates != 1) {
			cerr << "WARNING: Number of replicates was reset to 1." << endl;
			numberOfReplicates = 1;
		}
		if (foreseqs::predictionMode != foreseqs::PRED_NONE)
		{
			cerr << "WARNING: Ignored prediction algorithm." << endl;
			foreseqs::predictionMode = foreseqs::PRED_NONE;
		}
	}

#if(TEST_SIM)
	if (originalFileDefined && !foreseqs::Utils::existsFile(originalfile)) {
		cerr << "[ERROR] Alignment file (" << originalfile
				<< ") does not exist." << endl;
		exit(EX_IOERR);
	}

	if (originalFileDefined) {
		//TODO: try-catch
		originalAlignment = new Alignment(originalfile);
		// if (!originalAlignment) {
		// 	cerr
		// 			<< "[ERROR] There was an error parsing original alignment data."
		// 			<< endl;
		// 	exit(EX_IOERR);
		// }
		if (foreseqs::numberOfTaxa
				!= (unsigned int) originalAlignment->getSequenceCount()) {
			cerr << "[ERROR] Number of taxa in the original alignment ("
					<< originalAlignment->getSequenceCount()
					<< ") does not match the input file ("
					<< foreseqs::numberOfTaxa << ")." << endl;
			exit(EX_IOERR);
		}
		if (foreseqs::sequenceLength
				!= (unsigned int) originalAlignment->getSequenceLength()) {
			cerr << "[ERROR] Sequence length in the original alignment ("
					<< originalAlignment->getSequenceLength()
					<< ") does not match the input file ("
					<< foreseqs::sequenceLength << ")." << endl;
			exit(EX_IOERR);
		}
	}
#endif

	if (partitionsFileDefined) {
		/* Create partitions */
		if (!alignment->loadPartitionsFile(partitionsfile)) {
			cerr
					<< "[ERROR] There was an error parsing partitions data."
					<< endl;
			exit(EX_IOERR);
		}
	}

	/* Initialization done */

	/* Print header */

	cout << setfill('-') << setw(60) << "" << setfill(' ') << endl;
	char header[60];
	sprintf(header, " ForeSeqs v%s", PACKAGE_VERSION);
	int padding = (int) (30 + (strlen(header) / 2));
	cout << setw(padding) << header << endl;
	cout << setfill('-') << setw(60) << "" << setfill(' ') << endl;
	cout << setw(20) << left << "Input alignment:" << inputfile << endl;
	cout << setw(20) << left << "Input tree:" << treefile << endl;
	cout << setw(20) << left << "Partitions file:" << ((partitionsfile.length() > 0)?partitionsfile:"-") << endl;
	cout << setw(20) << left << "Output file:" << outputfile << endl;
	if (foreseqs::threshold > 0)
		cout << setw(20) << left << "Threshold:" << foreseqs::threshold << endl;
	cout << setw(20) << left << "Num.Taxa:" << foreseqs::numberOfTaxa << endl;
	cout << setw(20) << left << "Seq.Length:" << foreseqs::sequenceLength << endl;
	cout << setw(20) << left << "Branch lengths:";
	switch (foreseqs::branchLengthsMode) {
	case foreseqs::BL_AVERAGE:
		cout << "Average" << endl;
		break;
	case foreseqs::BL_DRAW:
		cout << "Drawn" << endl;
		break;
	case foreseqs::BL_SCALE:
		cout << "Average scaler" << endl;
		break;
	}
	cout << setw(20) << left << "Prediction mode:";
	switch (foreseqs::predictionMode) {
	case foreseqs::PRED_ANCSEQ:
		cout << "Ancestral sequences" << endl;
		break;
	case foreseqs::PRED_MAP:
		cout << "Marginal ancestral probabilities" << endl;
		break;
	case foreseqs::PRED_NONE:
		cout << "No prediction" << endl;
		break;
	}
	cout << setw(20) << left << "Gamma rates:";
	switch(foreseqs::categoriesMode) {
	case foreseqs::CAT_RANDOM:
		cout<< "Random" << endl;
		break;
	case foreseqs::CAT_AVERAGE:
			cout<< "Average" << endl;
			break;
	case foreseqs::CAT_ESTIMATE:
			cout<< "Estimated" << endl;
			break;
	}
	cout << setw(20) << left << "Random seed:" << randomNumberSeed << endl;
	cout << setw(20) << left << "#Replicates:" << numberOfReplicates << endl;
	cout << setfill('-') << setw(60) << "" << setfill(' ') << endl;
	cout << "ForeSeqs was called as follows:" << endl << "  ";
	for (int i=1; i<argc; i++) cout << argv[i] << " ";
	cout << endl << endl;

	/* Check partitions / Warn if needed */

	if (alignment->getNumberOfPartitions() == 1) {
		cout << endl << setfill('*') << setw(54) << "" << setfill(' ') << endl
				<< "WARNING: There is only one partition." << endl
				<< "         Branch lengths are taken from the input tree." << endl
				<< setfill('*') << setw(54) << "" << setfill(' ') << endl;
	}

	/* Start! */
	startTime = time(NULL);

	unsigned int testNumTaxa;
	tree = pll_utree_parse_newick(treefile.c_str(), &testNumTaxa);
	if (!tree) {
		cerr << "[ERROR] There was an error parsing newick file " << treefile << endl;
		return (EXIT_FAILURE);
	}

	if (foreseqs::numberOfTaxa != testNumTaxa)
	{
		cerr << "[ERROR] Number of taxa mismatch " << treefile << endl;
		return (EXIT_FAILURE);
	}

  foreseqs::taxaNames = alignment->getLabels();


	// foreseqs::taxaNames = pllTree->nameList;
	//
	// /* build translation array */
	// foreseqs::seqIndexTranslate = (unsigned int *) malloc ((foreseqs::numberOfTaxa+1)*sizeof(unsigned int));
	// for (unsigned int i=1; i<=foreseqs::numberOfTaxa; i++) {
	// 	for (unsigned int j=1; j<=foreseqs::numberOfTaxa; j++) {
	// 		if (!strcmp(pllAlignment->sequenceLabels[j], pllTree->tipNames[i])) {
	// 			foreseqs::seqIndexTranslate[i] = j;
	// 		}
	// 	}
	// }

	cout << "Initializing model " << endl;

	phylo = new foreseqs::Phylo(*alignment, tree, 4);

	for (size_t part =0; part < alignment->getNumberOfPartitions(); part++)
		double logLikelihood = phylo->optimizeModelParameters(part, 0.01, 0.01);

 	currentTime = time(NULL);
 	cout << "Initial inference done. It took " << currentTime - startTime << " seconds." << endl << endl;

// #ifdef PRINT_TRACE
// 	pllTreeToNewick(pllTree->tree_string, pllTree, pllPartitions,
// 			pllTree->start->back, true, true, true, false, false,
// 			PLL_SUMMARIZE_LH, false, false);
// 	cout << pllTree->tree_string << endl;
// 	for (size_t i = 0; i < missingBranches.size(); i++) {
// 		cout << "PARTITION " << i << "/" << missingBranches.size() << endl;
// 		cout << "  Branches: " << missingBranches[i].size() << endl;
// 		for (size_t j = 0; j < missingBranches[i].size(); j++) {
// 			cout << "    * " << missingBranches[i][j]->number << endl;
// 		}
// 	}
// #endif
//
// #if(TEST_SIM)
// 	cout << endl << "T(ini,avg): ";
// 	pllTreeToNewick(pllTree->tree_string, pllTree, pllPartitions,
// 					pllTree->start->back, true, true, true, false, false,
// 					PLL_SUMMARIZE_LH, false, false);
// 	cout << pllTree->tree_string << endl;
// #endif
//
// 	/* Predict sequences */
 	cout << "Predicting sequences..." << endl << endl;
 	for (unsigned int rep = 0; rep < numberOfReplicates; rep++) {
//
// 		if (numberOfReplicates > 1) {
// 			cout << "Replicate " << rep+1 << " of " << numberOfReplicates << endl;
// 		}
//
		for (unsigned int currentPartition = 0;
				currentPartition < (unsigned int) alignment->getNumberOfPartitions();
				currentPartition++)
		{
			foreseqs::Predictor sequencePredictor(tree, *phylo,
					*alignment, currentPartition);

		 	sequencePredictor.predictMissingSequences();
		}
//
// #if(TEST_SIM)
// 			cout << endl << "T(ini," << currentPartition << "): ";
// 			pllTreeToNewick(pllTree->tree_string, pllTree, pllPartitions,
// 					pllTree->start->back, true, true, true, false, false,
// 					currentPartition, false, false);
// 			cout << pllTree->tree_string << endl;
//
// 			if (originalFileDefined) {
// 				sequencePredictor.predictMissingSequences(originalAlignment);
// 				if (sequencePredictor.getMissingPartsCount()) {
// 					cout << "Similarity in partition "
// 							<< pllPartitions->partitionData[currentPartition]->partitionName
// 							<< ": "
// 							<< 100 * sequencePredictor.getSequenceSimilarity()
// 							<< "%" << endl;
// 				}
// 			} else {
// 				sequencePredictor.predictMissingSequences();
// 			}
//
// 			cout << endl << "T(end," << currentPartition << "): ";
// 			pllTreeToNewick(pllTree->tree_string, pllTree, pllPartitions,
// 					pllTree->start->back, true, true, true, false, false,
// 					currentPartition, false, false);
// 			cout << pllTree->tree_string << endl;
// 		}
//
// 		cout << endl << "T(end,avg): ";
// 		pllTreeToNewick(pllTree->tree_string, pllTree, pllPartitions,
// 				pllTree->start->back, true, true, true, false, false,
// 				PLL_SUMMARIZE_LH, false, false);
// 		cout << pllTree->tree_string << endl;
// #else
// 			sequencePredictor.predictMissingSequences();
// 		}
// #endif
//
 		currentTime = time(NULL);
 		cout << endl << "Prediction done. It took " << currentTime - startTime << " seconds." << endl << endl;
//
// 	#ifdef PRINT_TRACE
// 		pllAlignmentDataDumpConsole(pllAlignment);
// 	#endif
//
// 		cout << "Summarized tree (stolen branches): " << endl;
// 		pllTreeToNewick(pllTree->tree_string, pllTree, pllPartitions,
// 				pllTree->start->back, true, true, true, false, false,
// 				PLL_SUMMARIZE_LH, false, false);
// 		cout << pllTree->tree_string << endl;
//
// 		if (foreseqs::predictSequences) {
// 			stringstream rep_outputfile;
// 			rep_outputfile << outputfile;
// 			if (numberOfReplicates > 1) {
// 				rep_outputfile << ".R" << rep;
// 			}
// 			pllAlignmentDataDumpFile(pllAlignment, PLL_FORMAT_PHYLIP,
// 					rep_outputfile.str().c_str());
//
// 			cout << "Alignment dumped to " << rep_outputfile.str() << endl << endl;
// 		}
	}
//
// #if(TEST_SIM)
// 			if (originalFileDefined) {
// 				pllAlignmentDataDestroy(originalAlignment);
// 			}
// #endif
//
// 	free(foreseqs::seqIndexTranslate);
// 	pllAlignmentDataDestroy(pllAlignment);
// 	pllPartitionsDestroy(pllTree, &pllPartitions);
// 	pllDestroyInstance(pllTree);

	return (EX_OK);
}
