/*
 * Utils.cpp
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

#include "Utils.h"

#include <algorithm>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cmath>
#include <ctime>

#ifdef __SSE3__
#include <x86intrin.h>
#endif

using namespace std;

namespace foreseqs {

BLMode branchLengthsMode = BL_AVERAGE;
CatMode categoriesMode = CAT_ESTIMATE;
PredMode predictionMode = PRED_ANCSEQ;

/* Number of categories hardcoded to 4 */
unsigned int numberOfRateCategories = 4;
unsigned int numberOfTaxa, sequenceLength;
unsigned int numberOfThreads = 1;
double threshold = 0;

char ** taxaNames;
unsigned int * seqIndexTranslate;

bool predictSequences = true;

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

void Utils::transpose(double * m, size_t size)
{
    for (size_t i = 0; i < size; i++) {
        for (size_t j = i + 1; j < size; j++) {
            std::swap(m[i*size + j], m[j*size + i]);
        }
    }
}

void Utils::matrixMultiply(size_t ncols, size_t nrows, const double * A,
		double * B, double * result) {
	assert(ncols == 4 || ncols == 20);
		transpose(B, ncols);
		for (size_t i = 0; i < nrows; i++) {
			for (size_t j = 0; j < ncols; j++) {
#ifdef __SSE3__
				__m128d c = _mm_setzero_pd();

				for (size_t k = 0; k < ncols; k += 2) {
					c = _mm_add_pd(c, _mm_mul_pd(_mm_load_pd(&A[i*ncols + k]), _mm_load_pd(&B[j*ncols + k])));
				}
				c = _mm_hadd_pd(c, c);
				_mm_store_sd(&result[i*ncols + j], c);
#else
				double c = 0;
				for (size_t k = 0; k < ncols; k++) {
					c += A[i * ncols + k] * B[j * ncols + k];
				}
				result[i * ncols + j] = c;
#endif
			}
		}
		transpose(B, ncols);
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

vector< vector<unsigned int> > Utils::findMissingSequences( pll_partition_t ** partitions, size_t numberOfPartitions ) {

	vector< vector<unsigned int> > missingSeqs(numberOfPartitions);

	for (size_t part = 0; part < numberOfPartitions; part++) {

		pll_partition_t * partition = partitions[part];
		unsigned char undefinedSite = (partition->states==4)?15:22;

		int start = 0,
			  end = partition->sites;
		int count = end - start;

		int defined_sites;
		for (unsigned int i = 1; i <= (unsigned int) partition->tips; i++) {
			defined_sites = 0;
			for (int j = start; j < end; j++) {
				if (partition->tipchars[i][j] != undefinedSite)
					defined_sites++;
			}
			if (defined_sites <= (threshold * count)) {
				cout << "Add sequence " << i << " " << defined_sites << " " << count << " " << threshold * count << endl;
				missingSeqs.at(part).push_back(i);
			}
		}
	}

	return missingSeqs;
}

bool Utils::subtreeIsMissing( pll_utree_t * node,
	                            vector<unsigned int> * missingSequences,
		  	 										  vector<pll_utree_t *> * missingBranches,
															size_t numberOfTips ) {

	if (find(missingBranches->begin(), missingBranches->end(), node)
			!= missingBranches->end()) {
		return true;
	}
	if (node->clv_index >= numberOfTips) {
		if (subtreeIsMissing(node->next->back, missingSequences,
				missingBranches, numberOfTips)
				& subtreeIsMissing(node->next->next->back, missingSequences,
						missingBranches, numberOfTips)) {
			cout << "Push " << node->clv_index << " and " << node->back->clv_index << endl;
			missingBranches->push_back(node);
			missingBranches->push_back(node->back);
			return true;
		} else {
			return false;
		}
	} else {
		if (find(missingSequences->begin(), missingSequences->end(), node->clv_index)
				!= missingSequences->end()) {

			cout << "Push " << node->clv_index << " and " << node->back->clv_index << endl;
			missingBranches->push_back(node);
			missingBranches->push_back(node->back);

			/* remove visited taxon */
			cout << "Remove taxon " << node->clv_index << endl;
			missingSequences->erase(
					remove(missingSequences->begin(), missingSequences->end(),
							node->clv_index), missingSequences->end());

			return true;
		} else {
			return false;
		}
	}
}

pll_utree_t * Utils::findRootingNode( pll_utree_t * tree,
	                                    vector<unsigned int> * missingSequences,
														      	  vector<pll_utree_t *> * missingBranches,
														      	  size_t numberOfTips ) {

	pll_utree_t * currentNode;

	if (!missingSequences->size()) {
		return(0);
	}

	if (tree->back->clv_index < numberOfTips) {
		cout << "Erase " << tree->back->clv_index << endl;
		missingSequences->erase(
				remove(missingSequences->begin(), missingSequences->end(),
						tree->back->clv_index), missingSequences->end());
	}

	/* start searching in a random missing node */
	currentNode = tree;

	while (true) {
		bool missingRight = subtreeIsMissing(currentNode->next->back, missingSequences, missingBranches, numberOfTips);
		bool missingLeft = subtreeIsMissing(currentNode->next->next->back, missingSequences, missingBranches, numberOfTips);

		if (missingRight && missingLeft) {
			cerr << "ERROR: Everything is missing!!" << endl;
			exit(EX_IOERR);
		} else if (!(missingRight || missingLeft)) {
// #ifdef PRINT_TRACE
			cout << "TRACE: Found ancestor in " << currentNode->clv_index << endl;
// #endif

			missingBranches->push_back(currentNode);
			missingBranches->push_back(currentNode->back);

			return currentNode;
		} else {
			/* move to next position */
			missingBranches->push_back(currentNode);
			missingBranches->push_back(currentNode->back);
			currentNode =
					(!missingRight) ?
							currentNode->next->back :
							currentNode->next->next->back;
// #ifdef PRINT_TRACE
			cout << "TRACE: Moving node to " << currentNode->clv_index << endl;
// #endif
		}
	}

}

static bool compareNodes (pll_utree_t * a, pll_utree_t * b) {
	return (a->clv_index < b->clv_index);
}

vector< vector<pll_utree_t *> > Utils::findMissingBranches (
															  pll_utree_t ** tipNodes,
																vector< vector<unsigned int> > missingSequences,
																size_t numberOfTips,
															  size_t numberOfPartitions ) {

	vector< vector<pll_utree_t *> > missingBranches((size_t)numberOfPartitions);

	for (size_t part = 0; part < numberOfPartitions; part++) {
		cout << " " << missingSequences[part].size() << endl;
		 for (unsigned int ms=0; ms < missingSequences[part].size(); ms++)
		 	cout << " " << missingSequences[part][ms];
		 	cout << endl;
		while (missingSequences[part].size()) {
			cout << " Find rooting for " << tipNodes[part]->node_index << "." << tipNodes[part]->clv_index << " [" << part << "]" << endl;
			findRootingNode(tipNodes[part]->back,
				              &missingSequences[part],
											&missingBranches[part],
											numberOfTips);
		}
		cout << " " << missingSequences[part].size() << endl;
		 for (unsigned int ms=0; ms < missingSequences[part].size(); ms++)
		 	cout << " " << missingSequences[part][ms];
		 	cout << endl;
		sort(missingBranches[part].begin(), missingBranches[part].end(), compareNodes);
	}

	return missingBranches;
}

void * Utils::allocate(size_t n, size_t el_size) {
	void * mem;
#ifdef __SSE3__
	if (posix_memalign(&mem, 16, n * el_size)) {
		cerr << "Error allocating aligned memory" << endl;
		exit(EX_MEMORY);
	}
#else
	mem = malloc(n * el_size);
#endif
	return mem;
}

void Utils::optimizeModelParameters(pll_partition_t ** partitions,
																		pll_utree_t ** trees,
																		size_t numberOfPartitions) {
	// double lk = 0.0;
	// double epsilon = 0.01;
	// bool optimizeBranchLengths = pllPartitions->numberOfPartitions > 1;
	//
	// if (optimizeBranchLengths) {
	// 	/*
	// 	 * Optimize per-gene branch lengths.
	// 	 */
	// 	cout << "Optimizing per-gene branch lengths / model parameters " << endl;
	//
	// 	int smoothIterations = 64;
	// 	do {
	// 		lk = pllTree->likelihood;
	// 		pllOptimizeBranchLengths(pllTree, pllPartitions, smoothIterations);
	// 		pllOptimizeModelParameters(pllTree, pllPartitions, 0.1);
	// 	} while (fabs(lk - pllTree->likelihood) > epsilon);
	// } else {
	// 	/*
	// 	 * In case there is one single partition, we do not optimize the branch lengths.
	// 	 * Otherwise we would have weird results in the branches with missing data.
	// 	 */
	// 	cout << "Optimizing model parameters " << endl;
	//
	// 	pllOptRatesGeneric(pllTree, pllPartitions, 1.0,
	// 			pllPartitions->rateList);
	// 	pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
	// 			false);
	// 	pllOptBaseFreqs(pllTree, pllPartitions, 1.0, pllPartitions->freqList);
	// 	pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
	// 			false);
	// 	pllOptAlphasGeneric(pllTree, pllPartitions, 1.0,
	// 			pllPartitions->alphaList);
	// 	pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
	// 			false);
	// 	do {
	// 		lk = pllTree->likelihood;
	// 		pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
	// 				false);
	// 		pllOptRatesGeneric(pllTree, pllPartitions, 0.1,
	// 				pllPartitions->rateList);
	// 		pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
	// 				false);
	// 		pllOptBaseFreqs(pllTree, pllPartitions, 0.1,
	// 				pllPartitions->freqList);
	// 		pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
	// 				false);
	// 		pllOptAlphasGeneric(pllTree, pllPartitions, 0.1,
	// 				pllPartitions->alphaList);
	// 		pllEvaluateLikelihood(pllTree, pllPartitions, pllTree->start, true,
	// 				false);
	// 	} while (fabs(lk - pllTree->likelihood) > epsilon);
	// }
}

} /* namespace foreseqs */
