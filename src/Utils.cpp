/*
 * Utils.cpp
 *
 *  Created on: Oct 2, 2014
 *      Author: Diego Darriba
 *      E-mail: diego.darriba@h-its.org
 *
 *  This file is part of SeqPred.
 *
 *  SeqPred is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SeqPred is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with SeqPred.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Utils.h"

#include <algorithm>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cmath>

#ifdef __SSE3__
#include <x86intrin.h>
#endif

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

DataType Utils::getDataType(const partitionList * pllPartitions, size_t numberOfPartition) {
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

vector< vector<unsigned int> > Utils::findMissingSequences( pllInstance * pllTree, partitionList * pllPartitions ) {

	vector< vector<unsigned int> > missingSeqs((size_t)pllPartitions->numberOfPartitions);

	for (size_t part = 0; part < (size_t)pllPartitions->numberOfPartitions; part++) {

		unsigned char undefinedSite = (pllPartitions->partitionData[part]->states==4)?15:22;

		int start = pllPartitions->partitionData[part]->lower,
			  end = pllPartitions->partitionData[part]->upper;

		int missing;
		for (unsigned int i = 1; i <= (unsigned int) pllTree->mxtips; i++) {
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

boolean Utils::subtreeIsMissing( pllInstance * pllTree, vector<unsigned int> * missingSequences, 
	const nodeptr node, vector<nodeptr> * missingBranches ) {

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

nodeptr Utils::findRootingNode( pllInstance * pllTree, vector<unsigned int> * missingSequences, 
	nodeptr startingNode, vector<nodeptr> * missingBranches ) {

	nodeptr currentNode;

	if (!missingSequences->size()) {
		return(0);
	}

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
#ifdef PRINT_TRACE
			cout << "TRACE: Moving node to " << currentNode->number << endl;
#endif
		}
	}

}

static bool compareNodes (nodeptr a, nodeptr b) { return (a->number < b->number); }

std::vector< std::vector<nodeptr> > Utils::findMissingBranches ( pllInstance * pllTree, partitionList * pllPartitions, 
	vector< vector<unsigned int> > missingSequences ) {

	vector< vector<nodeptr> > missingBranches((size_t)pllPartitions->numberOfPartitions);

	for (size_t part = 0; part < (size_t)pllPartitions->numberOfPartitions; part++) {
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

} /* namespace seqpred */
