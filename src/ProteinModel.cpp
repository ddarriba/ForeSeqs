/*
 * ProteinModel.cpp
 *
 *  Created on: Oct 8, 2014
 *      Author: diego
 */

#include "ProteinModel.h"
#include "Utils.h"

#include <cassert>
#include <iostream>
#include <cstring>

/* for exit() */
#include <cstdlib>

using namespace std;

namespace seqpred {

ProteinModel::ProteinModel(partitionList * pllPartitions, int partitionIndex) :
				Model(pllPartitions, partitionIndex) {

			assert (partitionInfo->states == NUM_AA);
			assert (numberOfStates == NUM_AA);

			int numFreqs = NUM_AA;
			int numRates = (NUM_AA - 1) * NUM_AA / 2;
			frequencies.resize(numFreqs);
			substRates.resize(numRates);
			memcpy(&(frequencies[0]), partitionInfo->frequencies, numFreqs * sizeof(double));
			memcpy(&(substRates[0]), partitionInfo->substRates, numRates * sizeof(double));

			SetupGTR();

			cout << "That's all folks!" << endl;
			exit(1);
}

ProteinModel::~ProteinModel() {
	// TODO Auto-generated destructor stub
}

void ProteinModel::SetupGTR() {

	double fracchange = 1.0; //TODO: computeFracchange(frequencies, substRates);
	for (int i = 0; i < NUM_AA; i++) {
		Root[i] = -partitionInfo->EIGN[i] / fracchange;
	}

	for (int i = 0; i < NUM_AA; i++) {
		for (int j = 0; j < NUM_AA; j++) {
			for (int k = 0; k < NUM_AA; k++) {
				Cijk[i * SQNUM_AA + j * NUM_AA + k] = partitionInfo->EI[i * NUM_AA + k]
										* partitionInfo->EV[j * NUM_AA + k];
			}
		}
	}

}

char ProteinModel::getState(double * P) const {
	return '?';
}

void ProteinModel::setMatrix(double * matrix, double branchLength) const {
//	int i, j, k;
//	double expt[NUM_NUC];
//	double *P;
//
//	P = matrix;
//	if (branchLength < 1e-6) {
//		for (i = 0; i < NUM_NUC; i++) {
//			for (j = 0; j < NUM_NUC; j++) {
//				if (i == j)
//					*P = 1.0;
//				else
//					*P = 0.0;
//				P++;
//			}
//		}
//		return;
//	}
//
//	for (k = 1; k < NUM_NUC; k++) {
//		expt[k] = exp(branchLength * Root[k]);
//	}
//	for (i = 0; i < NUM_NUC; i++) {
//		for (j = 0; j < NUM_NUC; j++) {
//			(*P) = Cijk[i * SQNUM_NUC + j * NUM_NUC + 0];
//			for (k = 1; k < NUM_NUC; k++) {
//				(*P) += Cijk[i * SQNUM_NUC + j * NUM_NUC + k] * expt[k];
//			}
//			P++;
//		}
//	}
//
//	/* the rows are cumulative to help with picking one using
//	 a random number */
//	for (int i = 0; i < NUM_NUC; i++) {
//		for (int j = 1; j < NUM_NUC; j++) {
//			int nextIndex = NUM_NUC * i + j;
//			matrix[nextIndex] += matrix[nextIndex - 1];
//		}
//		assert(Utils::floatEquals(matrix[NUM_NUC * (i + 1) - 1], 1.0));
//	}

}

} /* namespace seqpred */
