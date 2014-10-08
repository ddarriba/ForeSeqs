/*
 * DnaModel.cpp
 *
 *  Created on: Oct 7, 2014
 *      Author: diego
 */

#include "DnaModel.h"
#include "Utils.h"

#include <iostream>
#include <cstring>
#include <cassert>
#include <cmath>

using namespace std;

namespace seqpred {

DnaModel::DnaModel(partitionList * pllPartitions, int partitionIndex) :
		Model(pllPartitions, partitionIndex) {

	assert (partitionInfo->states == NUM_NUC);
	assert (numberOfStates == NUM_NUC);

	int numFreqs = NUM_NUC;
	int numRates = (NUM_NUC - 1) * NUM_NUC / 2;
	frequencies.resize(numFreqs);
	substRates.resize(numRates);
	memcpy(&(frequencies[0]), partitionInfo->frequencies, numFreqs * sizeof(double));
	memcpy(&(substRates[0]), partitionInfo->substRates, numRates * sizeof(double));

	SetupGTR();
}

double computeFracchange(vector<double> freqs, vector<double> substRates) {
	/* convert rates into matrix */
	double r[NUM_NUC][NUM_NUC];
	int i = 0;
	for (int j = 0; j < (NUM_NUC-1); j++)
		for (int k = j + 1; k < NUM_NUC; k++)
			r[j][k] = substRates[i++];
	for (int j = 0; j < NUM_NUC; j++) {
		r[j][j] = 0.0;
		for (int k = 0; k < j; k++)
			r[j][k] = r[k][j];
	}
	/* evaluate fracchange */
	double fracchange = 0.0;
	for (int j = 0; j < NUM_NUC; j++)
		for (int k = 0; k < NUM_NUC; k++)
			fracchange += freqs[j] * r[j][k] * freqs[k];
	return fracchange;
}

char DnaModel::getState(double * P) const {
	int j;
	double r = Utils::genRand();
	for (j=0; r>(*P) && j<NUM_NUC-1; j++) P++;
	return (states[pow(2,j)]);
}

void DnaModel::SetupGTR() {

	double fracchange = computeFracchange(frequencies, substRates);
	for (int i = 0; i < NUM_NUC; i++) {
		Root[i] = -partitionInfo->EIGN[i] / fracchange;
	}

	for (int i = 0; i < NUM_NUC; i++) {
		for (int j = 0; j < NUM_NUC; j++) {
			for (int k = 0; k < NUM_NUC; k++) {
				Cijk[i * SQNUM_NUC + j * NUM_NUC + k] = partitionInfo->EI[i * NUM_NUC + k]
										* partitionInfo->EV[j * NUM_NUC + k];
			}
		}
	}

}

void DnaModel::setMatrix(double * matrix, double branchLength) const {
	int i, j, k;
	double expt[NUM_NUC];
	double *P;

	P = matrix;
	if (branchLength < 1e-6) {
		for (i = 0; i < NUM_NUC; i++) {
			for (j = 0; j < NUM_NUC; j++) {
				if (i == j)
					*P = 1.0;
				else
					*P = 0.0;
				P++;
			}
		}
		return;
	}

	for (k = 1; k < NUM_NUC; k++) {
		expt[k] = exp(branchLength * Root[k]);
	}
	for (i = 0; i < NUM_NUC; i++) {
		for (j = 0; j < NUM_NUC; j++) {
			(*P) = Cijk[i * SQNUM_NUC + j * NUM_NUC + 0];
			for (k = 1; k < NUM_NUC; k++) {
				(*P) += Cijk[i * SQNUM_NUC + j * NUM_NUC + k] * expt[k];
			}
			P++;
		}
	}

	/* the rows are cumulative to help with picking one using
	 a random number */
	for (int i = 0; i < NUM_NUC; i++) {
		for (int j = 1; j < NUM_NUC; j++) {
			int nextIndex = NUM_NUC * i + j;
			matrix[nextIndex] += matrix[nextIndex - 1];
		}
		assert(Utils::floatEquals(matrix[NUM_NUC * (i + 1) - 1], 1.0));
	}

}

DnaModel::~DnaModel() {
	// TODO Auto-generated destructor stub
}

} /* namespace seqpred */
