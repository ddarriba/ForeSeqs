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

#define pos(i,j,n)      ((i)*(n)+(j))

#define BASE        2    /* base of floating point arithmetic */
#define DIGITS     53    /* no. of digits to the base BASE in the fraction */
#define MAXITER    30    /* max. no. of iterations to converge */

typedef struct {
	double re, im;
} complex;

//static double Cijk[CUNUM_NUC], Root[NUM_NUC];

DnaModel::DnaModel(partitionList * pllPartitions, int partitionIndex) :
		Model(pllPartitions, partitionIndex) {

	int numFreqs = partitionInfo->states;
	int numRates = (partitionInfo->states - 1) * partitionInfo->states / 2;
	frequencies.resize(numFreqs);
	substRates.resize(numRates);
	memcpy(&(frequencies[0]), partitionInfo->frequencies, numFreqs * sizeof(double));
	memcpy(&(substRates[0]), partitionInfo->substRates, numRates * sizeof(double));

	SetupGTR();
}

double computeFracchange(double * freqs, double * substRates) {
	/* convert rates into matrix */
	double r[4][4];
	int i = 0;
	for (int j = 0; j < 3; j++)
		for (int k = j + 1; k < 4; k++)
			r[j][k] = substRates[i++];
	for (int j = 0; j < 4; j++) {
		r[j][j] = 0.0;
		for (int k = 0; k < j; k++)
			r[j][k] = r[k][j];
	}
	/* evaluate fracchange */
	double fracchange = 0.0;
	for (int j = 0; j < 4; j++)
		for (int k = 0; k < 4; k++)
			fracchange += freqs[j] * r[j][k] * freqs[k];
	return fracchange;
}

void DnaModel::SetupGTR() {

	double fracchange = computeFracchange(partitionInfo->frequencies, partitionInfo->substRates);
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
	double expt[numberOfStates];
	double *P;

	P = matrix;
	if (branchLength < 1e-6) {
		for (i = 0; i < numberOfStates; i++) {
			for (j = 0; j < numberOfStates; j++) {
				if (i == j)
					*P = 1.0;
				else
					*P = 0.0;
				P++;
			}
		}
		return;
	}

	for (k = 1; k < numberOfStates; k++) {
		expt[k] = exp(branchLength * Root[k]);
	}
	for (i = 0; i < numberOfStates; i++) {
		for (j = 0; j < numberOfStates; j++) {
			(*P) = Cijk[i * numberOfStates * numberOfStates + j * numberOfStates + 0];
			for (k = 1; k < numberOfStates; k++) {
				(*P) += Cijk[i * numberOfStates * numberOfStates + j * numberOfStates + k] * expt[k];
			}
			P++;
		}
	}

	/* the rows are cumulative to help with picking one using
	 a random number */
	for (int i = 0; i < numberOfStates; i++) {
		for (int j = 1; j < numberOfStates; j++) {
			int nextIndex = numberOfStates * i + j;
			matrix[nextIndex] += matrix[nextIndex - 1];
		}
		assert(Utils::floatEquals(matrix[numberOfStates * (i + 1) - 1], 1.0));
	}

}

DnaModel::~DnaModel() {
	// TODO Auto-generated destructor stub
}

} /* namespace seqpred */
