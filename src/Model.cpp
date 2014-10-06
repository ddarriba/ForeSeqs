/*
 * Model.cpp
 *
 *  Created on: Oct 1, 2014
 *      Author: diego
 */

#include "Model.h"
#include "Utils.h"

#include <iostream>
#include <cstring>
#include <cassert>
#include <cmath>

using namespace std;

namespace seqpred {

#define SQNUM_NUC 4*4
#define NUM_NUC 4
#define CUNUM_NUC 256
#define pos(i,j,n)      ((i)*(n)+(j))

#define BASE        2    /* base of floating point arithmetic */
#define DIGITS     53    /* no. of digits to the base BASE in the fraction */
#define MAXITER    30    /* max. no. of iterations to converge */

typedef struct {
	double re, im;
} complex;

static double Cijk[CUNUM_NUC], Root[NUM_NUC];

void Model::SetupGTR() {
	int i, j, k;
	double mr;
	double sum;
	double Qij[SQNUM_NUC];
	double U[SQNUM_NUC], V[SQNUM_NUC], T1[SQNUM_NUC], T2[SQNUM_NUC];

	k = 0;
	for (i = 0; i < NUM_NUC - 1; i++) {
		for (j = i + 1; j < NUM_NUC; j++) {
			Qij[i * NUM_NUC + j] = Qij[j * NUM_NUC + i] = substRates[k++];
		}
	}

	for (i = 0; i < NUM_NUC; i++) {
		for (j = 0; j < NUM_NUC; j++) {
			Qij[i * NUM_NUC + j] *= frequencies[j];
		}
	}

	mr = 0;
	for (i = 0; i < NUM_NUC; i++) {
		sum = 0;
		Qij[i * NUM_NUC + i] = 0;
		for (j = 0; j < NUM_NUC; j++) {
			sum += Qij[i * NUM_NUC + j];
		}
		Qij[i * NUM_NUC + i] = -sum;
		mr += frequencies[i] * sum;
	}

	Utils::abyx(1.0 / mr, Qij, SQNUM_NUC);

	if ((k = Utils::eigen(1, Qij, NUM_NUC, Root, T1, U, V, T2)) != 0) {
		fprintf(stderr, "\ncomplex roots in SetupGTR");
		assert(0);
	}

	Utils::xtoy(U, V, SQNUM_NUC);
	Utils::matinv(V, NUM_NUC, NUM_NUC, T1);

	for (i = 0; i < NUM_NUC; i++) {
		for (j = 0; j < NUM_NUC; j++) {
			for (k = 0; k < NUM_NUC; k++) {
				Cijk[i * SQNUM_NUC + j * NUM_NUC + k] = U[i * NUM_NUC + k]
						* V[k * NUM_NUC + j];
			}
		}
	}

	/**
	 * TODO: OUTPUT FOR THIS FUNCTION ARE Cijk AND Root
	 * 	     CHECK WHETHER THIS IS NECESSARY
	 */
//	cout << "Cijk: ";
//	Utils::printVector(Cijk, CUNUM_NUC);
//	cout << "Root ";
//	Utils::printVector(Root, NUM_NUC);
//	Utils::printVector(partitionInfo->EIGN, NUM_NUC);
//	cout << "U, V ";
//	Utils::printVector(U, SQNUM_NUC);
//	Utils::printVector(V, SQNUM_NUC);
//	cout << "EI, EV ";
//		Utils::printVector(partitionInfo->EI, SQNUM_NUC);
//		Utils::printVector(partitionInfo->EV, SQNUM_NUC);
}

void Model::setMatrix(double * matrix, double branchLength) {
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

Model::Model(partitionList * pllPartitions, int partitionIndex) :
		partitionInfo(pllPartitions->partitionData[partitionIndex]) {
	pInfo * part = pllPartitions->partitionData[partitionIndex];

	int numFreqs = part->states;
	int numRates = (part->states - 1) * part->states / 2;
	frequencies.resize(numFreqs);
	substRates.resize(numRates);
	memcpy(&(frequencies[0]), part->frequencies, numFreqs * sizeof(double));
	memcpy(&(substRates[0]), part->substRates, numRates * sizeof(double));

	SetupGTR();
}

Model::~Model() {
	/* nothing to do */
}

} /* namespace seqpred */
