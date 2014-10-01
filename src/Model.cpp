/*
 * Model.cpp
 *
 *  Created on: Oct 1, 2014
 *      Author: diego
 */

#include "Model.h"
#include <cstring>

namespace seqgen {

void Model::setMatrix(double * matrix, double branchLength)
{
        int i,j,k;
        double expt[4];
        double *P;

/* P(t)ij = SUM Cijk * exp{Root*t}
*/
        P=matrix;
        if (branchLength<1e-6) {
                for (i=0; i<4; i++) {
                        for (j=0; j<4; j++) {
                                if (i==j)
                                        *P=1.0;
                                else
                                        *P=0.0;
                                P++;
                        }
                }
                return;
        }

//        for (k=1; k<4; k++) {
////                expt[k]=exp(branchLength*Root[k]);
//                expt[k]=exp(branchLength*partitionInfo->EIGN[k]);
//        }
//        for (i=0; i<4; i++) {
//                for (j=0; j<4; j++) {
//                        (*P)=Cijk[i*4*4+j*4+0];
//                        for (k=1; k<4; k++) {
//                                (*P)+=Cijk[i*4*4+j*4+k]*expt[k];
//                        }
//                        P++;
//                }
//        }
//
//        CumulativeRows(matrix);
}


Model::Model(partitionList * pllPartitions, int partitionIndex) :
		partitionInfo(pllPartitions->partitionData[partitionIndex]) {
	pInfo * part = pllPartitions->partitionData[partitionIndex];

	int numFreqs = part->states;
	int numRates = (part->states-1)*part->states / 2;
	frequencies.resize(numFreqs);
	substRates.resize(numRates);
	memcpy(&(frequencies[0]), part->frequencies, numFreqs * sizeof(double));
	memcpy(&(substRates[0]), part->frequencies, numRates * sizeof(double));
}

Model::~Model() {
	/* nothing to do */
}

} /* namespace seqgen */
