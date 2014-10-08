/*
 * Utils.cpp
 *
 *  Created on: Oct 2, 2014
 *      Author: diego
 */

#include "Utils.h"

#include <climits>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "pll.h"

using namespace std;

namespace seqpred {

DataType dataType = DT_NUCLEIC;

int numberOfStates;
std::vector<char> states;
std::map<char, int> statesMap;
bool initialized = false;

int numberOfTaxa, sequenceLength;
char ** taxaNames;

map<string, int> protModelsMap;
int protModel = -1;

void Utils::init() {

	if (initialized)
		return;

	switch (dataType) {
	case DT_NUCLEIC:
		numberOfStates = 4;
		states.resize(16);
		states[1] = 'A'; //1; //A
		states[2] = 'C'; //2; //C
		states[4] = 'G'; //4; //G
		states[8] = 'T'; //8; //T
		statesMap['a'] = statesMap['A'] = 0;
		statesMap['c'] = statesMap['C'] = 1;
		statesMap['g'] = statesMap['G'] = 2;
		statesMap['t'] = statesMap['T'] = 3;
		break;
	case DT_PROTEIC:
		protModelsMap["DAYHOFF"]  = PLL_DAYHOFF;
		protModelsMap["DCMUT"] 	  = PLL_DCMUT;
		protModelsMap["JTT"]      = PLL_JTT;
		protModelsMap["JTTDCMUT"] = PLL_JTTDCMUT;
		protModelsMap["MTREV"]    = PLL_MTREV;
		protModelsMap["WAG"]      = PLL_WAG;
		protModelsMap["RTREV"]    = PLL_RTREV;
		protModelsMap["CPREV"]    = PLL_CPREV;
		protModelsMap["VT"]       = PLL_VT;
		protModelsMap["BLOSUM62"] = PLL_BLOSUM62;
		protModelsMap["MTMAM"]    = PLL_MTMAM;
		protModelsMap["LG"]       = PLL_LG;
		protModelsMap["MTART"]    = PLL_MTART;
		protModelsMap["MTZOA"]    = PLL_MTZOA;
		protModelsMap["PMB"]      = PLL_PMB;
		protModelsMap["HIVB"]     = PLL_HIVB;
		protModelsMap["HIVW"]     = PLL_HIVW;
		protModelsMap["FLU"]      = PLL_FLU;
		break;
	default:
		cerr << "ERROR: Unimplemented data type" << endl;
		exit(EX_UNIMPLEMENTED);
	}
	initialized = true;
}

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

} /* namespace seqpred */
