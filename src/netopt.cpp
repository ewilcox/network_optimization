//============================================================================
// Name        : netopt.cpp
// Author      : Eric Wilcox
// Version     :
// Copyright   : Analysis of Algorithms Course Project
// Description : Network Optimization project w/ graph generation and several
//			   		mixed routing algorithms and tests.
//============================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include "my_rand.hpp"

using namespace std;

const int MAXVERTICES = 10;  	// will be 5000
const int SPARSECONNECT = 6;	// number of sparse connections

void printGraph(vector <vector <int> > G) {
	unsigned int i, j;
	for (i=0; i<MAXVERTICES; i++) {
		cout << "v[" << i << "] = ";
		for (j=0; j<G.at(i).size(); j++)
			cout << G.at(i).at(j) << "  ";
		cout << endl;
	}
}
void makeSparseGraph(vector <vector <int> > &G) {
	for (int i=0; i<MAXVERTICES; i++) {
		for (int j=0; j<SPARSECONNECT; j++) {
			//G.at(i).push_back(getRand(0,MAXVERTICES));
		}
	}
}
void makeDenseGraph(vector <vector <int> > &G) {
	// TODO: make dense graph, return structure?
}

int main() {
	vector<vector <int> > G1(MAXVERTICES);
	vector<vector <int> > G2(MAXVERTICES);
	makeSparseGraph(G1);
	printGraph(G1);

	getRand(0,100);

	return 0;
}
