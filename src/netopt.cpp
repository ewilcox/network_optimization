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
using namespace std;

const int maxVertices = 10;  // will be 5000

void printGraph(vector <vector <int> > G) {
	unsigned int i, j;
	for (i=0; i<maxVertices; i++) {
		cout << "v[" << i << "] = ";
		for (j=0; j<G.at(i).size(); j++)
			cout << G.at(i).at(j) << ", ";
		cout << endl;
	}
}
void makeSparseGraph(vector <vector <int> > G) {

	for (int i=0; i<maxVertices; i++) {
		G.at(i).push_back(i*5);
		G.at(i).push_back(i*15);
	}
	G.at(3).push_back(175);
	G.at(3).push_back(215);
	G.at(3).push_back(225);
	G.at(3).push_back(25);
	printGraph(G);
}
void makeDenseGraph() {
	// TODO: make dense graph, return structure?
}

int main() {
	vector<vector <int> > G1(maxVertices);
	vector<vector <int> > G2(maxVertices);
	makeSparseGraph(G1);

	// want adjacency list or matrix for graph?
	// vectors?  need random number generation - will need a good one.

	return 0;
}
