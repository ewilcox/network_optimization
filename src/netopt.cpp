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

const int MAXVERTICES = 9;  	// will be 5000 (4999 due to starting at 0)
const int SPARSECONNECT = 6;	// number of sparse connections
const int MAXWEIGHT = 1000;		// max weight of graph edges

struct edges {
	int connection;
	int weight;
};
void printline(char c, int num) {
	for (int i=0; i<num; ++i) cout << c;
	cout << endl;
}
void printGraph(vector <vector <edges> > G) {
	unsigned int i, j;
	for (i=0; i<MAXVERTICES; i++) {
		cout << "v[" << i << "] = ";
		for (j=0; j<G.at(i).size(); j++) {
			cout<<'('<<G.at(i).at(j).connection<<'|'<<G.at(i).at(j).weight << ")  ";
		}
		cout << endl;
	}
}
// Finds duplicate number in the vector, return true if already exists in vector
bool duplicate(vector<edges> v, int vertex, int num) {
	if (vertex == num) return true;		// prevents cycle to same vertex
	for (u_int i=0; i<v.size(); ++i)
		if (v.at(i).connection == num) return true;
	return false;
}
// Makes sparse graph with number of connections per node = const SPARSECONNECT
void makeSparseGraph(vector <vector <edges> > &G) {
	edges e;
	for (int i=0; i<MAXVERTICES; i++) {
		for (int j=0; j<SPARSECONNECT; j++) {
			e.connection = getRand(0,MAXVERTICES);
			e.weight = getRand(1,MAXWEIGHT);
			while (duplicate(G.at(i),i,e.connection)) e.connection = getRand(0,MAXVERTICES);
			G.at(i).push_back(e);
		}
	}
}
void makeDenseGraph(vector <vector <edges> > &G) {
	edges e;
	for (int i=0; i<MAXVERTICES; i++) {
		for (int j=0; j<MAXVERTICES*0.2; j++) {
			e.connection = getRand(0,MAXVERTICES);
			e.weight = getRand(1,MAXWEIGHT);
			while (duplicate(G.at(i),i,e.connection)) e.connection = getRand(0,MAXVERTICES);
			G.at(i).push_back(e);
		}
	}
}

int main() {
	vector<vector <edges> > G1(MAXVERTICES);
	vector<vector <edges> > G2(MAXVERTICES);
	makeSparseGraph(G1);
	makeDenseGraph(G2);
	printGraph(G1);
	printline('-',100);
	printGraph(G2);

	return 0;
}
