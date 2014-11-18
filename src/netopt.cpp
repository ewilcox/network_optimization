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

const int MAXVERTICES = 20;
const int SPARSECONNECT = 6;	// number of sparse connections
const int MAXWEIGHT = 1000;		// max weight of graph edges

struct edge {
	int connection;
	int weight;
};
struct vertex {
	vector <edge> edges;
	int value;
};

//bool isFull(vector<edges> e, int full) {
//	if (e.size() >= full) return true;
//	return false;
//}
void printline(char c, int num) {
	for (int i=0; i<num; ++i) cout << c;
	cout << endl;
}
void printGraph(vector <vertex> G) {
	u_int i, j;
	for (i=0; i<MAXVERTICES; i++) {
		cout << "v[" << i << "] = ";
		for (j=0; j<G.at(i).edges.size(); j++) {
			cout<<'('<<G.at(i).edges.at(j).connection<<'|'<<G.at(i).edges.at(j).weight << ")  ";
		}
		cout << endl;
	}
}
// Finds duplicate number in the vector, return true if already exists in vector
bool duplicate(vector<edge> e, int vertex, edge newedge) {
	if (vertex == newedge.connection) return true;		// prevents cycle to same vertex
	for (u_int i=0; i<e.size(); ++i)
		if (e.at(i).connection == newedge.connection) return true;
	return false;
}
void printEdge(vector<edge> e) {
	for (u_int i=0; i<e.size(); ++i) {
		cout << e.at(i).connection << "  ";
	}
	cout << endl;
}
// Makes sparse graph with number of connections per node = const SPARSECONNECT
void makeSparseGraph(vector <vertex> &G) {
	edge newedge1, newedge2;
	int stop = 0;
	for (int i=0; i<MAXVERTICES; i++) {
		while (G.at(i).edges.size() < SPARSECONNECT && stop <= 2*MAXVERTICES) {
			newedge1.connection = getRand(0,MAXVERTICES-1);		// -1 accounting for 0 for possible connections
			newedge1.weight = getRand(1,MAXWEIGHT);
			newedge2.connection = i;
			newedge2.weight = newedge1.weight;
			if (!duplicate(G.at(i).edges, i, newedge1) && G.at(newedge1.connection).edges.size() < SPARSECONNECT) {
				G.at(i).edges.push_back(newedge1);
				G.at(newedge1.connection).edges.push_back(newedge2);
				stop = 0;
			}
			else ++stop;
			cout <<"v["<<i<<"] (size "<<G.at(i).edges.size()<<"): ";
			for (u_int k=0; k<G.at(i).edges.size(); ++k) {
				cout << G.at(i).edges.at(k).connection << "   ";
			}
			cout << "  Edges of [" << newedge1.connection << "]: ";
			printEdge(G.at(newedge1.connection).edges);
		}
	}
	if (stop >= 2*MAXVERTICES) cout << "***STOP REACHED DURING EXECUTION OF SPARSEGRAPH***" << endl;
}
void makeDenseGraph(vector <vertex> &G) {
	edge newedge1, newedge2;
	for (int i=0; i<=MAXVERTICES; i++) {
		while (G.at(i).edges.size() < SPARSECONNECT) {
			newedge1.connection = getRand(0,MAXVERTICES);
			newedge1.weight = getRand(1,MAXWEIGHT);
			newedge2.connection = i;
			newedge2.weight = newedge1.weight;
			if (!duplicate(G.at(i).edges, i, newedge1) && G.at(newedge1.connection).edges.size() < SPARSECONNECT) {
				G.at(i).edges.push_back(newedge1);
				G.at(newedge1.connection).edges.push_back(newedge2);
			}
		}
	}
	//	edge e;
//	for (int i=0; i<MAXVERTICES; i++) {
//		for (int j=0; j<MAXVERTICES*0.2; j++) {
//			e.connection = getRand(0,MAXVERTICES);
//			e.weight = getRand(1,MAXWEIGHT);
//			while (duplicate(G.at(i),i,e.connection)) e.connection = getRand(0,MAXVERTICES);
//			G.at(i).push_back(e);
//		}
//	}
}

int main() {
	vector<vertex> G1(MAXVERTICES);
	vector<vertex> G2(MAXVERTICES);
	makeSparseGraph(G1);
//	makeDenseGraph(G2);
	printGraph(G1);
	printline('-',100);
//	printGraph(G2);

	return 0;
}
