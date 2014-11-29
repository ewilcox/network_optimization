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

const int MAXVERTICES = 104;
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
	printline('-',MAXVERTICES);
}
// Print graph in matrix format for checking (row,col)
void printAsMatrix(vector<vertex> G) {
	u_int col, row = 0;
	printline('-',MAXVERTICES*10);
	int graph[MAXVERTICES][MAXVERTICES] = {0};
	for (row=0; row<MAXVERTICES; ++row) {
		for (auto&& v : G) {
			for (auto &col : v.edges) {
				graph[row][col.connection] = col.weight;
			}
			++row;
		}
	}
	cout << "      ";
	for (col=0; col<MAXVERTICES; ++col) cout << setw(5) << col;
	cout << endl;
	printline('-',MAXVERTICES*7);
	for (row=0; row<MAXVERTICES; ++row) {
		if (row>=100) cout << row << " - ";
		else if (row>=10) cout << " " << row << " - ";
		else cout << "  " << row << " - ";
		for (col=0; col<MAXVERTICES; ++col) {
			cout << setw(5) << graph[row][col];
		}
		cout << endl;
	}
	printline('-',MAXVERTICES*10);
}
// Finds duplicate number in the vector, return true if already exists in vector
bool duplicate(vector<edge> e, int vertex, edge newedge) {
	if (vertex == newedge.connection) return true;		// prevents cycle to same vertex
	for (u_int i=0; i<e.size(); ++i)
		if (e.at(i).connection == newedge.connection) return true;
	return false;
}
void printEdge(vector<edge> e) {
	for (u_int i=0; i<e.size(); ++i) cout << e.at(i).connection << "  ";
}
// Find matching edge and return index or -1 if not in edges list
int my_find(vector<edge> e, int match) {
	for (u_int i=0; i<e.size(); ++i) {
		if (match == e.at(i).connection) return i;
	}
	return -1;
}
// Makes sparse graph with number of connections per node
// No loops on 1 node, undirected (each edge = new edge in that node) and occasionally
// reached stale loop at end, so remove random element from start and replace
// This method creates column first order, but can end up short on total matrix
void makeGraph(vector<vertex> &G, u_int connections) {
	edge newedge1, newedge2;
	int badLoop = 0, temp;
	for (u_int j=0; j<connections; ++j) {
		for (int i=0; i<MAXVERTICES; i++) {
			if (G.at(i).edges.size() >= connections) continue;
			newedge1.connection = getRand(0,MAXVERTICES-1);
			newedge1.weight = getRand(1,MAXWEIGHT);
			newedge2.connection = i;
			newedge2.weight = newedge1.weight;
			while (duplicate(G.at(i).edges, i, newedge1) ||
					G.at(newedge1.connection).edges.size() >= connections) {		// keep random weight, get new connect
				newedge1.connection = getRand(0,MAXVERTICES-1);
				newedge2.connection = i;
				++badLoop;
				if (badLoop >= 2*MAXVERTICES) {
					temp = G.at(newedge1.connection).edges.front().connection;
					G.at(newedge1.connection).edges.erase(G[newedge1.connection].edges.begin());
					G[temp].edges.erase(G[temp].edges.begin()+my_find(G.at(temp).edges, newedge1.connection));
					--j;	// decrement outside loop per set delete to ensure enough add's to graph
				}
			}
			badLoop=0;
			G.at(i).edges.push_back(newedge1);
			G.at(newedge1.connection).edges.push_back(newedge2);
		}
	}
}

// Swapping function vertex
void swapem(vertex &a, vertex &b) {
	vertex temp;
	temp = a;
	a = b;
	b = temp;
}

// MINMUM function, returns top level element of the heap
vertex minHeap(vector<vertex> heap) {
	if (heap.size() > 0) return heap[0];
	else {
		cout << "Error - heap is empty when returning minHeap()\n";
		vertex error;
		error.value = -1;
		return error;
	}
}

void heapify(vector<vertex> &heap) {
	int current, parent;
	current = heap.size()-1;	// set current to last node (end of array)
	while (current >= 0) {
		parent = (current-1)/2;		// set parent to root of current
		if (heap[current].value < heap[parent].value) swapem(heap[current],heap[parent]);
		--current;					// cycle down in nodes until reach root
	}
}

// Insert element in correct place in heap so don't need sort
// vector takes care of size so don't need to maintain count
void insertHeap(vector<vertex> &heap, vertex v) {
	heap.push_back(v);
	heapify(heap);
}

void deleteHeap(vector<vertex> heap, vertex v) {
	heap[0] = heap.back();
//	heapsort(heap);
}

void printHeap(vector<vertex> heap) {
	for (u_int i=0; i<heap.size(); ++i) {
		cout << "H[" << i << "]=" << heap[i].value << "      ";
	}
	cout << endl;
}
int main() {
	vector<vertex> G1(MAXVERTICES);
	vector<vertex> G2(MAXVERTICES);
	makeGraph(G1,6);					// make sparse graph
	makeGraph(G2,MAXVERTICES*0.20);		// make dense graph
//	printGraph(G1);						// debug statements for graph and matrix views
//	printAsMatrix(G1);
//	printGraph(G2);
//	printAsMatrix(G2);					// to here.

	vertex v;							// create new struct for heap creation and testing.
	vector<vertex> heap;
	for (int i=0; i<20; ++i) {
		v.value = i;
		v.value = getRand(1,MAXWEIGHT);
		insertHeap(heap, v);
		printHeap(heap);
		printline('-',20*heap.size());
	}

	return 0;
}
