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

const int MAXVERTICES = 10;
const int SPARSECONNECT = 6;	// number of sparse connections
const int MAXWEIGHT = 1000;		// max weight of graph edges

enum {unseen, intree, fringe};

struct edge {
	int connection;
	int weight;
};
struct vertex {
	vector <edge> edges;
	int value;
	int status;
};
struct connect {		// struct used for search space exploration
	int from;
	int index;
	int to;
	int weight;
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
	for (auto v : e) cout << v.connection << ":" << v.weight << "  ";
	cout << endl;
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
// Swapping function edge
void swapem(edge &a, edge &b) {
	edge temp;
	a = b;
	b = temp;
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
// Method to print heap out in tree structure for debugging & evaluation
// This function not part of homework - taken directly from online site
// at http://xoax.net/comp_sci/crs/algorithms/lessons/Lesson9/
void printTree(vector<vertex> heap) {
	// Find the largest power of two, That is the depth
		int iDepth = 0;
		int iCopy = heap.size();
		while (iCopy > 0) {
			iCopy >>= 1;
			++iDepth;
		}
		int iMaxWidth = (1 << iDepth);
		int iCharWidth = 4*iMaxWidth;
		u_int iEntry = 0;
		for (int i = 0; i < iDepth; ++i) {
			int iPowerOf2 = (1 << i);
			for (int j = 0; j < iPowerOf2; ++j) {
				int iSpacesBefore = ((iCharWidth/(1 << (i + 1))) - 1);
				 // Spaces before number
				for (int k = 0; k < iSpacesBefore; ++k) {
					cout << " ";
				}
				 // Output an extra space if the number is less than 10
				if (iEntry < 10) {
					cout << " ";
				}
				 // Output the entry and the spaces after it
				cout << heap[iEntry].value;
				++iEntry;
				if (iEntry >= heap.size()) {
					cout << endl;
					return;
				}
				for (int k = 0; k < iSpacesBefore; ++k) {
					cout << " ";
				}
			}
			cout << endl << endl;
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
// Now sort up - for heapsort function
void heapifyUp(vector<vertex> &heap, int start) {
	int current, parent;
	current = heap.size()-1;	// set current to last node (end of array)
	while (current >= start) {
		parent = (current-1)/2;		// set parent to root of current
		if (heap[current].value < heap[parent].value) swapem(heap[current],heap[parent]);
		--current;					// cycle down in nodes until reach start
	}
}
// Now sort down - called from deleteHeap with index
void heapifyDown(vector<vertex> &heap, int i) {
	int parent, end, child;
	parent = i;
	end = heap.size()-1;
	while (parent * 2 + 1 <= end) {
		child = parent * 2 + 1;
		if (child+1 <= end && heap[child].value > heap[child+1].value) ++child;
		if (heap[parent].value > heap[child].value) {
			swapem(heap[parent],heap[child]);
			parent = child;
		}
		else return;//parent = child;
	}
}
// Insert element in correct place in heap so don't need sort
// vector takes care of size so don't need to maintain count
void insertHeap(vector<vertex> &heap, vertex v) {
	heap.push_back(v);
	heapify(heap);
}
void deleteHeap(vector<vertex> &heap, int idx) {
	if (heap.size() <= 0) {
		cout << "Error, heap empty on deleteHeap() call!!\n";
		return;
	}
	swapem(heap[idx],heap[heap.size()-1]);
	heap.pop_back();
	heapifyDown(heap, idx);
}
void printHeap(vector<vertex> heap) {
	for (u_int i=0; i<heap.size(); ++i) {
		cout << "H[" << i << "]=" << heap[i].value << "      ";
	}
	cout << endl;
}
// Heapsort function
vector<vertex> heapsort(vector<vertex> &heap) {
	vector<vertex> sortedHeap;
	heapify(heap);
	cout << "heap[0]: " << heap[0].value << endl;
	int end = heap.size()-1;
	while (end > 0) {
		sortedHeap.push_back(heap[0]);
		swapem(heap[0],heap[end]);		// not using deleteHeap here because we have to re-heapify anyway
		heap.pop_back();				// instead just push to new storage, swap to back and pop.
		heapify(heap);
		--end;
	}
	return sortedHeap;
}
int min(vector<vertex> V) {
	int answer = MAXVERTICES+1;
	for (u_int i=0; i<) {
		if (v.value < answer) answer = v.value;
	}
	return answer;
}
//connect min(vertex v) {
//	connect answer;
//	answer.weight = MAXWEIGHT+1;
//	answer.index = -1;
//	for (u_int i=0; i<v.edges.size(); ++i) {
//		if (v.edges[i].weight < answer.weight) {
//			answer.weight = v.edges[i].weight;
////			current = e.weight;
//			answer.to = v.edges[i].connection;
//			answer.index = i;
//		}
//	}
//	if (answer.index == -1) cout << "Error - returning -1 for edge in min() fct\n";
//	return answer;
//}

//void Dijkstra(vector<vertex> G, int start, int end) {
//	printEdge(G[start].edges);
//	connect n = min(G[start]);
//	cout << "node[" << n.index << "] "<<n.connection<<":"<<n.weight<<endl;
//
//}
//bool fringes(edge e) {
//	for (auto e : v.edges) {
//		if ()
//	}
//}
// Dijkstra's algorithm without heap structure (regular)
void Dijkstra(vector<vertex> G, int start, int end) {
	int i, n = 0;
	connect v;
	connect temp;
	bool done = false;
	int dad[MAXVERTICES] = {-1};
	for (auto v : G) v.status = unseen;
	G[start].status = intree;
	G[start].value = 0;
	dad[start] = NULL;
	for (auto e : G[start].edges) {
		G[e.connection].status = fringe;
		dad[e.connection] = start;
		G[e.connection].value = e.weight;
	}
	vector<connect> fringes;
	for (auto e : G[start].edges) {
		temp.from = start;
		temp.to = e.connection;
		temp.weight = e.weight;
		fringes.push_back(temp);
	}
	while (n < MAXVERTICES) {
		v = min(G[start]);
//		G[start].edges[v.connection].weight = MAXWEIGHT+1;
		cout << "v: " << v.to << "," << v.weight << "  \n";
		G[v.to].status = intree;
		cout << "G["<<v.to<<"] = intree, looping through edges:\n";
		for (i=0; i<G[v.to].edges.size(); ++i) {
			if (G[i].status == unseen) {
				cout << "G["<<i<<"] = unseen, .... = fringe and dad["<<v.to<<"] = "<<start<<":\n";
				G[v.to].status = fringe;
				dad[v.to] = start;
				G[v.to].value = G[start].value + v.weight;
			}
			else if (G[v.to].status == fringe && G[v.to].value < G[start].value + v.weight) {
				cout << "in else if...\n";
				G[v.to].value = G[start].value + v.weight;
				dad[v.to] = start;
			}
			++n;
			cout << "start was " << start << " and is now " << v.to << endl;
			start = v.to;

		}
//		if (v.index == -1) done = true;
//		printEdge(G[start].edges);
//		cout << "min = " << min(G[start]) << endl;
//		done = true;
	}
	printAsMatrix(G);
	for (int i=0; i<MAXVERTICES; ++i) {
		cout << "dad[" << i << "] " << dad[i] << endl;
	}
}
// Dijkstra's algorithm modified to use heap structure
void Dijkstra_heap(vector<vertex> G) {

}
// Dijkstra's algorithm with edges sorted by Heapsort
void Kruskal(vector<vertex> G) {

}
int main() {
	vector<vertex> G1(MAXVERTICES);
	vector<vertex> G2(MAXVERTICES);
	makeGraph(G1,6);					// make sparse graph
	makeGraph(G2,MAXVERTICES*0.20);		// make dense graph

	// Run 1,2,3
	// 1-Max Capacity Dijkstra's without heap structure
	// 2-Max Capacity Dijkstra's with heap structure
	// 3-Max Capacity Kruskal's with edges sorted by HeapSort
	switch (1) {
					case 1:
						cout << "Case 1: Max Capacity using Dijkstra's algorithm without the heap structure\n";
						Dijkstra(G1, 0, 5);
						break;
					case 2:
						cout << "case 2: Max Capacity using Dijkstra's algorithm with the heap structure modification\n";
						Dijkstra_heap(G1);
						break;
					case 3:
						cout << "case 3: Max Capacity using Kruskal's algorithm with edges sorted by Heapsort\n";
						Kruskal(G1);
						break;
	};

//	printGraph(G1);						// debug statements for graph and matrix views
//	printAsMatrix(G1);
//	printGraph(G2);
//	printAsMatrix(G2);

	// Testing for heap structures
//	vertex v;
//	vector<vertex> heap;
//	for (int i=0; i<20; ++i) {
//		v.value = i;
//		v.value = getRand(1,MAXWEIGHT);
//		heap.push_back(v);
////		insertHeap(heap, v);		// specific stmt for testing heapify'd structure as it was created - commented out for heapsort testing.
//	}

	// Testing data for heapsort function
//	printTree(heap);
//	printline('-',heap.size()*20);
//	heap = heapsort(heap);
//	printTree(heap);

/*	// Testing data for heap insert/delete/min stuff
	printTree(heap);
	printline('-',20*heap.size());
	deleteHeap(heap,0);
	printTree(heap);
	printline('-',20*heap.size());
	deleteHeap(heap,2);
	printTree(heap);
	printline('-',20*heap.size());
	deleteHeap(heap,2);
	printTree(heap);
	printline('-',20*heap.size());
	deleteHeap(heap,2);
	printTree(heap);
	printline('-',20*heap.size());
	deleteHeap(heap,2);
	printTree(heap);
	printline('-',20*heap.size());
	deleteHeap(heap,2);
	printTree(heap);
	printline('-',20*heap.size());
	deleteHeap(heap,2);
	printTree(heap);
	printline('-',20*heap.size());
	deleteHeap(heap,2);
	printTree(heap);
	printline('-',20*heap.size());
	deleteHeap(heap,2);
	printTree(heap);
	printline('-',20*heap.size());
	deleteHeap(heap,2);
	printTree(heap);
	printline('-',20*heap.size());
	deleteHeap(heap,2);
	printTree(heap);
	printline('-',20*heap.size());
	deleteHeap(heap,2);
	*/
	return 0;
}
