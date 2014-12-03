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
#include <algorithm>
#include <time.h>
#include <fstream>
#include "my_rand.hpp"

using namespace std;

const int MAXVERTICES = 5000;
const int SPARSECONNECT = 6;	// number of sparse connections
const int MAXWEIGHT = 1000;		// max weight of graph edges
long run;						// counter for steps on run
time_t starttime;				// timer for second evaluation on run
time_t endtime;

enum {unseen, intree, fringe};

struct edge {
	int from;
	int to;
	int weight;
};
struct vertex {
	vector <edge> edges;
	int index;
	int value;
	int status;
};
void printline(char c, int num) {
	for (int i=0; i<num; ++i) cout << c;
	cout << endl;
}
void printGraph(vector <vertex> G) {
	u_int i, j;
	for (i=0; i<G.size(); i++) {
		cout << "v[" << G[i].index << "] "<<G[i].value<<" = ";
		for (j=0; j<G[i].edges.size(); j++) {
			cout<<'('<< G[i].edges[j].from << '-' << G[i].edges[j].to<<'|'<<G[i].edges[j].weight << ")  ";
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
				graph[row][col.to] = col.weight;
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
	run++;										// this run only counted after graphs are made (not counting make)
	if (vertex == newedge.to) return true;		// prevents cycle to same vertex
	for (u_int i=0; i<e.size(); ++i) {
		run++;
		if (e.at(i).to == newedge.to) return true;
	}
	return false;
}
// Find matching edge and return index or -1 if not in edges list
int my_find(vector<edge> e, int match) {
	for (u_int i=0; i<e.size(); ++i) {
		run++;
		if (match == e.at(i).to) return i;
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
			newedge1.to = getRand(0,MAXVERTICES-1);
			newedge1.weight = getRand(1,MAXWEIGHT);
			newedge1.from = i;
			newedge2.to = i;
			newedge2.from = newedge1.to;
			newedge2.weight = newedge1.weight;
			while (duplicate(G.at(i).edges, i, newedge1) ||
					G.at(newedge1.to).edges.size() >= connections) {		// keep random weight, get new connect
				newedge1.to = getRand(0,MAXVERTICES-1);
				newedge2.to = i;
				newedge2.from = newedge1.to;
				++badLoop;
				if (badLoop >= 2*MAXVERTICES) {
					temp = G.at(newedge1.to).edges.front().to;
					G.at(newedge1.to).edges.erase(G[newedge1.to].edges.begin());
					G[temp].edges.erase(G[temp].edges.begin()+my_find(G.at(temp).edges, newedge1.to));
					--j;	// decrement outside loop per set delete to ensure enough add's to graph
				}
			}
			badLoop=0;
			G.at(i).edges.push_back(newedge1);
			G.at(newedge1.to).edges.push_back(newedge2);
		}
	}
	for (int i=0; i<MAXVERTICES; ++i) {
		G[i].index = i;
	}
}
// Function adds a path with random weights between start and end thus
// ensuring there is at least 1 connection between any 2 vertices.
// May add 2 paths (un-directed) to every vertex but 1st & last,
// won't add duplicate paths or cycles.
void addPath(vector<vertex> &G) {
	edge e1, e2;
	for (u_int i=0; i<MAXVERTICES-1; ++i) {
		e1.from = i;
		e1.to = i+1;
		e2.from = i+1;
		e2.to = i;
		e1.weight = getRand(1,MAXWEIGHT);
		e2.weight = e1.weight;
		run+=6;
		if (!duplicate(G[i].edges,i,e1)) {
//			cout << "adding to v["<<i<<"] - edge from["<<e1.from<<"] to ["<<e1.to<<"]\nadding to v["<<i+1<<"] - edge from ["<<e2.from<<"] to ["<<e2.to<<"]\n";
			G[i].edges.push_back(e1);
			G[i+1].edges.push_back(e2);
			run+=2;
		}
	}
}
// Swapping function for edges
void swapem(edge &a, edge &b) {
	edge temp;
	temp = a;
	a = b;
	b = temp;
	run+=3;
}
// Swapping function vertex
void swapem(vertex &a, vertex &b) {
	vertex temp;
	temp = a;
	a = b;
	b = temp;
	run+=3;
}
// Top heap return function, returns top level element of the heap
// This is min or max value depending on type of heap used.
vertex rootHeap(vector<vertex> heap) {
	run++;
	if (heap.size() > 0) return heap[0];
	else {
		cout << "Error - heap is empty when returning rootHeap()\n";
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
// Overloaded heapify for edges for Kruskal's algorithm
void heapify(vector<edge> &heap) {
	int current, parent;
	run++;
	current = heap.size()-1;	// set current to last node (end of array)
	while (current >= 0) {
		parent = (current-1)/2;		// set parent to root of current
		if (heap[current].weight > heap[parent].weight) swapem(heap[current],heap[parent]);
		--current;					// cycle down in nodes until reach root
		run+=3;
	}
}
// heapify for vector of vertex - for Dijkstra_heap function
void heapify(vector<vertex> &heap) {
	int current, parent;
	run++;
	current = heap.size()-1;	// set current to last node (end of array)
	while (current >= 0) {
		parent = (current-1)/2;		// set parent to root of current
		if (heap[current].value > heap[parent].value) swapem(heap[current],heap[parent]);//< heap[parent].value) swapem(heap[current],heap[parent]);
		--current;					// cycle down in nodes until reach root
		run+=3;
	}
}
// Now sort down - called from deleteHeap with index
void moveDown(vector<vertex> &heap, int i) {
	int parent, end, child;
	parent = i;
	end = heap.size()-1;
	run+=2;
	while (parent * 2 + 1 <= end) {
		run+=4;
		child = parent * 2 + 1;
		if (child+1 <= end && heap[child].value < heap[child+1].value) ++child;//> heap[child+1].value) ++child;
		if (heap[parent].value < heap[child].value) {//> heap[child].value) {
			swapem(heap[parent],heap[child]);
			parent = child;
		}
		else return;
	}
}
// Insert element in correct place in heap so don't need sort
// vector takes care of size so don't need to maintain count
void insertHeap(vector<vertex> &heap, vertex v) {
	heap.push_back(v);
	heapify(heap);
	run++;
}
void deleteHeap(vector<vertex> &heap, int idx) {
	if (heap.size() <= 0) {
		cout << "Error, heap empty on deleteHeap() call!!\n";
		return;
	}
	swapem(heap[idx],heap[heap.size()-1]);
	heap.pop_back();
	moveDown(heap, idx);
	run++;
}
// Overloaded to print heap list of edges
void printHeap(vector<edge> heap) {
	int i=0;
	for (auto e : heap) {
		cout <<"n["<<i++<< "]: from ["<<e.from<<"] to ["<<e.to<<"]: "<<e.weight<<endl;
	}
}
// Overloaded to print heap list of vertices
void printHeap(vector<vertex> heap) {
	for (u_int i=0; i<heap.size(); ++i) {
		cout << "H[" << i << "]=" << heap[i].value << "      ";
	}
	cout << endl;
}
// Heapsort function - designed to function for vector of type edge (heapsort edges)
vector<edge> heapsort(vector<edge> &heap) {
	vector<edge> sortedHeap;
	heapify(heap);
	int end = heap.size()-1;
	run++;
	while (end > 0) {
		sortedHeap.push_back(heap[0]);
		swapem(heap[0],heap[end]);		// not using deleteHeap here because we have to re-heapify anyway
		heap.pop_back();				// instead just push to new storage, swap to back and pop.
		heapify(heap);
		--end;
		run+=3;
	}
	return sortedHeap;
}
// returns MAX value of all vertices in vector v (the index of the max)
int max(vector<vertex> v) {
	int current = -1;
	int answer = -1;
	for (u_int i=0; i<v.size(); ++i) {
		if (v[i].value > current) {
			current = v[i].value;
			answer = i;
			run+=2;
		}
	}
	if (answer == -1) cout << "Error - returning -1 for edge in max() fct\n";
//	cout << "max() returning " << answer << " = " << current << endl;
	return answer;
}
void printDad(int dad[]) {
	for (int i=0; i<MAXVERTICES; ++i) {
		cout << "dad[" << i << "] " << dad[i] << endl;
	}
}
void printPath(int dad[], int start, int end, ofstream &outf) {
	cout << "Max Path from (" << start << ") to (" << end << ") is: ";
	outf << "Max Path from (" << start << ") to (" << end << ") is: ";
	vector<int> path;
	int current = end;
	while (current != start) {
		path.push_back(dad[current]);
		if (current == -1) {
			cout << "Error, reached -1 index in printPath()\n";
			outf << "Error, reached -1 index in printPath()\n";
			break;
		}
		current = dad[current];
	}
	while (!path.empty()) {
		cout << "["<< path.back() << "] -> ";
		outf << "["<< path.back() << "] -> ";
		path.pop_back();
	}
	cout << '[' << end << "].\n";
	outf << '[' << end << "].\n";
}
// Dijkstra's algorithm without heap structure (regular)
void Dijkstra(vector<vertex> G, int start, int end, ofstream &outf) {
//	printGraph(G);
	int v;
	vector<vertex> list;
	int dad[MAXVERTICES];
	for (int i=0; i<MAXVERTICES; ++i) {
		run++;
		dad[i] = -1;
	}
	for (auto &v : G) {
		v.status = unseen;
		v.value = 0;
		run+=2;
	}
	G[start].status = intree;
	G[start].value = MAXWEIGHT;
	dad[start] = start;
	run+=3;
	for (auto &e : G[start].edges) {
		G[e.to].status = fringe;
		dad[e.to] = start;
		G[e.to].value = e.weight;
		list.push_back(G[e.to]);
		run+=4;
	}
	while (list.size() > 0) {
		v = max(list);
//		cout << "v: " << v << " index: " << list[v].index<< "\nsetting G["<< list[v].index <<"].status to intree\n";
		G[list[v].index].status = intree;
		run++;
//		printGraph(list);
//		cout << "looping through G["<<list[v].index<<"] (v["<<list[v].index<<"]) edges:\n";
		for (auto &e : G[list[v].index].edges) {
//			cout << "if g["<<e.to<<"] unseen (state=" << G[e.to].status<<") then:\n";
//			if (G[e.to].status != unseen) { //output stmts for else below
//				cout << "G["<<e.to<<"].status =? fringe, G["<<e.to<<"].value="<<G[e.to].value<<" <? G["<<v<<"].value("<<G[v].value<<")+e.weight("<<e.weight<<")\n";
//			}
			run++;
			if (G[e.to].status == unseen) {
//				cout << "setting status G["<<e.to<<"] to fringe\n";
				G[e.to].status = fringe;
//				cout << "setting dad["<<e.to<<"] to ["<<e.from<<"]\n";
				dad[e.to] = e.from;
//				cout << "setting G["<<e.to<<"].value to G["<<v<<"].value ("<<G[v].value<<") + e.weight ("<< e.weight << ")\n";//list["<<v<<"].value("<<list[v].value<<") + e.weight:("<<e.weight<<")\n";
				G[e.to].value = e.weight;//G[v].value + e.weight;
				list.push_back(G[e.to]);
				run+=4;
			}
			else if (G[e.to].status == fringe && G[e.to].value < G[v].value + e.weight) {
//				cout << "fringe\n" << "erasing...begin + "<<v<<"...\n";
				list.erase(list.begin()+v);
				run+=list.size()-v+3;
				G[e.to].value = e.weight;//G[v].value + e.weight;
				dad[e.to] = e.from;
				list.push_back(G[e.to]);
			}
		}
//		printGraph(list);
//		cout << "erasing...begin + "<<v<<"...\n";
		list.erase(list.begin()+v);
//		printGraph(list);
		run+=list.size()-v+3;
	}
	cout << "Printing path for Dijkstra:\n";
	outf << "Printing path for Dijkstra:\n";
	printPath(dad,start,end,outf);
}
// Dijkstra's algorithm modified to use heap structure
void Dijkstra_heap(vector<vertex> G, int start, int end, ofstream &outf) {
//	printGraph(G);
	int v;
	vertex heapnode;
	vector<vertex> heap;
	int dad[MAXVERTICES];
	for (int i=0; i<MAXVERTICES; ++i) {
		run++;
		dad[i] = -1;
	}
	for (auto &v : G) {
		v.status = unseen;
		v.value = 0;
		run+=2;
	}
	G[start].status = intree;
	G[start].value = MAXWEIGHT;
	dad[start] = start;
	run+=3;
	for (auto &e : G[start].edges) {
		G[e.to].status = fringe;
		dad[e.to] = start;
		G[e.to].value = e.weight;
		insertHeap(heap,G[e.to]);
		run+=4;
	}
	while (!heap.empty()) {
		heapnode = rootHeap(heap);
		v = heapnode.value;
//		cout << "v: " << v << " index: " << heapnode.index<<"\nsetting G["<< heapnode.index <<"].status to intree\n";
		G[heapnode.index].status = intree;
//		printTree(heap);
//		cout << "looping through G["<<heapnode.index<<"] (v["<<heapnode.index<<"]) edges:\n";
		run+=2;
		for (auto &e : G[heapnode.index].edges) {
//			cout << "if g["<<e.to<<"] unseen (state=" << G[e.to].status<<") then:\n";
//			if (G[e.to].status != unseen) { //output stmts for else below
//				cout << "G["<<e.to<<"].status =? fringe, G["<<e.to<<"].value="<<G[e.to].value<<" <? G["<<v<<"].value("<<G[v].value<<")+e.weight("<<e.weight<<")\n";
//			}
			run++;
			if (G[e.to].status == unseen) {
//				cout << "setting status G["<<e.to<<"] to fringe\n";
				G[e.to].status = fringe;
//				cout << "setting dad["<<e.to<<"] to ["<<e.from<<"]\n";
				dad[e.to] = e.from;
//				cout << "setting G["<<e.to<<"].value to G["<<v<<"].value ("<<G[v].value<<") + e.weight ("<< e.weight << ")\n";//list["<<v<<"].value("<<list[v].value<<") + e.weight:("<<e.weight<<")\n";
				G[e.to].value = e.weight;//G[v].value + e.weight;
				insertHeap(heap,G[e.to]);
				run+=4;
			}
			else if (G[e.to].status == fringe && G[e.to].value < G[v].value + e.weight) {
//				cout << "fringe\n" << "erasing...begin + "<<v<<"...\n";
				deleteHeap(heap,0);
				G[e.to].value = e.weight;//G[v].value + e.weight;
				dad[e.to] = e.from;
				insertHeap(heap,G[e.to]);
				run+=4;
			}
		}
//		printTree(heap);
//		cout << "erasing...begin + "<<v<<"...\n";
		deleteHeap(heap,0);
		run++;
//		printTree(heap);
	}
	cout << "Printing path for Dijkstra_heap:\n";
	outf << "Printing path for Dijkstra_heap:\n";
	printPath(dad,start,end,outf);
}
// Find root function for Kruskal's algorithm
int findRoot(int dad[], int v) {
	while (dad[v] != -1) {
		run++;
		v = dad[v];
	}
	return v;
}
// Modified union function for Kruskals as given in class
void myUnion(int dad[], int rank[], int r1, int r2) {
	if (rank[r1] > rank[r2]) dad[r2] = r1;
	if (rank[r1] < rank[r2]) dad[r1] = r2;
	if (rank[r1] == rank[r2]) {
		dad[r2] = r1;
		rank[r1]++;
		run+=2;
	}
	run+=3;
}
// Dijkstra's algorithm with edges sorted by Heapsort
void Kruskal(vector<vertex> G, int start, int end, ofstream &outf) {
	vector<edge> edges;
	vector<edge> path;
	edge e;
	int dad[MAXVERTICES];
	int rank[MAXVERTICES];
	int fill_count = 0;		// counter to stop adding edges
	u_int i;
	for (i=0; i<MAXVERTICES; i++) {
		dad[i] = -1;		// root nodes listed as -1 instead of 0 - uses 0 as a node
		rank[i] = 0;
		run+=2;
	}
	for (auto v : G)
		for (auto e : v.edges) {
			edges.push_back(e);
			run++;
		}
	edges = heapsort(edges);
//	cout << "setting dad["<<edges[0].to<<"] = ["<<edges[0].from<<"]\n";
//	dad[edges[0].to] = edges[0].from;
	for (i=0; i<edges.size(); ++i) {
		e = edges[i];
		run++;
//		cout << "if findroot("<<e.from<<")="<<findRoot(dad,e.from)<<" != findroot("<<e.to<<")="<<findRoot(dad,e.to)<<" do union...\n";
//		cout << "setting "<<e.from<< " to "<<e.to<<endl;
		if (findRoot(dad,e.from) != findRoot(dad,e.to)) {
			if (findRoot(dad,e.from) == -1 || findRoot(dad,e.to) == -1) fill_count++;
//			cout << "calling union\n";
			path.push_back(e);
			myUnion(dad,rank,e.from,e.to);
			run+=2;
		}
		else if (dad[e.to] != -1 && dad[e.from] == -1){
//			cout << "else dad["<<e.to<<"] = "<<e.from<<"\n";
			dad[e.to] = e.from;
			run++;
		}

		if (fill_count > MAXVERTICES) break;	// stop processing - have added enough
//		cout << "dad["<<e.to<<"] = "<<dad[e.to]<<endl;
	}
//	printDad(dad);
//	for (auto p : path) {
//		cout << p.from << ":"<<p.to<<"|"<<p.weight<<endl;
//	}
	printPath(dad, start, end, outf);
}
int main() {
	vector<vertex> G1(MAXVERTICES);
	vector<vertex> G2(MAXVERTICES);
	ofstream outf;
	ofstream outf2;

	// Run 1 or 2 for switch statement below
	// 1-Run algorithms for Sparse Graph Structure
	// 2-Run algorithms for Dense Graph Structure
	//   Each runs 5 times for same graph structure, output to files and console

	int algorithm = 2;	// 1 for sparse, 2 for dense

	outf.open("sparse_output.txt");
	outf2.open("dense_output.txt");

	for (int i=1; i<=5; i++) {
		int start = getRand(0,MAXVERTICES-1);
		int end = getRand(0,MAXVERTICES-1);
		switch (algorithm) {
		case 1:
			outf << "Sparse Graph output:\nGraph vertices = "<<MAXVERTICES<<", finding path from node "<<start<<" to "<<end<<" - run #"<<i<<".\n";
			cout << "Sparse Graph output, also generated in 'sparse_output.txt' - run #"<<i<<"\n";
			cout << "Graph vertices = "<<MAXVERTICES<<", finding path from node "<<start<<" to "<<end<<".\n";
			makeGraph(G1,6);					// make sparse graph
			addPath(G1);
			run = 0;							// run count
			starttime = clock();
			Dijkstra(G1,start,end,outf);
			endtime = clock();
			cout << "Run count (steps): "<<run<<".\n";
			outf << "Run count (steps): "<<run<<".\n";
			cout << "Time of run: "<< (double)(endtime-starttime) << " seconds.\n";
			outf << "Time of run: "<< (double)(endtime-starttime) << " seconds.\n";
			run = 0;
			starttime = clock();
			Dijkstra_heap(G1,start,end,outf);
			endtime = clock();
			cout << "Run count (steps): "<<run<<".\n";
			outf << "Run count (steps): "<<run<<".\n";
			cout << "Time of run: "<< (double)(endtime-starttime) << " seconds.\n";
			outf << "Time of run: "<< (double)(endtime-starttime) << " seconds.\n";
			run = 0;
			starttime = clock();
			Kruskal(G1,start,end,outf);
			endtime = clock();
			cout << "Run count (steps): "<<run<<".\n";
			outf << "Run count (steps): "<<run<<".\n";
			cout << "Time of run: "<< (double)(endtime-starttime) << " seconds.\n";
			outf << "Time of run: "<< (double)(endtime-starttime) << " seconds.\n";
			break;
		case 2:
			outf2 << "Dense Graph output:\nGraph vertices = "<<MAXVERTICES<<", finding path from node "<<start<<" to "<<end<<" - run #"<<i<<".\n";
			cout << "Dense Graph output, also generated in 'dense_output.txt' - run #"<<i<<"\n";
			cout << "Graph vertices = "<<MAXVERTICES<<", finding path from node "<<start<<" to "<<end<<".\n";
			makeGraph(G2,MAXVERTICES*0.20);		// make dense graph
			addPath(G2);
			run = 0;							// run count
			starttime = clock();
			Dijkstra(G2,start,end,outf2);
			endtime = clock();
			cout << "Run count (steps): "<<run<<".\n";
			outf2 << "Run count (steps): "<<run<<".\n";
			cout << "Time of run: "<< (double)(endtime-starttime) << " seconds.\n";
			outf2 << "Time of run: "<< (double)(endtime-starttime) << " seconds.\n";
			run = 0;
			starttime = clock();
			Dijkstra_heap(G2,start,end,outf2);
			endtime = clock();
			cout << "Run count (steps): "<<run<<".\n";
			outf2 << "Run count (steps): "<<run<<".\n";
			cout << "Time of run: "<< (double)(endtime-starttime) << " seconds.\n";
			outf2 << "Time of run: "<< (double)(endtime-starttime) << " seconds.\n";
			run = 0;
			starttime = clock();
			Kruskal(G2,start,end,outf2);
			endtime = clock();
			cout << "Run count (steps): "<<run<<".\n";
			outf2 << "Run count (steps): "<<run<<".\n";
			cout << "Time of run: "<< (double)(endtime-starttime) << " seconds.\n";
			outf2 << "Time of run: "<< (double)(endtime-starttime) << " seconds.\n";
			break;
		};
	}
	cout << "\nDone.\n";
	outf.close();
	outf2.close();
	return 0;
}
