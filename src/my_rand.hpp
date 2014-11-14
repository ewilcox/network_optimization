/*
 * my_rand.hpp
 *
 *  Created on: Nov 13, 2014
 *      Author: Eric Wilcox - Custom random number generator
 */

#include <random>
using namespace std;

#ifndef MY_RAND_HPP_
#define MY_RAND_HPP_

/*
 * Generates a custom random number between integer values: min and max
 * Based on the C++11 random number implementation
 */
int getRand(int min, int max) {
	int i,j;
	int NUM = 20;
	random_device rd;
//	default_random_engine mt(rd());		// this works, is lightweight according to documentation
	mt19937 mt(rd());
	uniform_int_distribution<int> dist(0,NUM);
//	cout << "***** " << rd() << endl;   // this shows the random number the random_device generates
	int a[NUM];
	for (i=0; i<NUM; ++i) {
		a[i] = 0;
	}
	for (i=0; i<1000; ++i) {
		++a[dist(mt)];
	}
	for (i=0; i<NUM; ++i) {
		cout << "a[" << i << "] ";
		for (j=0; j<a[i]; ++j) {
			cout << '*';
		}
		cout << endl;
	}
	return 0;
}

#endif /* MY_RAND_HPP_ */
