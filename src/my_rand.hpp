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

random_device rd;	// seed for random generator mt19937 -> Mersenne Twister 19937
/*
 * Generates a custom random number between integer values: min and max
 * Based on the C++11 random number implementation
 */
int getRand(int min, int max) {
	mt19937 mt(rd());
	uniform_int_distribution<int> dist(min,max);
	return dist(mt);
}

#endif /* MY_RAND_HPP_ */
