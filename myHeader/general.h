
#ifndef GEN_H_
#define GEN_H_
//--------------------------------------------------
// #include <stdlib.h>	// standard C library function
//--------------------------------------------------
#include <cmath>	// math
#include <time.h>   // time
#include <fstream>	// file i/o
#include <iostream>	// i/o stream
#include <strstream>	// char *stream
#include <iomanip>	// i/o format
#include <string>  	// Standard string library
#include <map>   	// STL associative container
#include <limits>   	// Noeric limit
#include <algorithm> 	// STL algorithm
#include <vector>    	// STL vector
#include <functional>	// STL functional
#include <cmath>
#include <iterator>
#include <algorithm>
//--------------------------------------------------
// #include <unistd.h>
// #include <getopt.h>
// #include <stdio.h>
//--------------------------------------------------
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

// Section 1.2.	Constants
#define MAX(a,b) a>b? a:b
#define MIN(a,b) a>b? b:a
#define CONV (180./acos(-1.))
#define Rvdw 1.8	// maximum van der Waals radius
#define Dvdw 3.6	// maximum van der Waals overlaping distance
#define MAX_HOMOLOGY 100 // maximum homologous sequence number

// Section 1.4. Global functions
int my_lowercase(char *s);
int my_chopSpace(char *s);
void my_getopt(int argc, char *argv[]);
char compliment_nucleotide(char ch);
int* getRandNum(int upper_bound,int numbers);

template <class Type>
void myswap (Type& first,Type& second)
{
	Type tmp = first;
	first = second;
	second = tmp;
}

template <typename T>
T MAX3(T input1,T input2, T input3)
{
	T tmpmax1 = input1 >= input2 ? input1 : input2;
	T tmpmax2 = tmpmax1 >= input3 ? tmpmax1 : input3;
	return tmpmax2;
}

#endif /*GEN_H_*/

