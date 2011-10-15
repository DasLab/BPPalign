/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005,2006
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef WENDARRAY_H
#define WENDARRAY_H

#undef BOUNDS

#if defined (BOUNDS)
	#include <iostream>
#endif

#include "defines.h"

//This class encapsulates the w3 and w5 arrays for Dynalign
class wendarray {
private:
	short N1,N2;
	//short *itN2dN1; // array to store i times N2 divided by N1 for
                  // values of i
	short *Lowlimit;

#if defined (BOUNDS)
	//Bounds checking is on.
	short *Highlimit;
#endif

	public:
		integersize **array;
		wendarray();
		wendarray(short n1, short n2, short *lowlimit, short *highlimit);
		void allocate(short n1, short n2, short *lowlimit, short *highlimit);
		~wendarray();
		integersize &f(short i, short j);
};

inline integersize &wendarray::f(short i, short j) {

#if defined (BOUNDS)
	//bounds checking is on.  Check for valid bounds.
	if (i<0||i>N1+1) {
		std::cout << "Wend problem with i = "<<i<<"\n";

	}
	if (i<=N1) {
		if (j<Lowlimit[i]||j>Highlimit[i]) {
		 std::cout << "Wend problem with i = "<<i<<" and j = "<<j<<"\n";
		}
	}
	else if (j<0||j>N2+1) {
		 std::cout << "Wend problem with i = "<<i<<" and j = "<<j<<"\n";
	} 
#endif

  return array[i][j/*-itN2dN1[i] + M*/];
}

#endif
