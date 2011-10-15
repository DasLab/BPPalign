/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005,2006
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

//The wendarray class provides functionality for encapsulating the w3 and w5 arrays.

#include "wendarray.h"

#include "defines.h"


wendarray::wendarray() {
	

}

wendarray::wendarray(short n1, short n2, short *lowlimit, short *highlimit) {
	allocate(n1, n2, lowlimit, highlimit);

}

// Functions for wendarray.  Used for w3 and w5
void wendarray::allocate(short n1, short n2, short *lowlimit, short *highlimit) {
  N1 = n1;
  N2 = n2;
  Lowlimit = lowlimit;
  
#if defined (BOUNDS)
	//Bounds checking is on.
	Highlimit=highlimit;
#endif

  short Low,High;

  array = new integersize *[n1+2];
  for (int i = 0; i <= n1 + 1; ++i) {
    
    
	if (i<=N1) {
		Low=lowlimit[i];
		High=highlimit[i];
	}
	
	else {
		Low=0;
		High=N2+1;
	}
	array[i] = new integersize[High-Low+1/*2*m+2*/];
	for (int j = 0; j <= High-Low; ++j) {
      array[i][j]=DYNALIGN_INFINITY;
    }
	array[i]=array[i]-Low;//shift pointer for access
  
  }

  //itN2dN1 = new short[N1+2];
  //for(int i = 0; i <= N1 + 1; i++) {
  //  itN2dN1[i]=i*N2/N1;
  //}
}

wendarray::~wendarray() {
  int i;
	
  for (i=0;i<=N1+1;++i) {
	  //restore pointer:
	  if(i<=N1) array[i]=array[i]+Lowlimit[i];
	  //else array[i]=array[i];
	  delete[] array[i];
  }
  delete[] array;
  

}
