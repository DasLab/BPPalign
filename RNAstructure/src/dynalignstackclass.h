/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef DYNALIGNSTACKCLASS_H
#define DYNALIGNSTACKCLASS_H

#include "defines.h"

class dynalignstackclass {
	short **stack;
	int size, max;
	integersize *stackenergy;
	bool *openness;
	void allocate_stack();

public:
	dynalignstackclass(short int stacksize = 50);

	bool pull(short int *i,short int *j, short int *a, short int *b, 
            integersize *energy, bool *open);
	void push(short int i,short int j, short int a, 
            short int b, integersize energy, bool open=false);
	
	void delete_array();
	~dynalignstackclass();
};

#endif
