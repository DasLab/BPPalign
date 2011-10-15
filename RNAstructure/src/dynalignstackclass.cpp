/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#include "dynalignstackclass.h"

#include "defines.h"

void dynalignstackclass::allocate_stack() {
  stackenergy = new integersize[max];
  stack=new short int *[max];
  for (short i = 0; i < max; i++) {
    stack[i] = new short int[4];
  }
  openness = new bool [max];
}

dynalignstackclass::dynalignstackclass(short int stacksize) {
  max = stacksize;
  size = 0;
  allocate_stack();
}

bool dynalignstackclass::pull(short int *i,short int *j, short int *a, short int *b, 
                              integersize *energy, bool *open) {
  if (size==0) {
    return false;
  }

  size--;
  *i = stack[size][0];
  *j = stack[size][1];
  *a = stack[size][2];
  *energy = stackenergy[size];
  *b = stack[size][3];
  *open = openness[size];

  return true;
}
	
void dynalignstackclass::push(short int i,short int j, short int k, short int l, 
                              integersize energy, bool open){
  short index;

  if (size == max) {
    //allocate more space:
    dynalignstackclass *temp;
    temp = new dynalignstackclass(max);
    for (index=0;index<max;index++) {
      temp->push(stack[index][0],stack[index][1],stack[index][2],stack[index][3],stackenergy[index],openness[index]);
    }
    delete_array();
    max = 2*max;

    allocate_stack();
    for (index=0;index<(max/2);index++) {
      temp->pull(&stack[index][0],&stack[index][1],&stack[index][2],&stack[index][3],&stackenergy[index],&openness[index]);
    }
  }
  stack[size][0] = i;
  stack[size][1] = j;
  stack[size][2] = k;
  stackenergy[size] = energy;
  stack[size][3] = l;
  openness[size] = open;
  size++;
}
	
void dynalignstackclass::delete_array() {
  for (short i = 0; i < max; i++) {
    delete[] stack[i];
  }
  delete[] stack;

  delete[] stackenergy;
  delete[] openness;
}

dynalignstackclass::~dynalignstackclass() {
  delete_array();
}
