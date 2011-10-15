//Bimol provides bimolecular secondary structure prediction where no intramolecular pairing is allowed.
//Written by Laura DiChiacchio in 2008.

#include "algorithm.h"
#include "rna_library.h"
#include "structure.h"
#include "stackstruct.h"
#include "stackclass.h"
#include "platform.h"
#include "defines.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
using namespace std;

void bimoltracebackV(int i, int j, int jct3, int k, int N1, int N2, int maxloop, short **V, structure *ct3, datatable *data);
void bimoltracebackVp(int i, int j, int jct3, int k, int N1, int N2, int maxloop, short **Vp, structure *ct3, datatable *data);


bool inc[6][6]={{false,false,false,false,false,false}, {false,false,false,false,true,false},
                {false,false,false,true,false,false}, {false,false,true,false,true,false}, {false,true,false,true,false,true}, {false,false,false,false,false,false}};

                
void bimol(structure *ct1, structure *ct2, structure *ct3, int maxloop, int maxtracebacks, int percent, int windowsize, 
           datatable *data) {
  
  int i, j, ip, jp, jct3, jpct3;
  int big, biggest, end;
  int N1 = ct1->numofbases;
  int N2 = ct2->numofbases;
  int maxseq = ct1->numofbases + ct2->numofbases + 3;
    
  int a, b, k, k1, k2, k3;
  int vmin, crit, num, numbp;
  int q, up, ir, c, cur;
  int iret, jret, sort;
  int *energy, *heapi, *heapj; 
  
  ct3->allocate(maxseq);
  //ct3->allocatestructure();
  ct3->numofbases = maxseq;

  
  
  //ct3 must contain a sequence that concatonates the sequences of ct1 and ct2  
  strcpy(ct3->ctlabel[1], ct1->ctlabel[1]);
  ct3->ctlabel[1][strlen(ct3->ctlabel[1])-1] = '_';
  strcat(ct3->ctlabel[1], ct2->ctlabel[1]); 
  
  for (i=1; i<=N1; i++) {
    ct3->numseq[i] = ct1-> numseq[i];
    ct3->nucs[i] = ct1->nucs[i];
    ct3->hnumber[i] = i;
  }
  
  //add intermolecular linker
  for (i=N1+1; i<(N1+4); i++) {
    ct3->numseq[i] = 5;
    ct3->nucs[i] = 'I';
    ct3->hnumber[i] = 0;
  }
  
  for (i=N1+4, j=1; i<(N1+N2+4); i++, j++) {
    ct3->numseq[i] = ct2->numseq[j];
    ct3->nucs[i] = ct2->nucs[j];
    ct3->hnumber[i] = j;
  }
        
  //dynamic allocation of memory for V[i][j], the mfe from 1 to i in sequence 1 and j to N2 in sequence 2, where i is paired to j.
  short **V;
  V = new short *[N1+1];
  for (int k=0; k<=N1; k++) {
    V[k] = new short [N2+1];
  }
   
  //dynamic allocation of memory for Vp[i][j], the mfe from i to N1 in sequence 1 and 1 to j in sequence 2, where i is paired to j.
  short **Vp;
  Vp = new short *[N1+1];
  for (int k = 0; k <= N1; k++) {
    Vp[k] = new short[N2+1];
  }
  
  //main loops over i and j to fill the V array
  for (i=N1; i >= 1; i--) {
    for (j=1, jct3 = N1+4; j <= N2; j++, jct3++) {

      //initialize V[i][j]
      V[i][j] = INFINITE_ENERGY;
      
	  //first, check if this is a valid pair
      if (inc[ct3->numseq[i]][ct3->numseq[jct3]]) {
        
        //check case for a terminal mismatch on a starting pair
        if (i < N1 && j > 1) {
          V[i][j] = min(V[i][j], data->tstack[ct3->numseq[i]][ct3->numseq[jct3]][ct3->numseq[i+1]][ct3->numseq[jct3-1]] +   
                        penalty(i, jct3, ct3, data));
        } 
        //i+1 dangles if nothing can stack on j
        else if (i < N1 && j == 1) {
          V[i][j] = min(V[i][j], erg4(i, jct3, i+1, 1, ct3, data, false) + penalty(i, jct3, ct3, data));
        }  
        //j-1 dangles if nothing can stack on i
        else if (i == N1 && j > 1) {
          V[i][j] = min(V[i][j], erg4(i, jct3, jct3-1, 2, ct3, data, false) + penalty(i, jct3, ct3, data));
        }
		//here nothing can stack because this pair involves terminal nucleotides
        else (V[i][j] = min(V[i][j], penalty(i, jct3, ct3, data)));
        
		//now loop over all possible closing pairs at the far end of a base pair stack, bulge loop, or interior loop
        biggest = N1-i+j-3; 
        biggest = min(biggest, maxloop);
        
        for (int a = 0; a <= biggest; a++) {
          end = i+a+1;
          big = min(end, N1);
          for (ip = max(i+1, i-j+a+3); ip <= big; ip++) {
            jp = j+ip-i-a-2;
            jpct3 = jp+N1+3;
             
			//make sure this pair at the far end is valid
            if (inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
              
              //base pair stack
			  if (a==0){//ip == i+1 && jp == j-1) {
                V[i][j] = min(V[i][j], erg1(i, jct3, ip, jpct3, ct3, data)+V[ip][jp]);
              }
              
              //either internal loop or bulge
              else V[i][j] = min(V[i][j], erg2(i, jct3, ip, jpct3, ct3, data, 0, 0)+V[ip][jp]);
              
              
            }
          }//end loop over ip
        
        }//end loop over a
      }//end check on valid pair
    }//end loop over j
  }//end loop over i
    
  
  //main loops over i and j for Vp, V calculations in reverse direction 
  for (i=1; i <= N1; i++) {
    for (j = N2, jct3 = maxseq; j >= 1; j--, jct3--) {

	  //initialize Vp
      Vp[i][j] = INFINITE_ENERGY; 
      
	  //Make sure this is a valid pair
      if (inc[ct3->numseq[i]][ct3->numseq[jct3]]) {
    
        //consider a terminal mm
        if (i > 1 && j < N2) {
          Vp[i][j] = min(Vp[i][j], data->tstack[ct3->numseq[jct3]][ct3->numseq[i]][ct3->numseq[jct3+1]][ct3->numseq[i-1]] + 
                         penalty(jct3, i, ct3, data));
        }   
        //i-1 dangles if nothing can stack on j
        else if (i > 1 && j == N2) {
          Vp[i][j] = min(Vp[i][j], erg4(jct3, i, i-1, 2, ct3, data, false) + penalty(jct3, i, ct3, data));
        }
        //j+1 dangles if nothing can stack on i
        else if (i == 1 && j < N2) {
          Vp[i][j] = min(Vp[i][j], erg4(jct3, i, jct3+1, 1, ct3, data, false) + penalty(jct3, i, ct3, data));
        }     
        //terminal bp, so no stacks are allowed
        else (Vp[i][j] = min(Vp[i][j], penalty(jct3, i, ct3, data)));
        

		//Consider an internal loop, bulge loop, or helical stack closed at the far end by another pair
        biggest = N2-j+i-3;
        biggest = min(biggest, maxloop);
     
        for (int a = 0; a <= biggest; a++) {
          end = j+a+1;
          big = min(end, N2); 
          for (jp = max(j+1, j-i+a+3), jpct3 = jp+N1+3; jp <= big; jp++, jpct3++) {
            ip = i+jp-j-a-2;
            
			//check that the pair at the far end is valid
            if (inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
            
              //bp stack
			  if (a==0){//ip == i-1 && jp == j+1) {
                Vp[i][j] = min(Vp[i][j], erg1(jct3, i, jpct3, ip, ct3, data)+Vp[ip][jp]);
              }
              
              //internal loop or bulge
              else Vp[i][j] = min(Vp[i][j], erg2(jct3, i, jpct3, ip, ct3, data, 0, 0)+Vp[ip][jp]);
            }//end valid pair at far end
        
          }//end loop over ip
        }//end loop over a
      }//end valid pair
      
    }//end loop over j
      
  }//end loop over i

   
  //begin traceback
  //for (k = 1; k <= maxtracebacks; k++) {

    //initialize vmin
    vmin = 0;
    
    
    //find lowest dG, set vmin = min dG
    for (i = N1; i >= 1; i--) {
      for (j = 1; j <= N2; j++) {
        
        if ((V[i][j] + Vp[i][j] + data->init) < vmin) {
          vmin = V[i][j] + Vp[i][j] + data->init;
        }
      }
    }
    
    //mark will keep track of bps that have been seen before
    bool **mark;
    mark = new bool *[maxseq+1];
    for (i = 0; i <= maxseq; i++) {
      mark[i] = new bool [maxseq+1]; 
    }
    
    //initialize mark
    for (int count = 1; count <= maxseq; count++) {
      for (int count2 = 1; count2 <= maxseq; count2++) {
        mark[count][count2] = false;
      }
    }
    
    //set critical value for vmin - suboptimal structures within percentage difference of lowest dG set by user
    crit = vmin - (int) ((((float) percent)/100.0)*((float) vmin));
    
    sort = 90000;
    
    energy = new int[sort+1];
    heapi = new int[sort+1];
    heapj = new int[sort+1];
    
    num = 0;
    for (i = N1; i >= 1; i--) {
      for (j = 1; j <= N2; j++) {
        if (num == sort) {
          
          //allocate more space for heapi, heapj, energy because it is about to overflow
          delete[] heapi;
          delete[] heapj;
          delete[] energy;
          sort = 10*sort;
          heapi = new int[sort+1];
          heapj = new int[sort+1];
          energy = new int[sort+1];
          i = N1;
          j = 1;
          num = 0;
        
        }
      
        //add critical values to heap
        if (V[i][j] + Vp[i][j] + data->init <= crit) {
      
        num++;
        heapi[num] = i;
        heapj[num] = j;
        energy[num] = V[i][j] + Vp[i][j] + data->init;
        
        }
      }
    }
    
    //make a heap
    for (q = 2; q <= num; q++) {
      
      cur = q;
      up = cur/2;
      
      while ((energy[cur] < energy[up] && up >= 1)) {
        
        swap(&heapi[cur], &heapi[up]);
        swap(&heapj[cur], &heapj[up]);
        swap(&energy[cur], &energy[up]);
        cur = cur/2;
        up = up/2;
        
      }
    }
    
    //sort the heap
    for (ir = num-1; ir >= 1; ir--) {
      swap(&heapi[ir+1], &heapi[1]);
      swap(&heapj[ir+1], &heapj[1]);
      swap(&energy[ir+1], &energy[1]);
      
      up = 1;
      c = 2;
      
      while (c <= ir) {
        
        if (c != ir) {
          if (energy[c+1] < energy[c]) c++;
        }
        if (energy[c] < energy[up]) {
          swap(&heapi[c], &heapi[up]);
          swap(&heapj[c], &heapj[up]);
          swap(&energy[c], &energy[up]);
          up = c;
          c = 2*c;
          
        }
        else c = ir+1;
      }
    }
    
	ct3->numofstructures = 0; //initialize the number of structures

    //for max number of structures input by user
    for (int cntr = num; cntr > 0 && ((ct3->numofstructures) < maxtracebacks); cntr--) {
      
      //find next unmarked bp
      if (!mark[heapi[cntr]][N1+3+heapj[cntr]]) {
        
        iret = heapi[cntr];
        jret = heapj[cntr];
        ct3->numofstructures++;
        ct3->checknumberofstructures();        
        jct3 = jret + N1 + 3;

		//reset all pairs in the current structure to zero
        for (i = 0; i <= maxseq; i++) {
          ct3->basepr[ct3->numofstructures][i] = 0; 
        }
        
        //perform traceback
        bimoltracebackV(iret, jret, jct3, ct3->numofstructures, N1, N2, maxloop, V, ct3, data);
        bimoltracebackVp(iret, jret, jct3, ct3->numofstructures, N1, N2, maxloop, Vp, ct3, data);
      
        ct3->energy[ct3->numofstructures] = energy[cntr];
        
        numbp = 0;

		//check if this structure is different enough from previous structures using the idea of a window
        for (k1 = 1; k1 <= maxseq; k1++) {
          
          if (k1 < (ct3->basepr[ct3->numofstructures][k1])) {
            if (!(mark[k1][ct3->basepr[ct3->numofstructures][k1]])) numbp++;
          }
        }
        for (k1 = 1; k1 <= maxseq; k1++) {
          if (k1 < (ct3->basepr[ct3->numofstructures][k1])) {
            mark[k1][ct3->basepr[ct3->numofstructures][k1]] = true;
            if (windowsize > 0) {
              for (k2 = max(1, k1 - windowsize); k2 <= min(maxseq, k1 + windowsize); k2++) {
                for (k3 = max(k2, (ct3->basepr[ct3->numofstructures][k1] - windowsize)); k3 <= min(maxseq, (ct3->basepr[ct3->numofstructures][k1] + windowsize)); k3++) {
                      
                      mark[k2][k3] = true;
                      
                    }
              }
            }
          }
        }
        
        if (numbp <= windowsize && ct3->numofstructures > 1) {
          ct3-> numofstructures--;
        }
        else {
          
          if (ct3->numofstructures != 1) strcpy(ct3->ctlabel[ct3->numofstructures], ct3->ctlabel[1]);
        } 
        
      }
    }
        
    de_allocate (mark, maxseq+1);
    delete[] energy;
    delete[] heapi;
    delete[] heapj;
    
  //}
    
    for (int l = 0; l <= N1; l++) {
      delete[] V[l];
      delete[] Vp[l];
    }
    
    delete[] V;
    delete[] Vp;
    
}



void bimoltracebackV(int i, int j, int jct3, int k, int N1, int N2, int maxloop, short **V, structure *ct3, datatable *data) {
  
  bool done, found;
  int jpct3, ip, jp;
  int biggest, big, end; 
  
  //done will switch to true when the end of this traceback is complete
  done = false; 

  while (!done) {

    //register the current base pair and find the next, if there are any
     ct3->basepr[k][i] = jct3;
     ct3->basepr[k][jct3] = i;
      
	 
	 
      
	 if (i!=N1&&j!=1) {
		 if (V[i][j] == data->tstack[ct3->numseq[i]][ct3->numseq[jct3]][ct3->numseq[i+1]][ct3->numseq[jct3-1]] + penalty(i, jct3, ct3, data)) {
          
			done = true;
		 }
          
     }
	 else if (i!=N1) {
		  
		 if (V[i][j] == erg4(i, jct3, i+1, 1, ct3, data, false) + penalty(i, jct3, ct3, data)) {
          
			done = true;
          
		}
	 }
	 else if (j!=1){
		 if (V[i][j] == erg4(i, jct3, jct3-1, 2, ct3, data, false) + penalty(i, jct3, ct3, data)) {
        
			done = true;
        
		 }
	 }
	 else {
		 if (V[i][j] == penalty(i, jct3, ct3, data)) {
        
			done = true;
		 }
        
	 }

	 if (!done) {
		//look for the next pair because we are not yet done tracing back  

	    //found will change to true when the next pair is found
        found = false;
        
        biggest = N1-i+j-3; 
        biggest = min(biggest, maxloop);
        
        for (int a = 0; a <= biggest && !found; a++) {
          end = i+a+1;
          big = min(end, N1);
          for (ip = max(i+1, i-j+a+3); ip <= big  && !found ; ip++) {
            jp = j+ip-i-a-2;
            jpct3 = jp+N1+3;
            
            if (ip == i+1 && jp == j-1 && V[i][j] == erg1(i, jct3, ip, jpct3, ct3, data)+V[ip][jp] && 
                 inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
              
              i++;
              j--;
              jct3--;
              found = true;
                         
            }
            else if (V[i][j] == erg2(i, jct3, ip, jpct3, ct3, data, 0, 0)+V[ip][jp] &&  
                     inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
              
              i = ip;
              j = jp;
              jct3 = jpct3;
              found = true;
            }
          }
          
        }
        
        if (found == false) {
          found = true;
          cerr << "Error in tracebackV at "<< i << " " << jct3 << V[i][j] << "\n";
          done = true;//throw in the towel because an error was found 
        }
      }//end if !done
	}//end while !done
}
      
void bimoltracebackVp(int i, int j, int jct3, int k, int N1, int N2, int maxloop, short **Vp, structure *ct3, datatable *data) { 
  
  bool done, found;
  int jpct3, ip, jp;
  int biggest, big, end;
  
  //done will switch to true when this traceback is complete
  done = false;
  
  while (!done) {
    
    //record the current base pair and then go look for the next
    ct3->basepr[k][i] = jct3;
    ct3->basepr[k][jct3] = i;
    
	if (i!=1&&j!=N2) {
		if(Vp[i][j] == data->tstack[ct3->numseq[jct3]][ct3->numseq[i]][ct3->numseq[jct3+1]][ct3->numseq[i-1]] + penalty(jct3, i, ct3, data)) {
      
			done = true;
		}
      
    }
	else if (i!=1) {
		if (Vp[i][j] == erg4(jct3, i, i-1, 2, ct3, data, false) + penalty(jct3, i, ct3, data)) {
        
			done = true;
		}
        
    }
	else if (j!=N2) {
		if (Vp[i][j] == erg4(jct3, i, jct3+1, 1, ct3, data, false) + penalty(jct3, i, ct3, data)) {
      
			done = true;
		}
      
    }
	else {
		if (Vp[i][j] == penalty(jct3, i, ct3, data)) {
    
			done = true;
		}
    }
	if (!done) {
      
      found = false;
      

      biggest = N2-j+i-3;
      biggest = min(biggest, maxloop);
     
      for (int a = 0; a <= biggest && !found; a++) {
        end = j+a+1;
        big = min(end, N2); 
        for (jp = max(j+1, j-i+a+3), jpct3 = jp+N1+3; jp <= big && !found; jp++, jpct3++) {
          ip = i+jp-j-a-2;
          
          if (ip == i-1 && jp == j+1 && Vp[i][j] == erg1(jct3, i, jpct3, ip, ct3, data)+Vp[ip][jp] && 
              inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
            
            i--;
            j++;
            jct3++;
            found = true;
            
          }
            
          else if (Vp[i][j] == erg2(jct3, i, jpct3, ip, ct3, data, 0, 0)+Vp[ip][jp] && 
                   inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
            
            i = ip;
            j = jp;
            jct3 = jpct3;
            found = true;
            
          }
        }
      }
      
      if (found == false) {
        found = true;
        cerr << "Error in tracebackVp at " << i << " " << jct3 << " " << Vp[i][j] << "\n";
        done = true; 
      } 
    
    }//end if !done
      
  }//end while !done
}


