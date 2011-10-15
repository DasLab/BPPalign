/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005,2006
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 * Contributors: Arif Harmanci and Gaurav Sharma, 2006, 2007
 */

#include "dynalign.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "algorithm.h"
#include "arrayclass.h"
#include "defines.h"
#include "dynalignarray.h"
#include "dynalignheap.h"
#include "dynalignstackclass.h"
#include "forceclass.h"
#ifdef _WINDOWS
	#include "../RNAstructure_windows_interface/platform.h"
#endif //_WINDOWS
#include "rna_library.h"
#include "structure.h"

#ifdef _WINDOWS
#include "../RNAstructure_windows_interface/TProgressDialog.h"
#else

#ifdef _JAVA_GUI
#include "../RNAstructure_java_interface/SWIG/TProgressDialog.h"
#else
#include "TProgressDialog.h"
#endif // JAVA GUI

#endif //WINDOWS

#include "varray.h"
#include "wendarray.h"
#include "align_services.h"

#ifdef DYNALIGN_SMP
#include "observingtextprogressbar.h"
#include "rankconsumer.h"
#include "rankmanager.h"
#include "rankproducer.h"
#endif //DYNALIGN_SMP

using namespace std;

//#define timer //flag to turn on a timer that writes a file
#undef timer

//maximum size of unpaired nucs on each side of an internal loop
#define maxloop 20

//specify the maximum size of an internal loop
#define maxinternal 30

short int edangle5noforce(int i,int j,int ip,structure* ct,datatable* data) {
  return data->dangle[ct->numseq[j]][ct->numseq[i]][ct->numseq[ip]][2];
}

short int edangle3noforce(int i,int j,int ip,structure* ct,datatable* data) {
  return data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][1];
}

short int edangle5force(int i,int j,int ip,structure* ct,datatable* data) {
  if (ct->fcedbl[ip]) {
    return DYNALIGN_INFINITY;
  } else {
    return data->dangle[ct->numseq[j]][ct->numseq[i]][ct->numseq[ip]][2];
  }
}

short int edangle3force(int i,int j,int ip,structure* ct,datatable* data) {
  if (ct->fcedbl[ip]) {
    return DYNALIGN_INFINITY;
  } else {
    return data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][1];
  }
}

// Return the most 5' nucleotide in sequence 2 that can be aligned to
// nucleotide i in sequence 1, using the M constraint:
short lowlimit(short i, short M, short N1, short N2) {
  if (i<=N1) {
    return i*N2/N1-M;
  } else {
    return (i-N1)*N2/N1+N2-M;
  }
}
//overload to use a bool array of allowed alignments:
short lowlimit(short i, bool **allowed_alignments, short N1, short N2) {
	int j;
	
	if (i==0) return 0;

	if (i<=N1) {

		for (j=1;j<=N2;j++) {
			if (allowed_alignments[i][j]) {
				return j;
			}

		}
		return min(i,N2);//return something for program execution
	}
	else {
		for (j=1;j<=N2;j++) {
			if (allowed_alignments[i-N1][j]) return j+N2;
		}
		return min(i+N1,2*N2);
	}
	


}

// Return the most 3' nucleotide in sequence 2 that can be aligned to
// nucleotide i in sequence 1, using the M constraint:
short highlimit(short i, short M, short N1, short N2) {
  if (i<=N1) {
    return i*N2/N1+M;
  } else {
    return (i-N1)*N2/N1+N2+M;
  }
}
//overload to use a bool array of allowed alignments:
short highlimit(short i, bool **allowed_alignments, short N1, short N2) {
	int j;
	
    if (i==0) return N2;

	if (i<=N1) {
		for (j=N2;j>0;j--) {
			if (allowed_alignments[i][j]) {
				return j;
			}
		}
		return min(i,N2);//return something for program execution
	}
	else {
		for (j=N2;j>0;j--) {
			if (allowed_alignments[i-N1][j]) return j+N2;
		}
		return min(i+N1,2*N2);//return something for program execution

	}
}

short reference(short i, short k, short N, short N2, short M) {
  return k-i*N2/N+M;
}

void dynalignfceunpaired(structure *ct,char **fce,int nopair) {

  int i;

  for (i=nopair+1;i<nopair+(ct->numofbases);++i) {
    fce[jref(nopair,i,ct->numofbases)][iref(nopair,i,ct->numofbases)]=fce[jref(nopair,i,ct->numofbases)][iref(nopair,i,ct->numofbases)]|SINGLE;
  }
  for (i=1;i<nopair;++i) {
    fce[nopair][i]=fce[nopair][i]|SINGLE;
  }
  for (i=nopair+1;i<=ct->numofbases;++i) {
    fce[jref(i,nopair+ct->numofbases,ct->numofbases)][iref(i,nopair+ct->numofbases,ct->numofbases)]=
      fce[jref(i,nopair+ct->numofbases,ct->numofbases)][iref(i,nopair+ct->numofbases,ct->numofbases)]|SINGLE;
  }

}

void dynforcedbl(int dbl,structure* ct,char **fce,bool *lineardbl) {
  int i,j;

  lineardbl[dbl]=true;
  lineardbl[dbl+ct->numofbases]=true;

  for(i=dbl+1;i<=ct->numofbases;++i) {
    for (j=1;j<dbl;++j) {
      fce[i][j] = fce[i][j]|DUBLE;
    }
  }
  for(j=(dbl+(ct->numofbases)-1);j>ct->numofbases;j--) {
    for (i=dbl+1;i<=ct->numofbases;++i) {
      fce[jref(i,j,ct->numofbases)][iref(i,j,ct->numofbases)] = fce[jref(i,j,ct->numofbases)][iref(i,j,ct->numofbases)]|DUBLE;
    }
  }
}

void dynforcepair(int x,int y,structure *ct,char **v) {
  int i,j;
  //v->f(x,y) = v->f(x,y)|PAIR;
  //v->f(y,x+ct->numofbases)=v->f(y,x+ct->numofbases)|PAIR;
  for (i=y+1;i<=x-1+ct->numofbases;++i) {
    v[jref(x,i,ct->numofbases)][iref(x,i,ct->numofbases)] = v[jref(x,i,ct->numofbases)][iref(x,i,ct->numofbases)]|NOPAIR;
  }
  for (i=x;i<=y-1;++i) {
    v[jref(x,i,ct->numofbases)][iref(x,i,ct->numofbases)] = v[jref(x,i,ct->numofbases)][iref(x,i,ct->numofbases)]|NOPAIR;
  }
  for (i=1;i<=x-1;++i) {
    v[jref(i,y,ct->numofbases)][iref(i,y,ct->numofbases)] = v[jref(i,y,ct->numofbases)][iref(i,y,ct->numofbases)]|NOPAIR;
  }
  for (i=x+1;i<=y;++i) {
    v[jref(i,y,ct->numofbases)][iref(i,y,ct->numofbases)] = v[jref(i,y,ct->numofbases)][iref(i,y,ct->numofbases)]|NOPAIR;
  }
  for (i=1;i<=x-1;++i) {
    v[jref(i,x,ct->numofbases)][iref(i,x,ct->numofbases)] = v[jref(i,x,ct->numofbases)][iref(i,x,ct->numofbases)]|NOPAIR;
  }
  for (i=y+1;i<=ct->numofbases;++i) {
    v[jref(i,y+ct->numofbases,ct->numofbases)][iref(i,y+ct->numofbases,ct->numofbases)]=v[jref(i,y+ct->numofbases,ct->numofbases)][iref(i,y+ct->numofbases,ct->numofbases)]|NOPAIR;
  }
  for (i=y;i<=x-1+(ct->numofbases);++i) {
    v[jref(y,i,ct->numofbases)][iref(y,i,ct->numofbases)] = v[jref(y,i,ct->numofbases)][iref(y,i,ct->numofbases)]|NOPAIR;
  }
  for (i=(ct->numofbases)+x+1;i<=(ct->numofbases)+y-1;++i) {
    v[jref(y,i,ct->numofbases)][iref(y,i,ct->numofbases)] = v[jref(y,i,ct->numofbases)][iref(y,i,ct->numofbases)]|NOPAIR;
  }
  for (i=x+1;i<=y-1;++i) {
    v[jref(i,x+ct->numofbases,ct->numofbases)][iref(i,x+ct->numofbases,ct->numofbases)] = v[jref(i,x+ct->numofbases,ct->numofbases)][iref(i,x+ct->numofbases,ct->numofbases)]|NOPAIR;
  }
  for (i=y+1;i<=ct->numofbases;++i) {
    v[jref(i,x+ct->numofbases,ct->numofbases)][iref(i,x+ct->numofbases,ct->numofbases)] = v[jref(i,x+ct->numofbases,ct->numofbases)][iref(i,x+ct->numofbases,ct->numofbases)]|NOPAIR;
  }
  for (i=1;i<=x-1;++i) {
    for (j = x+1;j<=y-1;++j){
      v[jref(i,j,ct->numofbases)][iref(i,j,ct->numofbases)] = v[jref(i,j,ct->numofbases)][iref(i,j,ct->numofbases)] |NOPAIR;
    }
  }
  for (i=x+1;i<=y-1;++i) {
    for (j=y+1;j<=(ct->numofbases)+x-1;++j) {
      v[jref(i,j,ct->numofbases)][iref(i,j,ct->numofbases)] = v[jref(i,j,ct->numofbases)][iref(i,j,ct->numofbases)]|NOPAIR;
    }
  }
  for (i=y+1;i<=ct->numofbases;++i) {
    for (j=(ct->numofbases)+x+1;j<=(ct->numofbases)+y-1;++j) {
      v[jref(i,j,ct->numofbases)][iref(i,j,ct->numofbases)] = v[jref(i,j,ct->numofbases)][iref(i,j,ct->numofbases)]|NOPAIR;
    }
  }
}

void dynforcepairg(int x,structure *ct,char **fce) {
  int i,j;


  for (j=x+1;j<x+ct->numofbases;++j) {
    if (ct->numseq[j]!=3) fce[jref(x,j,ct->numofbases)][iref(x,j,ct->numofbases)]=fce[jref(x,j,ct->numofbases)][iref(x,j,ct->numofbases)]|NOPAIR;
  }

  for (j=x+ct->numofbases+1;j<2*ct->numofbases;++j) {
    if (ct->numseq[j]!=3) fce[jref(x+ct->numofbases,j,ct->numofbases)][iref(x+ct->numofbases,j,ct->numofbases)]=fce[jref(x+ct->numofbases,j,ct->numofbases)][iref(x+ct->numofbases,j,ct->numofbases)]|NOPAIR;
  }

  for (i=x-1;i>0;--i) {
    if (ct->numseq[i]!=3) fce[jref(i,x,ct->numofbases)][iref(i,x,ct->numofbases)]=fce[jref(i,x,ct->numofbases)][iref(i,x,ct->numofbases)]|NOPAIR;
  }

  for (i=x+ct->numofbases-1;i>x;--i) {
    if (ct->numseq[i]!=3) fce[jref(i,x+ct->numofbases,ct->numofbases)][iref(i,x+ct->numofbases,ct->numofbases)]=fce[jref(i,x+ct->numofbases,ct->numofbases)][iref(i,x+ct->numofbases,ct->numofbases)]|NOPAIR;
  }
}

void dynalignforce(structure *ct1, structure *ct2/*, short ****fce*/, char **fce1, char **fce2, bool *dbl1, 
                   bool *dbl2, bool *mod1, bool *mod2) {
  int count;
  //fill the fce array for a dynalign calculation using the following key:

  //SINGLE applies to any fce(i,j) s.t. i or j should be single stranded
  //PAIR applies to any fce(i,j) where i is paired to j
  //NOPAIR applies to any fce(i,j) where either i or j is paired to
  //    another nucleotide or i and j are forbidden to pair
  //DUBLE applies to any fce(i.j) where an nuc, k, i<k<j is double stranded
  //INTER applies to any fce(i,j) such that some k, i<k<j, is the virtual linker
  //    used for intermolecular folding

  //refer to fce[j][i][a][b]
  //where a = k-i+maxsep and i+maxsep>=k>=i-maxsep and i<j
  //and for l<N and j<N2 b = l-j+maxsep and j+maxsep>=l>=j-maxsep
  //for l>N and j>N2, b = l-j+maxsep+N-N2 so that maxsep is "reset" beyond the ends of the sequences when they are of different lengths 


  //start with unpaired nucs in seq1
  for (count=1;count<=ct1->nnopair;++count) {
    //dynalignforceunpaired1(ct1, ct2, fce, maxseparation, ct1->nopair[count]);
    dynalignfceunpaired(ct1,fce1,ct1->nopair[count]);
                
  }

  //now encode unpaired nucs in sequence 2
  for (count=1;count<=ct2->nnopair;++count) {
    //dynalignforceunpaired2(ct1,ct2,fce,maxseparation,ct2->nopair[count]);
    dynalignfceunpaired(ct2,fce2,ct2->nopair[count]);

  }

  //force nucleotides double-stranded in seq1
  for(count=1;count<=ct1->ndbl;++count) {
    dynforcedbl(ct1->dbl[count],ct1,fce1,dbl1);
  }

  //force nucleotides double-stranded in seq2
  for(count=1;count<=ct2->ndbl;++count) {
    dynforcedbl(ct2->dbl[count],ct2,fce2,dbl2);
  }


  //now handle pairs that are required:
  for (count = 1; count <=ct1->npair; ++count) {
                
    dynforcepair(ct1->pair[count][0],ct1->pair[count][1],ct1,fce1);
    dynforcedbl(ct1->pair[count][0],ct1,fce1,dbl1);
    dynforcedbl(ct1->pair[count][1],ct1,fce1,dbl1);


  }

  for (count = 1; count <=ct2->npair;++count) {
    dynforcepair(ct2->pair[count][0],ct2->pair[count][1],ct2,fce2);
    dynforcedbl(ct2->pair[count][0],ct2,fce2,dbl2);
    dynforcedbl(ct2->pair[count][1],ct2,fce2,dbl2);

  }


  //now handle FMN cleavage (U in GU pair)
  for (count=0;count<ct1->ngu;++count) {
    dynforcedbl(ct1->gu[count],ct1,fce1,dbl1);
    dynforcepairg(ct1->gu[count],ct1,fce1);

  }

  for (count=0;count<ct2->ngu;++count) {
    dynforcedbl(ct2->gu[count],ct2,fce2,dbl2);
    dynforcepairg(ct2->gu[count],ct2,fce2);

  }

  //now handle prohibited base pairs
  for (count=0;count<ct1->nforbid;++count) {
    fce1[jref(ct1->forbid[count][0],ct1->forbid[count][1],ct1->numofbases)][iref(ct1->forbid[count][0],ct1->forbid[count][1],ct1->numofbases)]=
      fce1[jref(ct1->forbid[count][0],ct1->forbid[count][1],ct1->numofbases)][iref(ct1->forbid[count][0],ct1->forbid[count][1],ct1->numofbases)]|NOPAIR;
                
    fce1[jref(ct1->forbid[count][1],ct1->forbid[count][0]+ct1->numofbases,ct1->numofbases)][iref(ct1->forbid[count][1],ct1->forbid[count][0]+ct1->numofbases,ct1->numofbases)]=
      fce1[jref(ct1->forbid[count][1],ct1->forbid[count][0]+ct1->numofbases,ct1->numofbases)][iref(ct1->forbid[count][1],ct1->forbid[count][0]+ct1->numofbases,ct1->numofbases)]|NOPAIR;
  }

  for (count=0;count<ct2->nforbid;++count) {
    fce2[jref(ct2->forbid[count][0],ct2->forbid[count][1],ct2->numofbases)][iref(ct2->forbid[count][0],ct2->forbid[count][1],ct2->numofbases)]=
      fce1[jref(ct2->forbid[count][0],ct2->forbid[count][1],ct2->numofbases)][iref(ct2->forbid[count][0],ct2->forbid[count][1],ct2->numofbases)]|NOPAIR;
                
    fce2[jref(ct2->forbid[count][1],ct2->forbid[count][0]+ct2->numofbases,ct2->numofbases)][iref(ct2->forbid[count][1],ct2->forbid[count][0]+ct2->numofbases,ct2->numofbases)]=
      fce1[jref(ct2->forbid[count][1],ct2->forbid[count][0]+ct2->numofbases,ct2->numofbases)][iref(ct2->forbid[count][1],ct2->forbid[count][0]+ct2->numofbases,ct2->numofbases)]|NOPAIR;
  }

  //now handle chemical modification
  for (count=1;count<=ct1->nmod;++count) {

    if (ct1->mod[count]!=1&&ct1->mod[count]!=ct1->numofbases) {
      mod1[ct1->mod[count]]=true;
      mod1[ct1->mod[count]+ct1->numofbases]=true;
    }
  }
  for (count=1;count<=ct2->nmod;++count) {

    if (ct2->mod[count]!=1&&ct2->mod[count]!=ct2->numofbases) {
      mod2[ct2->mod[count]]=true;
      mod2[ct2->mod[count]+ct2->numofbases]=true;
    }
  }
}


//maxseparation is the M parameter that limits the alignment space
// maxseparation < 0 implies that the alignment constraints are instead coming from allowed_alignments
//return an int that indicates whether an error occurred.  0 = no error
int dynalign(structure *ct1, structure *ct2, short **alignment,
              short int maxseparation, short int gapincrease, datatable *data, 
              bool singleinsert, short maxtracebacks, short window, short awindow, short percentsort, short **forcealign,
              bool **allowed_alignments, TProgressDialog *progress, const char *Savefile, bool energyonly, 
			  bool local, bool force, short int numProcessors) {
  
  short int i,j,k,l,a,b,c,d,/*maxsep,*/gap,N,N2,Ndiff;
  short ip,kp,jp,lp,cp,dp;
  integersize en1,lowest;
  short int I;
  int flag;
  bool ikincrement,jldecrement,jldecrement2;
  integersize imin, jmax, kmin, lmax, wval;
  int error;

#ifdef DYNALIGN_SMP
  bool removeprogress;
#endif
  
  
  dynalignarray *w;
  varray *v;
        
  wendarray *w5,*w3;
        
  short int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
                       {0,1,0,1,0,0},{0,0,0,0,0,0}};
        
  dynalignarray *vmod;
        

  short *lowend,*highend;//store the limits calculated from lowlimit and highlimit

  // Select the version of edangle to use
  
  short int (*edangle5)(int i, int j, int ip, structure *ct, datatable *data) = edangle5noforce;
  short int (*edangle3)(int i, int j, int ip, structure *ct, datatable *data) = edangle3noforce;

  if (force) {
    edangle5 = &edangle5force;
    edangle3 = &edangle3force;
  }

  // Local symbols edangle3 and edangle5 are now function pointers to
  // actual functions edangle3noforce and edangle5noforce or
  // edangle3force and edangle5force, respectively, if force is true
  
#ifdef timer
#include <time.h>
  ofstream timeout;
  int seconds;
  char timerstring[100];
  char timelength[10];
  strcpy(timerstring,"time_pf_");
  sprintf(timelength,"%i",ct1->numofbases);
  strcat(timerstring,timelength);
  strcat(timerstring,".out");

  timeout.open(timerstring);
  timeout<<time(NULL)<<"\n";
  seconds = time(NULL);
#endif //timer

  
  //here is a 4-d array to store data about folding constraints that involve two sequences
  //dynalignarray *fce;
  //short ****fce;

  //fce1 and fce2 contain configuration constraints for single sequences
  char **fce1;
  char **fce2;

  //dbl1 and dbl2 flag nucleotides that must be double-stranded in seq1 and seq2, respectively
  //false means a nucleotide is not forced double
  //true indicates a nucleotide is forced double
  bool *dbl1,*dbl2;

  //mod indicates whether a nucleotide is chemically modified
  //mod1 for seq1 and mod2 or seq2
  bool *mod1,*mod2;
        

  //the following are used for chemical modfication cases
  bool modification;
  bool alignmentforced;


  short crit;
  crit = DYNALIGN_INFINITY;

  //store the number of bases in ct1 and ct2 in register shorts
  N = ct1->numofbases;//the length of sequence 1
  N2 = ct2->numofbases;//the length of sequence 2
  Ndiff = N-N2;//the difference in sequence lengths
  //maxsep = maxseparation;

  

  lowest = DYNALIGN_INFINITY;

  
  
        
  
  //double up the sequences as part of suboptimal traceback support:
  for (i=1;i<=N;i++) {
    ct1->numseq[(N)+i] = ct1->numseq[i];
  }
  for (i=1;i<ct2->numofbases;i++) ct2->numseq[ct2->numofbases+i]=ct2->numseq[i];

	
  //fill the low and highend arrays:
  //allocate space:
  lowend = new short [2*N];
  highend = new short [2*N];

  if (maxseparation>=0) {
	//Using the traditional M parameter to constrain the alignment space
	for (i=0;i<2*N;++i) {
		lowend[i]=lowlimit(i,maxseparation,N,N2);
		highend[i]=highlimit(i,maxseparation,N,N2);

		//printf("%d -> (%d, %d)\n", i, lowend[i], highend[i]);
	}
  }
  else {
	//allowed_alignments must be defined to constrain the alignment space
    //need to determine lowend and highend
	  for (i=0;i<2*N;++i) {
		lowend[i]=lowlimit(i,allowed_alignments,N,N2);
		highend[i]=highlimit(i,allowed_alignments,N,N2);

		//printf("%d -> (%d, %d)\n", i, lowend[i], highend[i]);
	  }


  }

  //Allocate the 2-d and 4-d arrays for storing mfe's for subfragments:
  w5 = new wendarray(N,N2,lowend,highend);
  w3 = new wendarray(N,N2,lowend,highend);

  v = new varray(N,N2,lowend,highend,ct1->tem,energyonly);
  w = new dynalignarray(N,N2,lowend,highend,energyonly);
  
  
	
  if (force) {
    
    //allocate the fce array for constraints
    //The following definitions are bitwise applied to the fce array:
    //SINGLE applies to any fce(i,j) s.t. i or j should be single stranded
    //PAIR applies to any fce(i,j) where i is paired to j
    //NOPAIR applies to any fce(i,j) where either i or j is paired to
    //  another nucleotide or i and j are forbidden to pair
    //DUBLE applies to any fce(i.j) where an nuc, k, i<k<j is double stranded
    //INTER applies to any fce(i,j) such that some k, i<k<j, is the virtual linker
    //  used for intermolecular folding
    //INTER applies to any fce(i,i) s.t. i is the center of the virtual linker
    //The above terms are defined in define.h

    //note that INTER is not currently supported in Dynalign
    //fce = new dynalignarray(N,N2,maxseparation);


    //allocate space for single sequence configuration in fce1 and fce2

    fce1 = new char *[2*N+1];

    for (i=0;i<=2*N;++i) {
      if (i<=N) {
        fce1[i] = new char [i+1];
        for (j=0;j<=i;++j) fce1[i][j]=0;
      }
      else {
        fce1[i] = new char [2*N-i+1];
        for (j=0;j<=2*N-i;++j) fce1[i][j]=0;
      }
    }

    fce2 = new char *[2*N2+1];

    for (i=0;i<=2*N2;++i) {
      if (i<=N2) {
        fce2[i] = new char [i+1];
        for (j=0;j<=i;++j) fce2[i][j]=0;
      }
      else {
        fce2[i] = new char [2*N2-i+1];
        for (j=0;j<=2*N2-i;++j) fce2[i][j]=0;
      }
    }

    dbl1 = new bool [2*N];
    mod1=new bool [2*N];
    for (i=0;i<2*N;++i) {
      dbl1[i]=false;
      mod1[i]=false;
    }
    dbl2 = new bool [2*N2];
    mod2 = new bool [2*N2];
    for (i=0;i<2*N2;++i) {
      dbl2[i]=false;
      mod2[i]=false;
    }

    //store the locations of dbl1 and dbl2 in ct1 and ct2, respectively
    //this facilitates checks done in functions edangle5 and edangle3
    ct1->fcedbl=dbl1;
    ct2->fcedbl=dbl2;

    //now assign values to fce according to the above key using the function dynalignforce
    dynalignforce(ct1,ct2/*,fce*/,fce1,fce2,dbl1,dbl2,mod1,mod2);

    if (ct1->nmod>0||ct2->nmod>0) {
      modification = true;
      //For chemical modification, a second v array, vmod is required
      vmod = new dynalignarray(N,N2,lowend,highend);
                
    }
    else {
      modification =false;
      vmod = NULL;
    }

    alignmentforced = (forcealign != NULL);

  } else {
    vmod = NULL;
	dbl1=NULL;
	dbl2=NULL;
	alignmentforced=false;
	fce1=NULL;fce2=NULL;
	modification=false;
	mod1=NULL;
	mod2=NULL;
  }

  //maxsep = maxseparation; //place maxseparation into a register int
  gap = gapincrease;//place the gap penalty in a register short
        

  // look for alignment of every i paired to j in sequence 1 to k
  // paired to l in sequence 2
  
  //for (size=minloop;size<=ct1->numofbases;size++) {

#ifndef DYNALIGN_SMP
  for (j = minloop; j <= N; ++j) {
    if (progress != NULL) {
      if (!energyonly) progress->update((100 * j) / (2*N));
	  else progress->update((100 * j) / (N));
    }
    
    for (i = min(N,j-1); i >= 1; --i) {
      
      kmin = max(lowend[i],1);
      
      for (k = min(highend[i], N2); k >= kmin; k--) {
        
        lmax = min(highend[j], ct2->numofbases);
        
        for (l = max(max(lowend[j],1),k+minloop); l <= lmax; ++l) {
          
          if (((j - i) > minloop) && ((l - k) > minloop)) {
            dynalignstep(ct1, ct2, data,
                         v, w, w5, w3, 
                         vmod, mod1, mod2, modification,
                         fce1, fce2,
                         alignmentforced, forcealign,
                         dbl1, dbl2,
                         i, j, k, l, N, N2,
                         lowend, highend, gap, singleinsert,
                         force,local,
                         edangle5, edangle3);
          }
        }
      }
    }
  }
#else //DYNALIGN_SMP defined
  pthread_t consumerids[numProcessors];
  pthread_t producerid;


  //If progress is the NULL pointer, go ahead and allocate it here...
  if (progress==NULL) {
	  progress = dynamic_cast<TProgressDialog*> (new ObservingTextProgressBar());
	  removeprogress = true;
  }
  else removeprogress = false;
   dynamic_cast<ObservingTextProgressBar*>(progress)
    ->setMaximum(getNumInternalRanks(N, N2) +
                 getNumExternalRanks(N, N2));
  
  RankManager internalManager(getNumInternalRanks(N, N2),
                              numProcessors);

  internalManager.subscribe(dynamic_cast<ObservingTextProgressBar*>(progress));

  rankconsumerargs rca = {&internalManager, data, vmod, w, ct1, ct2, v, w5, w3,
                          mod1, mod2, dbl1, dbl2, lowend, highend,
                          fce1, fce2, forcealign,
                          edangle5, edangle3,
                          gap, modification,
                          alignmentforced, singleinsert, force,local};
  rankproducerargs rpa = {&internalManager, lowend, highend,
                          N, N2, minloop};


  for (int i = 0; i < numProcessors; i++) {
    pthread_create(&(consumerids[i]), NULL, &rankconsumer, static_cast<void *>(&rca));
  }

  pthread_create(&producerid, NULL, &internalrankproducer, static_cast<void *>(&rpa));
  pthread_join(producerid, NULL);

  for (int i = 0; i < numProcessors; i++) {
    pthread_join(consumerids[i], NULL);
  }
#endif //DYNALIGN_SMP

  //calculate w3 and w5
  //w3[c][d] = lowest free energy of fragment c->N and d->N2
  //w5[c][d] = lowest free energy of fragment 1->c and 1->d
  //initialize w5[i][k]:
  if (local) {
    for (ip=0;ip<=N;++ip){
                  
      //kp = highend[ip];
      //kp = lowend[j];
      //kp = lowend[ip];
      for (kp=lowend[ip]/*0*/;kp<=highend[ip]/*<2*maxsep+2*/;++kp){
        w5->f(ip,kp)=0;//w5[ip][kp] = 0;
      }
    }
  } else {
							
    for (ip=0;ip<=N;ip++){	
		
		for (kp = max(0, lowend[ip]); kp <= min(highend[ip], N2); ++kp) {
			
			w5->f(ip,kp) = gap* (abs(ip - kp));
		}
	}
  }
							
  for (ip=1;ip<=N;++ip) {
    cp = min(N2,highend[ip]);//min(N2,ip+maxsep);
    for (kp=max(1,lowend[ip])/*max(1,ip-maxsep)*/;kp<=cp;++kp) {
      //ap = kp-ip+maxsep;

      //although this is really a decrement...
      ikincrement = kp-1 <= highend[ip-1] && kp-1 >= lowend[ip-1];

      if (ikincrement &&
          (!force ||
           (!dbl1[ip] && !dbl2[kp]))) {
        en1=w5->f(ip-1,kp-1);//w5[ip-1][ap]; //adding one nuc on both that is unpaired and unstacked
      } else {
        en1 = DYNALIGN_INFINITY;
      }

      if (kp<=highend[ip-1] && kp>=lowend[ip-1] &&
          (!force ||
           (!dbl1[ip]))) {
        en1 = min(en1,w5->f(ip-1,kp)/*w5[ip-1][ap+1]*/+gap);//add a gapped nuc, ip
      }
                  
      if (kp-1>=lowend[ip] && kp-1<=highend[ip] &&
          (!force ||
           (!dbl2[kp]))) {
        en1 = min(en1,w5->f(ip,kp-1)/*w5[ip][ap-1]*/+gap);//add a gapped nuc, kp (aka ap)
      }
                  
      for (jp=0;jp+minloop<ip;++jp) {
        dp = min(highend[jp],N2);
        dp = min(dp,kp-minloop-1);
        for (lp=max(0,lowend[jp]);lp<=dp;++lp) {
          	

          //check whether mine[i][a] is split so that the lowest free energy is
          //an exterior fragment from 1 to j and a helix from j+1 to i; 

          //check all possible alignments (in index a(k) and b(l))
					
          //must consider whether an unpaired nucleotide is stacked
          //is stacked onto each of the nucleotides in a helix
          //There are 16 cases:
          //consider them in this order (0 unstacked, 1 stacked)
          //		j+1	i	l+1	k
          //1		0	0	0	0
          //2		0	0	0	1
          //3		0	0	1	0
          //4		0	0	1	1
          //5		0	1	0	0
          //6		0	1	0	1
          //7		0	1	1	0
          //8		0	1	1	1
          //9		1	0	0	0
          //10		1	0	0	1
          //11		1	0	1	0
          //12		1	0	1	1
          //13		1	1	0	0
          //14		1	1	0	1
          //15		1	1	1	0
          //16		1	1	1	1

          //note that for these exterior loops:
          //j<i
          //l<k
          
          //although this is an increment...
          jldecrement  = lp+1 >= lowend[jp+1] && lp+1 <= highend[jp+1];
          //although this is an increment...
          jldecrement2 = lp+2 >= lowend[jp+2] && lp+2 <= highend[jp+2];


          //no stacking
          //case 1: 0	0	0	0
          if(jldecrement) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip,lp+1,kp)/*v[ip][jp+1][bp][ap]*/+
                                    penalty(ip,jp+1,ct1,data)+penalty(kp,lp+1,ct2,data));

          //case 6: 0	1	0	1


          if(ikincrement&&jldecrement) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip-1,lp+1,kp-1)/*v[ip-1][jp+1][bp][ap]*/+edangle3(ip-1,jp+1,ip,ct1,data)+
                                                 edangle3(kp-1,lp+1,kp,ct2,data)+
                                                 penalty(ip-1,jp+1,ct1,data)+penalty(kp-1,lp+1,ct2,data));

          //case 11	1	0	1	0
          if (jldecrement2) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip,lp+2,kp)/*v[ip][jp+2][bp][ap]*/+edangle5(jp+2,ip,jp+1,ct1,data)+
                                      edangle5(lp+2,kp,lp+1,ct2,data)+penalty(ip,jp+2,ct1,data)+penalty(lp+2,kp,ct2,data));


          //case 16	1	1	1	1
          if (ikincrement&&jldecrement2) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip-1,lp+2,kp-1)/*v[ip-1][jp+2][bp][ap]*/+edangle5(jp+2,ip-1,jp+1,ct1,data)+
                                                   edangle3(ip-1,jp+2,ip,ct1,data)+edangle5(lp+2,kp-1,lp+1,ct2,data)+
                                                   edangle3(kp-1,lp+2,kp,ct2,data)+penalty(ip-1,jp+2,ct1,data)+penalty(lp+2,kp-1,ct2,data));

					
					
          if (kp-1>=lowend[ip]) {
            //case 2: 0	0	0	1
            if (jldecrement) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip,lp+1,kp-1)/*v[ip][jp+1][bp][ap-1]*/+edangle3(kp-1,lp+1,kp,ct2,data)+
                                       penalty(ip,jp+1,ct1,data)+penalty(kp-1,lp+1,ct2,data)+gap);

            //case 12	1	0	1	1
            if (jldecrement2) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip,lp+2,kp-1)/*v[ip][jp+2][bp][ap-1]*/+edangle5(jp+2,ip,jp+1,ct1,data)+
                                        edangle5(lp+2,kp-1,lp+1,ct2,data) + edangle3(kp-1,lp+2,kp,ct2,data)+
                                        penalty(ip,jp+2,ct1,data)+penalty(kp-1,lp+2,ct2,data)+gap);

            if (lp+2<=highend[jp+1]/*bp+1<2*maxsep+2*/&&lp+2>=lowend[jp+1]) {
              //case 4
              en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip,lp+2,kp-1)/*v[ip][jp+1][bp+1][ap-1]*/+edangle3(kp-1,lp+2,kp,ct2,data)+
                        edangle5(lp+2,kp-1,lp+1,ct2,data)+penalty(ip,jp+1,ct1,data)+penalty(kp-1,lp+2,ct2,data)+2*gap);
            }
          }
					
          if (lp+2<=highend[jp+1]&&lp+2>=lowend[jp+1]) {
            //case 3: 0	0	1	0
            en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip,lp+2,kp)/*v[ip][jp+1][bp+1][ap]*/+edangle5(lp+2,kp,lp+1,ct2,data)+
                      penalty(ip,jp+1,ct1,data)+penalty(lp+2,kp,ct2,data)+gap);

            //case 8:		0	1	1	1
            if(ikincrement) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip-1,lp+2,kp-1)/*v[ip-1][jp+1][bp+1][ap]*/+edangle3(ip-1,jp+1,ip,ct1,data)+
                                      edangle5(lp+2,kp-1,lp+1,ct2,data)+edangle3(kp-1,lp+2,kp,ct2,data)+
                                      penalty(ip-1,jp+1,ct1,data)+penalty(lp+2,kp-1,ct2,data)+gap);

            if (kp<=highend[ip-1]&&kp>=lowend[ip-1]) {
              //case 7		0	1	1	0
              en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip-1,lp+2,kp)/*v[ip-1][jp+1][bp+1][ap+1]*/+edangle5(lp+2,kp,lp+1,ct2,data)+
                        edangle3(ip-1,jp+1,ip,ct1,data)+penalty(ip-1,jp+1,ct1,data)+penalty(lp+2,kp,ct2,data)+2*gap);

            }
          }

          if (kp<=highend[ip-1]&&kp>=lowend[ip-1]) {
            //case5: 0	1	0	0
            if (jldecrement) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip-1,lp+1,kp)/*v[ip-1][jp+1][bp][ap+1]*/+edangle3(ip-1,jp+1,ip,ct1,data)+
                                       penalty(ip-1,jp+1,ct1,data)+penalty(kp,lp+1,ct2,data)+gap);

            //case 15	1	1	1	0
            if (jldecrement2) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip-1,lp+2,kp)/*v[ip-1][jp+2][bp][ap+1]*/ + edangle5(jp+2,ip-1,jp+1,ct1,data)+
                                        edangle3(ip-1,jp+2,ip,ct1,data)+edangle5(lp+2,kp,lp+1,ct2,data)+
                                        penalty(ip-1,jp+2,ct1,data)+penalty(lp+2,kp,ct2,data)+gap);

            if (lp+1>=lowend[jp+2]&&lp+1<=highend[jp+2]/*bp>0*/) {
              //case 13	1	1	0	0
              en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip-1,lp+1,kp)/*v[ip-1][jp+2][bp-1][ap+1]*/+edangle5(jp+2,ip-1,jp+1,ct1,data)+
                        edangle3(ip-1,jp+2,ip,ct1,data)+
                        penalty(ip-1,jp+2,ct1,data)+penalty(kp,lp+1,ct2,data)+2*gap);

            }


          }
          if (lp+1>=lowend[jp+2]&&lp+1<=highend[jp+2]/*bp>0*/) {
            //case 9		1	0	0	0
            en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip,lp+1,kp)/*v[ip][jp+2][bp-1][ap]*/+edangle5(jp+2,ip,jp+1,ct1,data)+
                      penalty(ip,jp+2,ct1,data)+penalty(kp,lp+1,ct2,data)+gap);

					
            //case 14	1	1	0	1
            if (ikincrement) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip-1,lp+1,kp-1)/*v[ip-1][jp+2][bp-1][ap]*/+edangle5(jp+2,ip-1,jp+1,ct1,data)+
                                       edangle3(ip-1,jp+2,ip,ct1,data)+edangle3(kp-1,lp+1,kp,ct2,data)+
                                       penalty(ip-1,jp+2,ct1,data)+penalty(kp-1,lp+1,ct2,data)+gap);

            if (kp-1>=lowend[ip]/*ap>0*/) {
              //case 10	1	0	0	1
              en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip,lp+1,kp-1)/*v[ip][jp+2][bp-1][ap-1]*/+edangle5(jp+2,ip,jp+1,ct1,data)+
                        edangle3(kp-1,lp+1,kp,ct2,data)+penalty(ip,jp+2,ct1,data)+penalty(kp-1,lp+1,ct2,data)+2*gap);
            }
          }
        }
      }
      w5->f(ip,kp) /*w5[ip][ap]*/ = min(en1,w5->f(ip,kp));
      if (local) {
        if (en1<lowest) {
          lowest = en1;
        }
      } else if (/*ip==N&&*/abs((N-ip)-(N2-kp))*gap+en1/*en1+gap*abs(ap-maxsep)*/<lowest) {
        lowest = abs((N-ip)-(N2-kp))*gap+en1;
      }
    }
  }
							
  if (!energyonly) { // calculate w3 if proceeding with extra
                     // calculations beyond energy
    //initialize w3[i][k]:
    if (local) {
      //local alignment calculation
	  //for (kp=N2;kp>=0;kp--){
      //    w3->f(N2+1,kp) = 0;
      //}
		for (kp=N2+1;kp>=0;kp--){
			w3->f(N+1,kp) = 0;
		}
		//for (ip=N+1;ip>0;ip--){//Unallocated space!
		//	w3->f(ip,N2+1) = 0;
		//}
      for (ip=N;ip>0;ip--){
		
        for (kp=highend[ip];kp>=lowend[ip]/*0*/;kp--){
          w3->f(ip,kp)=0;
        }
      }
    }
    else {
      //global alignment calculation
		for (kp=N2+1;kp>=0;kp--){
			w3->f(N+1,kp) = gap*(abs(N2-kp+1));
		}
		//for (ip=N+1;ip>0;ip--){ //Unallocated space!
		//	w3->f(ip,N2+1) = gap*(abs(N-ip+1));
		//}
      for (ip=N;ip>0;ip--){
        for (kp=highend[ip];kp>=lowend[ip];kp--){
          w3->f(ip,kp) = gap*(abs((N-ip)-(N2-kp)));
        }
      }
    }
							
    for (ip=N;ip>=1;ip--) {
      cp = max(1,lowend[ip]);
      for (kp=min(N2,highend[ip]);kp>=cp;--kp) {
        
        ikincrement = kp+1 >= lowend[ip+1] && kp+1 <= highend[ip+1];

        if(ikincrement &&
           (!force ||
            (!dbl1[ip] && !dbl2[kp]))) {
          en1 = w3->f(ip+1,kp+1);//adding one nuc on both that is unpaired and unstacked
        } else {
          en1 = DYNALIGN_INFINITY;
        }
									
        if (kp>=lowend[ip+1]&&kp<=highend[ip+1] &&
            (!force ||
             (!dbl1[ip]))) {
          en1 = min(en1,w3->f(ip+1,kp)/*w3[ip+1][ap-1]*/+gap);//add a gapped nuc, ip
        }

        if (kp+1<=highend[ip] /*ap+1<2*maxsep+2*/ &&
            (!force ||
             (!dbl2[kp]))) {
          en1 = min(en1,w3->f(ip,kp+1)/*w3[ip][ap+1]*/+gap);//add a gapped nuc, kp (aka ap)
        }
			
        for (jp=N+1;jp>=ip+minloop;jp--) {
          dp = max(lowend[jp]/*jp-maxsep-1*/,kp+minloop);
          for (lp=min(N2,highend[jp]/*jp+maxsep*/);lp>=dp;lp--) {
            //bp = lp-jp+maxsep;	

            //check whether mine[i][a] is split so that the lowest free energy is
            //an exterior fragment from 1 to j and a helix from j+1 to i; 

            //check all possible alignments (in index a(k) and b(l))
					
            //must consider whether an unpaired nucleotide is stacked
            //is stacked onto each of the nucleotides in a helix
            //There are 16 cases:
            //consider them in this order (0 unstacked, 1 stacked)
            //		jp+1	ip	lp+1	k
            //1		0		0	0		0
            //2		0		0	0		1
            //3		0		0	1		0
            //4		0		0	1		1
            //5		0		1	0		0
            //6		0		1	0		1
            //7		0		1	1		0
            //8		0		1	1		1
            //9		1		0	0		0
            //10	1		0	0		1
            //11	1		0	1		0
            //12	1		0	1		1
            //13	1		1	0		0
            //14	1		1	0		1
            //15	1		1	1		0
            //16	1		1	1		1

			  if (lp+1>=lowend[jp+1]) {
				  
				  if ((lp+1<=highend[jp+1])) {
						if (lp+1==N2+1) {
							//This is a case where there is no extension on sequence 2, so handle that correctly:
							if (local) wval = 0;
							else {
								wval = gap*(N-jp);

							}

						}
						else wval=w3->f(jp+1,lp+1);
				  }
				  
				  else wval = DYNALIGN_INFINITY;//In the future, this could probably be a condition around these cases...


				  jldecrement = lp-1 <= highend[jp-1] && lp-1 >= lowend[jp-1];

				  //no stacking
				  //case 1: 0		0	0		0
				  en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp,kp,lp)/*v[jp][ip][ap][bp]*/+
							penalty(ip,jp,ct1,data)+penalty(kp,lp,ct2,data));

				  //case 6: 0		1	0		1

						
				  if (ikincrement) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp,kp+1,lp)/*v[jp][ip+1][ap][bp]*/+edangle5(ip+1,jp,ip,ct1,data)+
											 edangle5(kp+1,lp,kp,ct2,data)+
											 penalty(ip+1,jp,ct1,data)+penalty(kp+1,lp,ct2,data));

				  //case 11	1	0	1	0
				  if (jldecrement) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp-1,kp,lp-1)/*v[jp-1][ip][ap][bp]*/+edangle3(jp-1,ip,jp,ct1,data)+
											 edangle3(lp-1,kp,lp,ct2,data)+penalty(ip,jp-1,ct1,data)+penalty(lp-1,kp,ct2,data));


				  //case 16	1	1	1	1
				  if(ikincrement&&jldecrement) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp+1,lp-1)/*v[jp-1][ip+1][ap][bp]*/+edangle3(jp-1,ip+1,jp,ct1,data)+
														 edangle5(ip+1,jp-1,ip,ct1,data)+edangle3(lp-1,kp+1,lp,ct2,data)+
														 edangle5(kp+1,lp-1,kp,ct2,data)+penalty(ip+1,jp-1,ct1,data)+penalty(lp-1,kp+1,ct2,data));

						

						
				  if (kp+1<=highend[ip]/*ap+1<2*maxsep+2*/) {
					//case 2: 0		0	0		1
					en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp,kp+1,lp)/*v[jp][ip][ap+1][bp]*/+edangle5(kp+1,lp,kp,ct2,data)+
							  penalty(ip,jp,ct1,data)+penalty(kp+1,lp,ct2,data)+gap);

					//case 12	1	0	1	1
					if (jldecrement) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp-1,kp+1,lp-1)/*v[jp-1][ip][ap+1][bp]*/+edangle3(jp-1,ip,jp,ct1,data)+
											   edangle3(lp-1,kp,lp,ct2,data) + edangle5(kp+1,lp-1,kp,ct2,data)+
											   penalty(ip,jp-1,ct1,data)+penalty(kp+1,lp-1,ct2,data)+gap);

					if (lp-1>=lowend[jp]/*bp>0*/) {
					  //case 4 0		0	1		1
					  en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp,kp+1,lp-1)/*v[jp][ip][ap+1][bp-1]*/+edangle5(kp+1,lp-1,kp,ct2,data)+
								edangle3(lp-1,kp+1,lp,ct2,data)+penalty(ip,jp,ct1,data)+penalty(kp+1,lp-1,ct2,data)+2*gap);

					}

				  }

						
				  if (lp-1>=lowend[jp]/*bp>0*/) {
					//case 3: 0		0	1		0
					en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp,kp,lp-1)/*v[jp][ip][ap][bp-1]*/+edangle3(lp-1,kp,lp,ct2,data)+
							  penalty(ip,jp,ct1,data)+penalty(lp-1,kp,ct2,data)+gap);

					//case 8:		0	1	1	1
					if (ikincrement) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp,kp+1,lp-1)/*v[jp][ip+1][ap][bp-1]*/+edangle5(ip+1,jp,ip,ct1,data)+
											   edangle3(lp-1,kp+1,lp,ct2,data)+edangle5(kp+1,lp-1,kp,ct2,data)+
											   penalty(ip+1,jp,ct1,data)+penalty(lp-1,kp+1,ct2,data)+gap);

					if (kp>=lowend[ip+1]&&kp<=highend[ip+1]) {
					  //case 7		0	1	1	0
					  en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp,kp,lp-1)/*v[jp][ip+1][ap-1][bp-1]*/+edangle3(lp-1,kp,lp,ct2,data)+
								edangle5(ip+1,jp,ip,ct1,data)+penalty(ip+1,jp,ct1,data)+penalty(lp-1,kp,ct2,data)+2*gap);

					}


				  }


				  if (kp>=lowend[ip+1]&&kp<=highend[ip+1]) {
					//case5: 0 1 0 0
					en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp,kp,lp)/*v[jp][ip+1][ap-1][bp]*/+edangle5(ip+1,jp,ip,ct1,data)+
							  penalty(ip+1,jp,ct1,data)+penalty(kp,lp,ct2,data)+gap);

					//case 15	1	1	1	0
					if (jldecrement) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp,lp-1)/*v[jp-1][ip+1][ap-1][bp]*/ + edangle3(jp-1,ip+1,jp,ct1,data)+
											   edangle5(ip+1,jp-1,ip,ct1,data)+edangle3(lp-1,kp,lp,ct2,data)+
											   penalty(ip+1,jp-1,ct1,data)+penalty(lp-1,kp,ct2,data)+gap);

					if (lp<=highend[jp-1]&&lp>=lowend[jp-1]) {
					  //case 13	1	1	0	0
					  en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp,lp)/*v[jp-1][ip+1][ap-1][bp+1]*/+edangle3(jp-1,ip+1,jp,ct1,data)+
								edangle5(ip+1,jp-1,ip,ct1,data)+
								penalty(ip+1,jp-1,ct1,data)+penalty(kp,lp,ct2,data)+2*gap);

					}


				  }
				  if (lp<=highend[jp-1]&&lp>=lowend[jp-1]) {
					//case 9		1	0	0	0
					en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp-1,kp,lp)/*v[jp-1][ip][ap][bp+1]*/+edangle3(jp-1,ip,jp,ct1,data)+
							  penalty(ip,jp-1,ct1,data)+penalty(kp,lp,ct2,data)+gap);

						
					//case 14	1	1	0	1
					if (ikincrement) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp+1,lp)/*v[jp-1][ip+1][ap][bp+1]*/+edangle3(jp-1,ip+1,jp,ct1,data)+
											   edangle5(ip+1,jp-1,ip,ct1,data)+edangle5(kp+1,lp,kp,ct2,data)+
											   penalty(ip+1,jp-1,ct1,data)+penalty(kp+1,lp,ct2,data)+gap);


					if (kp+1<=highend[ip]/*ap<2*maxsep+2*/) {
					  //case 10	1	0	0	1
					  en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp-1,kp+1,lp)/*v[jp-1][ip][ap+1][bp+1]*/+edangle3(jp-1,ip,jp,ct1,data)+
								edangle5(kp+1,lp,kp,ct2,data)+penalty(ip,jp-1,ct1,data)+penalty(kp+1,lp,ct2,data)+2*gap);

					}
				  }
				} 
          }
        }
        w3->f(ip,kp)/*w3[ip][ap]*/ = min(en1,w3->f(ip,kp));
      }
    }

    // w3 is now filled, proceed with extra filling operations on v
    // and w
	
#ifndef DYNALIGN_SMP
    jmax = 2*N - 1;
    for (; j <= jmax; ++j) {
      if (progress != NULL) {
        progress->update((100 * j) / (2*N));
      }
      
      imin = j - N + 1;
      for (i = min(N, j - 1); i >= imin; --i) {

        kmin = max(lowend[i], 1);
        for (k = min(highend[i], N2); k >= kmin; k--) {

          lmax = min(highend[j],  N2+k);//2 * N2 - 1);
          for (l = max(lowend[j], N2+1); l <= lmax; ++l) {
            //fprintf(stderr, "('v', %d, %d, %d, %d),\n", i, j, k, l);
            dynalignstep(ct1, ct2, data,
                         v, w, w5, w3, 
                         vmod, mod1, mod2, modification,
                         fce1, fce2,
                         alignmentforced, forcealign,
                         dbl1, dbl2,
                         i, j, k, l,
                         N, N2, lowend, highend, gap, singleinsert,
                         force,local,
                         edangle5, edangle3);
          }
        }
      }
    }

    if (progress != NULL) {
      progress->update(100);
    }
#else
    RankManager externalManager(getNumExternalRanks(N, N2),
                                numProcessors);
    
    externalManager.subscribe(dynamic_cast<ObservingTextProgressBar*>(progress));

    rca.manager = &externalManager;
    rpa.manager = &externalManager;

    for (int i = 0; i < numProcessors; i++) {
      pthread_create(&(consumerids[i]), NULL, &rankconsumer, static_cast<void *>(&rca));
    }

    pthread_create(&producerid, NULL, &externalrankproducer, static_cast<void *>(&rpa));
    pthread_join(producerid, NULL);

    for (int i = 0; i < numProcessors; i++) {
      pthread_join(consumerids[i], NULL);
    }
#endif
	
  }//if (!energyonly) ...

  if (Savefile!=NULL) {
    ofstream sav(Savefile,ios::binary);
	
		
    //write the save file information so that the sequence can be re-folded,
    //include thermodynamic data to prevent traceback errors

    //write a flag as to whether vmod will be needed
    //also include whether the data needed for suboptimal structure prediction has been included
    //or not.  For backwards compatibility, use one digit.
    //0 - suboptimal, no modifications
    //1 - suboptimal, with modifications
    //2 - optimal only, no modifications
    //3 - optimal only, with modifications

    if (ct1->nmod>0||ct2->nmod>0) {
      if (energyonly) flag=3;
      else flag = 1;

    }
    else {
      if (energyonly) flag = 2;
      else flag = 0;

    }
    write(&sav,&flag);


    //write the 3 parameters needed for array allocation
    write(&sav,&(ct1->numofbases));
    write(&sav,&(ct2->numofbases));
    write(&sav,&maxseparation);

    //start with structure information
		
    write(&sav,&(ct1->npair));
    for (i=0;i<=ct1->npair;++i) {
      write(&sav,&(ct1->pair[i][0]));
      write(&sav,&(ct1->pair[i][1]));
    }
    for (i=0;i<=ct1->numofbases;++i) {
      write(&sav,&(ct1->hnumber[i]));
      sav.write(&(ct1->nucs[i]),1);
    }

    for (i=0;i<=2*ct1->numofbases;i++) {
      write(&sav,&(ct1->numseq[i]));
    }
    
    write(&sav,&(ct1->ndbl));
    for (i=0;i<=ct1->ndbl;i++) {
      write(&sav,&(ct1->dbl[i]));
    }
    
    write(&sav,&(ct1->intermolecular));
    if (ct1->intermolecular) {
      for (i=0;i<3;i++) {
        write(&sav,&(ct1->inter[i]));
      }
    }

    write(&sav,&(ct1->nnopair));
    for (i=0;i<=ct1->nnopair;i++) {
      write(&sav,&(ct1->nopair[i]));
    }
    
    write(&sav,&(ct1->nmod));
    for (i=0;i<=ct1->nmod;i++) {
      write(&sav,&(ct1->mod[i]));
    }
    
    write(&sav,&(ct1->ngu));
    for (i=0;i<=ct1->ngu;i++) {
      write(&sav,&(ct1->gu[i]));
    }
    
    write(&sav,ct1->ctlabel[1]);

    write(&sav,&(ct1->templated));
    if (ct1->templated) {
      for (i=0;i<=ct1->numofbases;i++) {
        for (j=0;j<=i;j++) {
          write(&sav,&(ct1->tem[i][j]));
        }
      }
    }

		
		
    write(&sav,&(ct2->npair));
    for (i=0;i<=ct2->npair;++i) {
      write(&sav,&(ct2->pair[i][0]));
      write(&sav,&(ct2->pair[i][1]));
    }
    
    for (i=0;i<=ct2->numofbases;++i) {
      write(&sav,&(ct2->hnumber[i]));
      sav.write(&(ct2->nucs[i]),1);
    }

    for (i=0;i<=2*ct2->numofbases;i++) {
      write(&sav,&(ct2->numseq[i]));
    }
    
    write(&sav,&(ct2->ndbl));
    for (i=0;i<=ct2->ndbl;i++) {
      write(&sav,&(ct2->dbl[i]));
    }
    
    write(&sav,&(ct2->intermolecular));
    if (ct2->intermolecular) {
      for (i=0;i<3;i++) {
        write(&sav,&(ct2->inter[i]));
      }	
    }

    write(&sav,&(ct2->nnopair));
    for (i=0;i<=ct2->nnopair;i++) {
      write(&sav,&(ct2->nopair[i]));
    }
    
    write(&sav,&(ct2->nmod));
    for (i=0;i<=ct2->nmod;i++) {
      write(&sav,&(ct2->mod[i]));
    }
    
    write(&sav,&(ct2->ngu));
    for (i=0;i<=ct2->ngu;i++) {
      write(&sav,&(ct2->gu[i]));  
    }
    write(&sav,ct2->ctlabel[1]);

    write(&sav,&(ct2->templated));
    if (ct2->templated) {
      for (i=0;i<=ct2->numofbases;i++) {
        for (j=0;j<=i;j++) {
          write(&sav,&(ct2->tem[i][j]));
        }
      }
    }

    //now write the thermodynamic data:
    for (i=0;i<5;++i) {
      write(&sav,&(data->poppen[i]));
    }

    write(&sav,&(data->maxpen));
    for (i=0;i<11;++i) {
      write(&sav,&(data->eparam[i]));
    }
    for (i=0;i<31;++i) {
      write(&sav,&(data->inter[i]));
      write(&sav,&(data->bulge[i]));
      write(&sav,&(data->hairpin[i]));
    }
    
    for (i=0;i<6;++i) {
      for (j=0;j<6;++j) {
        for (k=0;k<6;++k) {
          for (l=0;l<3;++l) {
            write(&sav,&(data->dangle[i][j][k][l]));
          }
          for (l=0;l<6;++l) {
            write(&sav,&(data->stack[i][j][k][l]));
            write(&sav,&(data->tstkh[i][j][k][l]));
            write(&sav,&(data->tstki[i][j][k][l]));
            write(&sav,&(data->coax[i][j][k][l]));
            write(&sav,&(data->tstackcoax[i][j][k][l]));
            write(&sav,&(data->coaxstack[i][j][k][l]));
            write(&sav,&(data->tstack[i][j][k][l]));
            write(&sav,&(data->tstkm[i][j][k][l]));
            write(&sav,&(data->tstki23[i][j][k][l]));
            write(&sav,&(data->tstki1n[i][j][k][l]));
            for (a=0;a<6;++a) {
              for (b=0;b<6;++b) {
                write(&sav,&(data->iloop11[i][j][k][l][a][b]));
                for (c=0;c<6;++c) {
                  if (inc[i][j]&&inc[b][c]) {
                    write(&sav,&(data->iloop21[i][j][k][l][a][b][c]));
                  }
                  for (d=0;d<6;++d) {
                    if (inc[i][k]&&inc[j][l]) {
                      write(&sav,&(data->iloop22[i][j][k][l][a][b][c][d]));
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    
    write(&sav,&(data->numoftloops));
    for (i=0;i<=data->numoftloops;++i) {
      for (j=0;j<2;j++) {
        write(&sav,&(data->tloop[i][j]));
      }
    }

    write(&sav,&(data->numoftriloops));
    for (i=0;i<=data->numoftriloops;++i) {
      for (j=0;j<2;j++) {
        write(&sav,&(data->triloop[i][j]));
      }
    }
    
    write(&sav,&(data->numofhexaloops));
    for (i=0;i<=data->numofhexaloops;++i) {
      for (j=0;j<2;j++) {
        write(&sav,&(data->hexaloop[i][j]));
      }
    }
    
    write(&sav,&(data->auend));
    write(&sav,&(data->gubonus));
    write(&sav,&(data->cint));
    write(&sav,&(data->cslope));
    write(&sav,&(data->c3));
    write(&sav,&(data->efn2a));
    write(&sav,&(data->efn2b));
    write(&sav,&(data->efn2c));
    write(&sav,&(data->init));
    write(&sav,&(data->mlasym));
    write(&sav,&(data->strain));
    write(&sav,&(data->prelog));
    write(&sav,&(data->singlecbulge));
	
    //now write the folding data:
		
    write(&sav,&gap);
    write(&sav,&lowest);
    write(&sav,&singleinsert);


	//write allowed_alignments if maxseparation < 0
	if (maxseparation < 0) {
		for (i=0;i<=ct1->numofbases;++i) {
			for (j=0;j<=ct2->numofbases;++j) {
				write(&sav,&(allowed_alignments[i][j]));
			}
		}
	}

    short int kmax, lmax;
    for (i=0;i<=N;++i) {
      if (energyonly) {
        I = N;
      } else {
        I = i+N-1;
      }
      
      kmax = highend[i];//limit(i, maxsep, N, N2);
      for (j=i;j<=I;++j) {
        lmax = highend[j];//limit(j, maxsep, N, N2);
        for (k = lowend[i]/*limit(i, maxsep, N, N2)*/; k <= kmax; k++) {
          for (l = lowend[j]/*limit(j, maxsep, N, N2)*/; l <= lmax; l++) {
            if (j > N) {
              b = i;
              a = j-N;
            } else {
              b = j;
              a = i;
            }

            if (ct1->tem[b][a]) {
              write(&sav,&(v->f(i,j,k,l)));
            }

            write(&sav,&(w->f(i,j,k,l)));

            if (modification) {
              write(&sav,&(vmod->f(i,j,k,l)));
            }
          }
        }
      }
    }

    for (i=0;i<=N;++i) {
      for (j=lowend[i]/*0*/;j<=highend[i]/*2*maxsep+2*/;++j) {
        if (!energyonly) {
          write(&sav,&(w3->f(i,j))/*w3[i][j]*/);
        }
        write(&sav,&(w5->f(i,j))/*w5[i][j]*/);
      }
    }
    if (!energyonly) {
      for (j=0/*lowend[N+1]/*0*/;j<=N2+1/*highend[N+1]/*2*maxsep+2*/;++j) {
        write(&sav,&(w3->f(N+1,j))/*w3[N+1][j]*/);
      }
    }

    if (local) {
      flag=1;
    }
    else flag=0;
		
    write(&sav,&flag);

    //for (i=1;i<=2*ct1->numofbases;i++) {
    //	write(&sav,&pair[0][i]);			

    //}
    //for (i=1;i<=2*ct2->numofbases;i++) {
    //	write(&sav,&pair[1][i]);

    //}
	
    sav.close();
  }

  if (energyonly) {
	  //traceback the lowest free energy structure
     

	//set the basepr arrays to zero
	for (int index=1;index<=ct1->numofbases;index++) ct1->basepr[1][index]=0;
	for (int index=1;index<=ct2->numofbases;index++) ct2->basepr[1][index]=0;

	//set the alignment to zero
	for (int index=0;index<=ct1->numofbases;index++) alignment[0][index]=0;

	//one structure for each sequence
	ct1->numofstructures=1;
	ct2->numofstructures=1;

	//calculate the total energy 
	ct1->energy[1] = lowest;
	ct2->energy[1] = lowest;

	error = dyntrace(1, N, 1, N2, ct1, ct2, 1, alignment[0],
		  w, v, w3, w5, lowend, highend,data, gap, vmod, local, true);

  } else {
    error = dyntraceback(maxtracebacks,window,awindow,percentsort,
                 v,w,w3,w5,
                 ct1,ct2,
                 alignment,lowend,highend,gapincrease,data,
                 singleinsert,lowest,vmod,
                 local);
  }

  delete w;
  delete v;
  

  delete w3;
  delete w5;



  if (force) {
    


	
    for (i=0;i<2*N;++i) {
      delete[] fce1[i];	
    }

    delete[] fce1;

	

    for (i=0;i<2*N2;++i) {	
      delete[] fce2[i];
    }

    delete[] fce2;
    delete[] dbl1;
    delete[] dbl2;

    if (ct1->nmod>0||ct2->nmod>0) {
      //delete vmod
      delete vmod;
      

    }
  }

  delete[] lowend;
  delete[] highend;
  
#ifdef timer
  timeout << time(NULL)<<"\n";
  timeout << time(NULL) - seconds;
  timeout.close();
#endif	


#ifdef DYNALIGN_SMP
  //If progress was allocated in this function, remove it.
  if (removeprogress) delete progress;
#endif
  //return the errorcode, which was returned during traceback
  return error;


}



//calculate a single point in the v and w arrays - allowing constraints if force is true
void dynalignstep(structure *ct1, structure *ct2, datatable *data,
                  varray *v, dynalignarray *w, wendarray *w5, wendarray *w3, 
                  dynalignarray *vmod, bool *mod1, bool *mod2, bool modification,
                  char **fce1, char **fce2, bool alignmentforced, short **forcealign,
                  bool *dbl1, bool *dbl2, 
                  int i, int j, int k, int l, int N, int N2,
                  short *lowend, short *highend, int gap, bool singleinsert,
                  bool force,bool local,
                  short int (*edangle5)(int i, int j, int ip, structure *ct, datatable *data),
                  short int (*edangle3)(int i, int j, int ip, structure *ct, datatable *data)) {

  short int c,d,e,f;
  integersize idangle,jdangle,ldangle,kdangle,ijdangle,kldangle;
  bool ikincrement,jldecrement,ikincrement2,jldecrement2;
  short int startd,endd,starte,ende,startf,endf;

  short int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
                       {0,1,0,1,0,0},{0,0,0,0,0,0}};

  integersize einternal[maxloop+1][maxloop+1]; //precalculate internal loop free energies
  integersize einternal2[maxloop+1][maxloop+1]; //precalculate internal loop free energies

  integersize en1,en2,en3;



  //use a and b when testing whether a recursion is valid
  //a = k-i+maxsep;
  //if (l<=N2) b = l-j+maxsep;
  //else b = l-j+maxsep+Ndiff;

  //calculate some values for bounds checking
  /*kmin = lowlimit(i,maxsep,N,N2);
    kmax = highlimit(i,maxsep,N,N2);

    lmin = lowlimit(j,maxsep,N,N2);
    lmax = highlimit(j,maxsep,N,N2);*/


  //filter out isolated base pairs
  if (i+1<N&&j-1!=N&&i+1<j-1) e = inc[ct1->numseq[i+1]][ct1->numseq[j-1]];
  else e = 0;
  if (i-1>0&&j+1<2*N&&j!=N) e+=inc[ct1->numseq[i-1]][ct1->numseq[j+1]];


  if (k+1<ct2->numofbases&&l-1!=ct2->numofbases&&k+1<l-1) f = inc[ct2->numseq[k+1]][ct2->numseq[l-1]];
  else f = 0;
  if (k-1>0&&l+1<2*N2&&l!=N2) f += inc[ct2->numseq[k-1]][ct2->numseq[l+1]];


  //declare some variables needed for processing calculations with constraints
  short g;
  integersize en1m,en2m,en3m;
  bool helicalstack1,helicalstack2;

  if (force) {
    //if i,j,k, or l should be single stranded, don't allow the pair
    if (fce1[jref(i,j,N)][iref(i,j,N)]&SINGLE||fce2[jref(k,l,N2)][iref(k,l,N2)]&SINGLE||
        fce1[jref(i,j,N)][iref(i,j,N)]&NOPAIR||fce2[jref(k,l,N2)][iref(k,l,N2)]&NOPAIR) {
      e=0;
      f=0;
    }
    if (alignmentforced) {
      //if alignment is being forced, be sure that the pairs satisfy the requirement

      if (forcealign[0][i]>0&&forcealign[0][i]!=k) {
        e=0;
        f=0;

      }
      else if (forcealign[0][j]>0&&forcealign[0][j]!=l) {
        e=0;
        f=0;

      }
      else if (forcealign[1][k]>0&&forcealign[1][k]!=i) {
        e=0;
        f=0;

      }
      else if (forcealign[1][l]>0&&forcealign[1][l]!=j) {
        e=0;
        f=0;

      }
      else {
        //scan through to make sure that the base pair isn't precluding an alignment
        for (g=i+1;g<j;g++) {
          if (forcealign[0][g]>0) {

            if (forcealign[0][g]<k||forcealign[0][g]>l) {
              e=0;
              f=0;
            }
          }
        }

        for (g=k+1;g<l;g++) {
          if (forcealign[1][g]>0) {

            if (forcealign[1][g]<i||forcealign[1][g]>j) {
              e=0;
              f=0;
            }
          }
        }
      }
    }
  }

  // Use template information to decide whether to proceed with
  // calculating v(i,j,k,l)
  if (ct1->templated) {
    if (j<=N) {
      c=j;
      d=i;
    }
    else {
      d=j-N;
      c=i;
    }
    
    if (!ct1->tem[c][d]) {
      //forbind this pair
      e = 0;	
    }
  }
  if (ct2->templated) {
    if (l<=N2) {
      c=l;
      d=k;
    }
    else {
      d=l-N2;
      c=k;
    }

    if (!ct2->tem[c][d]) {
      //forbind this pair
      f = 0;	
    }
  }

  //test whether it is safe to increment or decrement the aligned nucleotides
						
  
  ikincrement = k+1 >= lowend[i+1] && k+1 <= highend[i+1];
  jldecrement = l-1 <= highend[j-1] && l-1 >= lowend[j-1];
												
  ikincrement2 = k+2 >= lowend[i+2] && k+2 <= highend[i+2];
  jldecrement2 = l-2 <= highend[j-2] && l-2 >= lowend[j-2];
						

  //Fill v, the best energy with i-j paired, k-l paired,
  //and the pair of basepairs aligned:
  if (inc[ct1->numseq[i]][ct1->numseq[j]]&&
      inc[ct2->numseq[k]][ct2->numseq[l]]&&e&&f)
    {
      //Now consider hairpins for both
      if (j<=N) v->f(i,j,k,l) =
                    erg3(i, j, ct1, data, force ? fce1[jref(i,j,N)][iref(i,j,N)] : 0) +
                    erg3(k, l, ct2, data, force ? fce2[jref(k,l,N2)][iref(k,l,N2)] : 0) +
                    gap*abs(j-i-l+k);
							
        en1 = DYNALIGN_INFINITY;

		if (j > N ) {//&& jldecrement && (ikincrement):This needs to be more tailored for each case
          
			//Consider the exterior loop closed by i-j (actually j-N paired to i):
          //must consider whether an unpaired nucleotide is stacked
          //is stacked onto each of the nucleotides in a helix
          //There are 16 cases:
          //consider them in this order (0 unstacked, 1 stacked)
          //		i	j	k	l
          //1		0	0	0	0
          //2		0	0	0	1
          //3		0	0	1	0
          //4		0	0	1	1
          //5		0	1	0	0
          //6		0	1	0	1
          //7		0	1	1	0
          //8		0	1	1	1
          //9		1	0	0	0
          //10	1	0	0	1
          //11	1	0	1	0
          //12	1	0	1	1
          //13	1	1	0	0
          //14	1	1	0	1
          //15	1	1	1	0
          //16	1	1	1	1

          //to save time, only do the function calls once:
          idangle = edangle3(i,j,i+1,ct1,data);
          jdangle = edangle5(j,i,j-1,ct1,data);
          if (k<N2-1) kdangle = edangle3(k,l,k+1,ct2,data);
          else kdangle = DYNALIGN_INFINITY;
          if (l>N2+1) ldangle = edangle5(l,k,l-1,ct2,data);
          else ldangle=DYNALIGN_INFINITY;

								
		  //allow an exception for jldecrement2 if this reaches the end of the sequence
		  if (j-2==N) jldecrement2=true;


          //case 1, no stacks
		  if (jldecrement) {
			  if (ikincrement) en1 = min(en1,w5->f(j-1-N,l-1-N2)+w3->f(i+1,k+1));
			  else if (k==N2&&local) en1 = min(en1,w5->f(j-1-N,l-1-N2));
			  else if (k==N2) en1 = min(en1,w5->f(j-1-N,l-1-N2)+gap*(N-i));

		  }

          //case 16, all stack
          if ((j-2-N>=0)&&(l-2-N2>0)) if ((i<N)&&(j>N+1)&&(k<N2)&&jldecrement2&&ikincrement2) en1 = min(en1,w5->f(j-2-N,l-2-N2)+w3->f(i+2,k+2)/*w5[j-2-N][b]+w3[i+2][a]*/+jdangle+ldangle
                                                                    +idangle+kdangle);

          //case 6, j and l stacked
		  if (j-2-N>=0&&l-2-N2>0) if (j>N+1&&jldecrement2) {
			if (ikincrement) en1 = min(en1,w5->f(j-2-N,l-2-N2)+w3->f(i+1,k+1)/*w5[j-2-N][b]+w3[i+1][a]*/+jdangle+ldangle);
			else if (k==N2&&local)  en1 = min(en1,w5->f(j-2-N,l-2-N2)+jdangle+ldangle); 
			else if (k==N2) en1 = min(en1,w5->f(j-2-N,l-2-N2)+jdangle+ldangle+gap*(N-i)); 

		  }

          //case 11, i and k stacked
          if ((i<N)&&ikincrement2&&(k<N2)&&jldecrement) en1 = min(en1,w5->f(j-1-N,l-1-N2)+w3->f(i+2,k+2)/*w5[j-1-N][b]+w3[i+2][a]*/+idangle+kdangle);

          if (l-2-N2>0)if (/*b>0*/l-2-N2>=lowend[j-1-N]&&l-2-N2<=highend[j-1-N]) {
            //case 2, l stacked
            if (ikincrement) en1 = min(en1,w5->f(j-1-N,l-2-N2)+w3->f(i+1,k+1)/*w5[j-1-N][b-1]+w3[i+1][a]*/+ldangle+gap);
			else if (k==N2&&local) en1 = min(en1,w5->f(j-1-N,l-2-N2)+ldangle+gap);
			else if (k==N2) en1 = min(en1,w5->f(j-1-N,l-2-N2)+ldangle+gap+gap*(N-i));

            //case 12, i, k, and l stacked
            if ((i<N)&&ikincrement2&&(k<N2)) en1 = min(en1,w5->f(j-1-N,l-2-N2)+w3->f(i+2,k+2)/*w5[j-1-N][b-1]+w3[i+2][a]*/+ldangle+idangle
                                               +kdangle+gap);

            if (k+1<highend[i+1]&&k+1>=lowend[i+1]&&(k<N2)) {
              //case 4, l and k stacked
              en1 = min(en1,w5->f(j-1-N,l-2-N2)+w3->f(i+1,k+2)/*w5[j-1-N][b-1]+w3[i+1][a+1]*/+ldangle+kdangle+2*gap);
            }

            if (k+1>=lowend[i+2]&&(i<N)) {
              //case 10, i and l stacked
              if (k+1<=highend[i+2]) en1 = min(en1,w5->f(j-1-N,l-2-N2)+w3->f(i+2,k+1)/*w5[j-1-N][b-1]+w3[i+2][a-1]*/+ldangle+idangle+2*gap);
			  else if (k==N2 && local) en1 = min(en1,w5->f(j-1-N,l-2-N2)+ldangle+idangle+2*gap);
			  else if (k==N2) en1 = min(en1,w5->f(j-1-N,l-2-N2)+ldangle+idangle+2*gap+gap*(N-i-1));
            }

									

          }

          if (k+2<=highend[i+1]&&k+2>=lowend[i+1]&&(k<N2)) {
            //case 3, k alone stacked
            if (jldecrement) en1 = min(en1,w5->f(j-1-N,l-1-N2)+w3->f(i+1,k+2)/*w5[j-1-N][b]+w3[i+1][a+1]*/+kdangle+gap);

			if (j-2-N>=0) {
            //case 8, j, k, and l stacked
				if (l-2-N2>0) if ((j>N+1)&&jldecrement2) en1 = min(en1,w5->f(j-2-N,l-2-N2)+w3->f(i+1,k+2)/*w5[j-2-N][b]+w3[i+1][a+1]*/+kdangle
													 +jdangle+ldangle+gap);
										
				if (l-1-N2<=highend[j-2-N]&&l-1-N2>=lowend[j-2-N]/*b<2*maxsep+2*/) {
				  //case 7, j and k stacked
				  if ((j>N+1)) en1 = min(en1,w5->f(j-2-N,l-1-N2)+w3->f(i+1,k+2)/*w5[j-2-N][b+1]+w3[i+1][a+1]*/+kdangle+jdangle+2*gap);
				}
			}
          }

		  if (j-2-N>=0) {
			  if (l-1-N2<=highend[j-2-N]&&l-1-N2>=lowend[j-2-N]) {
				//case 5, j stacked
				  if ((j>N+1)) {
					  
					  if (ikincrement) en1 = min(en1,w5->f(j-2-N,l-1-N2)+w3->f(i+1,k+1)/*w5[j-2-N][b+1]+w3[i+1][a]*/+jdangle+gap);
						else if (k==N2&&local) en1 = min(en1,w5->f(j-2-N,l-1-N2)+jdangle+gap);
						else if (k==N2) en1 = min(en1,w5->f(j-2-N,l-1-N2)+jdangle+gap+gap*(N-i));
				  }
				//case 15, i, j, and k stacked
				if ((i<N)&&(j>N+1)&&ikincrement2&&k+2<=N2) en1 = min(en1,w5->f(j-2-N,l-1-N2)+w3->f(i+2,k+2)/*w5[j-2-N][b+1]+w3[i+2][a]*/+jdangle
															+idangle+kdangle+gap);
										
				if (k+1>=lowend[i+2]&&(i<N)&&(j>N+1)) {
				  //case 13, j and i stacked
					if (k+1<=highend[i+2]) {
						en1 = min(en1,w5->f(j-2-N,l-1-N2)+w3->f(i+2,k+1)/*w5[j-2-N][b+1]+w3[i+2][a-1]*/+jdangle+idangle+2*gap);
					}
					else if (k==N2&&local) en1 = min(en1,w5->f(j-2-N,l-1-N2)+jdangle+idangle+2*gap); 
					else if (k==N2)en1 = min(en1,w5->f(j-2-N,l-1-N2)+jdangle+idangle+2*gap+gap*(N-i-1)); 
				}
			  }
		  }

          if (k+1>=lowend[i+2]) {
            //case 9, i stacked
			  if ((i<N)&&jldecrement) {
				  if (k+1<=highend[i+2]) en1 = min(en1,w5->f(j-1-N,l-1-N2)+w3->f(i+2,k+1)/*w5[j-1-N][b]+w3[i+2][a-1]*/+idangle+gap);
				  else if (k==N2&&local) en1 = min(en1,w5->f(j-1-N,l-1-N2)+idangle+gap);
				  else if (k==N2) en1 = min(en1,w5->f(j-1-N,l-1-N2)+idangle+gap+gap*(N-i-1));
			  }
            //case 14, i, j, and l stacked
			  if ((i<N)&&(j>N+1)&&jldecrement2&&(j-2-N>=0)&&(l-2-N2>0) ) {
				  if (k+1<=highend[i+2]) en1 = min(en1,w5->f(j-2-N,l-2-N2)+w3->f(i+2,k+1)/*w5[j-2-N][b]+w3[i+2][a-1]*/+idangle
                                                        +jdangle+ldangle+gap);
				  else if (k==N2&&local) en1 = min(en1,w5->f(j-2-N,l-2-N2)+idangle
                                                        +jdangle+ldangle+gap);
				  else if (k==N2) en1 = min(en1,w5->f(j-2-N,l-2-N2)+idangle
                                                        +jdangle+ldangle+gap+gap*(N-i-1));

			  }
          }


		  //reset jl decrement2 for remainder of code... i.e. undo the exception above
		  if (j-2==N) jldecrement2 = l-2 <= highend[j-2] && l-2 >= lowend[j-2];
								
		}

		

        if (i!=N&&j!=N+1&&k!=N2&&l!=N2+1&&ikincrement&&jldecrement) {
          //Now consider multiloop:
          //junction closed by i-j pair aligned with a and b
          //calculate the free energy of 2 fragments merged:

          idangle = edangle3(i,j,i+1,ct1,data);
          jdangle = edangle5(j,i,j-1,ct1,data);
          kdangle = edangle3(k,l,k+1,ct2,data);
          ldangle = edangle5(l,k,l-1,ct2,data);


          for (c=i+minloop+1;c<j-minloop;++c) {
            //if (c>N) kp = Ndiff;
            //else kp=0;
            startd=max(k+minloop+1,lowend[c]);
            endd=min(highend[c],l-minloop);
            for (d=startd/*c-maxsep-kp*/;
                 d<endd/*d<l-minloop&&d<=c+maxsep-kp*/;++d) {
											
              //e = d-c+maxsep+kp;
											
              if (((c<N)&&(d<N2)||((c>N)&&(d>N2)))&&d+1>=lowend[c+1]&&d+1<=highend[c+1]) {
												
									
                //must consider whether an unpaired nucleotide 
                //is stacked onto each of the nucleotides in a helix
                //There are 16 cases:
                //consider rthem in this order (0 unstacked, 1 stacked)
                //		i	j	k	l
                //1		0	0	0	0
                //2		0	0	0	1
                //3		0	0	1	0
                //4		0	0	1	1
                //5		0	1	0	0
                //6		0	1	0	1
                //7		0	1	1	0
                //8		0	1	1	1
                //9		1	0	0	0
                //10	1	0	0	1
                //11	1	0	1	0
                //12	1	0	1	1
                //13	1	1	0	0
                //14	1	1	0	1
                //15	1	1	1	0
                //16	1	1	1	1

												

                //case 1 - no stacks
                en1 = min(en1,w->f(i+1,c,k+1,d)+w->f(c+1,j-1,d+1,l-1)/*w[c][iref(i+1,c,N)][a][e]*/
                          /*w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]);
								

                //case 16 - all four stacked
												
                if (((i+1)!=N)&&((j-1)!=N)&&((k+1)!=N2)&&((l-1)!=N2)&&ikincrement2&&jldecrement2) {
                  en1 = min(en1,w->f(i+2,c,k+2,d)+w->f(c+1,j-2,d+1,l-2)/*w[c][iref(i+2,c,N)][a][e]
                                                                         +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                            +4*data->eparam[6]+ldangle
                            +kdangle
                            +jdangle
                            +idangle);
                }

								
                //case 6 - j and l stacked:
                if (((j-1)!=N)&&((l-1)!=N2)&&jldecrement2) {
                  en1 = min(en1,w->f(i+1,c,k+1,d)+w->f(c+1,j-2,d+1,l-2)/*w[c][iref(i+1,c,N)][a][e]
                                                                         +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                            +2*data->eparam[6]+jdangle
                            +ldangle);
                }

									
                //case 11 - i and k stacked
                if (((i+1)!=N)&&((k+1)!=N2)&&ikincrement2) {
                  en1 = min(en1,w->f(i+2,c,k+2,d)+w->f(c+1,j-1,d+1,l-1)/*w[c][iref(i+2,c,N)][a][e]
                                                                         +w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                            +2*data->eparam[6]+kdangle+
                            idangle);
                }



                if (l-2>=lowend[j-1]&&l-2<=highend[j-1]&&((l-1)!=N2)) {
                  //case 2 - stack on l
                  en1 = min(en1,w->f(i+1,c,k+1,d)+w->f(c+1,j-1,d+1,l-2)/*w[c][iref(i+1,c,N)][a][e]
                                                                         +w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                            +data->eparam[6]+ldangle+gap);
							
                  if ((i+1)!=N) {
                    //case 12 - i, k, and l stacked
                    if ((k+1)!=N2&&ikincrement2) {
                      en1 = min(en1,w->f(i+2,c,k+2,d)+w->f(c+1,j-1,d+1,l-2)/*w[c][iref(i+2,c,N)][a][e]
                                                                             +w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                                +3*data->eparam[6]+ldangle
                                +idangle+kdangle+gap);
                    }


                    if (k+1>=lowend[i+2]&&k+1<=highend[i+2]) {
                      //case 10 - l and i stacked
                      en1 = min(en1,w->f(i+2,c,k+1,d)+w->f(c+1,j-1,d+1,l-2)/*w[c][iref(i+2,c,N)][a-1][e]
                                                                             +w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                                +2*data->eparam[6]+ldangle
                                +idangle+2*gap);

									
                    }
                  }
                  if (k+2<=highend[i+1]&&k+2>=lowend[i+1]&&((k+1)!=N2)) {
                    //case 4 - k and l stacked
                    en1 = min(en1,w->f(i+1,c,k+2,d)+w->f(c+1,j-1,d+1,l-2)/*w[c][iref(i+1,c,N)][a+1][e]
                                                                           +w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                              +2*data->eparam[6]+ldangle
                              +kdangle+2*gap);
                  }
                }
                if (k+1>=lowend[i+2]&&k+1<=highend[i+2]&&(i+1!=N)) {

                  //case 9 - i stacked
                  en1 = min(en1,w->f(i+2,c,k+1,d)+w->f(c+1,j-1,d+1,l-1)/*w[c][iref(i+2,c,N)][a-1][e]
                                                                         +w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                            +data->eparam[6]+idangle+gap);

                  if (j-1!=N) {
                    //case 14 - i, j, and l stacked
                    if (l-1!=N2&&jldecrement2) {
                      en1 = min(en1,w->f(i+2,c,k+1,d)+w->f(c+1,j-2,d+1,l-2)/*w[c][iref(i+2,c,N)][a-1][e]
                                                                             +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                                +3*data->eparam[6]+ldangle
                                +jdangle
                                +idangle+gap);
                    }
			      
                    if (l-1>=lowend[j-2]&&l-1<=highend[j-2]/*b+1<2*maxsep+2*/) {
                      //case 13 - i and j stacked
                      en1 = min(en1,w->f(i+2,c,k+1,d)+w->f(c+1,j-2,d+1,l-1)/*w[c][iref(i+2,c,N)][a-1][e]
                                                                             +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                                +2*data->eparam[6]+jdangle
                                +idangle+2*gap);
                    }
                  }
                }
									
                if (k+2<=highend[i+1]&&k+2>=lowend[i+1]/*(a+1<2*maxsep+2)*/&&(k+1!=N2)) {
                  //case 3 - stack on k 
                  en1 = min(en1,w->f(i+1,c,k+2,d)+w->f(c+1,j-1,d+1,l-1)/*w[c][iref(i+1,c,N)][a+1][e]*/+2*data->eparam[10]
                            +/*w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b]*/+2*data->eparam[5]
                            +data->eparam[6]+kdangle+gap);

                  if (j-1!=N) {
                    //case 8 - j, k, and l stacked
                    if (l-1!=N2&&jldecrement2) {
                      en1 = min(en1,w->f(i+1,c,k+2,d)+w->f(c+1,j-2,d+1,l-2)/*w[c][iref(i+1,c,N)][a+1][e]
                                                                             +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                                +3*data->eparam[6]+jdangle
                                +ldangle+kdangle+gap);
                    }			
                    if (l-1<=highend[j-2]&&l-1>=lowend[j-2]) {
                      //case 7 - j and k stacked
                      en1 = min(en1,w->f(i+1,c,k+2,d)+w->f(c+1,j-2,d+1,l-1)/*w[c][iref(i+1,c,N)][a+1][e]
                                                                             +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                                +2*data->eparam[6]+jdangle
                                +kdangle+2*gap);																		
                    }	
                  }
										

                }
                if (l-1<=highend[j-2]&&l-1>=lowend[j-2]&&(j-1!=N)) {
                  //case 15 - i, j, and k stacked
                  if ((i+1!=N)&&(k+1!=N2)&&ikincrement2) {
                    en1 = min(en1,w->f(i+2,c,k+2,d)+w->f(c+1,j-2,d+1,l-1)/*w[c][iref(i+2,c,N)][a][e]
                                                                           +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                              +3*data->eparam[6]+kdangle
                              +jdangle
                              +idangle+gap);
                  }
						
                  //case 5 - j stacked
                  en1 = min(en1,w->f(i+1,c,k+1,d)+w->f(c+1,j-2,d+1,l-1)/*w[c][iref(i+1,c,N)][a][e]
                                                                         +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                            +data->eparam[6]+jdangle+gap);
                }
              }
            }
          }
        }
        en1 = en1 + penalty(i,j,ct1,data) + penalty(k,l,ct2,data);
        v->f(i,j,k,l)=min(en1,v->f(i,j,k,l));


        en3 = DYNALIGN_INFINITY;
        en3m = DYNALIGN_INFINITY; // Only used if force is on, but it's
        // simpler just to do the assign than to
        // check force.


        //first precalculate some internal loop free energies:
        for (c=i+1;c<=i+maxloop&&c<j/*-minloop*/&&c<=N;++c) {
          for (d=j-1;d>=j-maxloop&&d>c/*+minloop*/;--d) {
            if (c-i+j-d-2<maxinternal) {
              einternal[c-i][j-d] = erg2(i, j, c, d, ct1, data,
                                         force ? fce1[jref(i,c,N)][iref(i,c,N)] : 0,
                                         force ? fce1[jref(d,j,N)][iref(d,j,N)] : 0);
            }
          }
        }

        for (c=k+1;c<=k+maxloop&&c<l/*-minloop*/&&c<=N2;++c) {
          for (d=l-1;d>=l-maxloop&&d>c/*+minloop*/;--d) {
            if (c-k+l-d-2<maxinternal) {
              einternal2[c-k][l-d] = erg2(k, l, c, d, ct2, data,
                                          force ? fce2[jref(k,c,N2)][iref(k,c,N2)] : 0,
                                          force ? fce2[jref(d,l,N2)][iref(d,l,N2)] : 0);
            }
          }
        }

        //Now consider internal loops/one stacked and other not
        for (c=i+1;c<=i+maxloop&&c<j/*-minloop*/&&c<=N;++c) {
          for (d=j-1;d>=j-maxloop&&d>c/*+minloop*/;--d) {
            //if (d>N) ap = Ndiff;
            //else ap=0;
            starte=max(k+1,lowend[c]);
            ende=min(min(min(k+maxloop,l),highend[c]),N2);
            for (e=starte/*max(k+1,c-maxsep)*/;e<=ende/*k+maxloop&&e<l&&(e-c)<=maxsep&&e<=N2*/;++e) {
              endf=min(l-1,highend[d]);
              startf=max(max((l-maxloop),e+1),lowend[d]);
              for (f=endf/*min(l-1,d+maxsep-ap)*/;f>=startf/*l-maxloop&&f>e&&f>=d-maxsep-ap*/;--f) {

                if (c-i+j-d-2<maxinternal&&e-k+l-f-2<maxinternal&&(((d>N)&&(f>N2))||((d<=N)&&(f<=N2)))) {	

                  {
                    // Only used if force is on, but it's simpler just to do
                    // the assign than to check force.
                    helicalstack1=false; 
                    helicalstack2=false;
                  }
              
                  if (c==i+1&&d==j-1&&(i!=N)&&(j-1!=N)) {
                    //seq1 is helical stack
                    en1 = erg1(i,j,i+1,j-1,ct1, data);
                
                    if (force && modification) {
                      //calculate en1m
                      if (mod1[i]||mod1[j]) {
                        //check for a GU exception
                        if ((ct1->numseq[i]==3&&ct1->numseq[j]==4)||
                            (ct1->numseq[i]==4&&ct1->numseq[j]==3)||
                            (ct1->numseq[i+1]==3&&ct1->numseq[j-1]==4)||
                            (ct1->numseq[i+1]==4&&ct1->numseq[j-1]==3)||
                            (ct1->numseq[i-1]==3&&ct1->numseq[j+1]==4)||
                            (ct1->numseq[i-1]==4&&ct1->numseq[j+1]==3))
                          en1m=en1;
                        else en1m=DYNALIGN_INFINITY;
                      }
                      else en1m=en1;
                      helicalstack1=true;
                    }
                  }
                  else {
                    //internal loop
                    en1 = einternal[c-i][j-d];
                    //erg2(i,j,c,d,ct1,data,fce1[jref(i,c,N)][iref(i,c,N)],fce1[jref(d,j,N)][iref(d,j,N)]);
                    //erg2(i,j,c,d,ct1,data,0,0);

                    if (force && modification) {
                      en1m=en1;
                    }
                  }

                  if (e==k+1&&f==l-1&&k!=N2&&(l-1!=N2)) {
                    //seq2 is helical stack
                    en2 = erg1(k,l,k+1,l-1,ct2, data);
                    if (force && modification) {
                      if (mod2[k]||mod2[l]) {
                        //check for a GU exception
                        if ((ct2->numseq[k]==3&&ct2->numseq[l]==4)||
                            (ct2->numseq[k]==4&&ct2->numseq[l]==3)||
                            (ct2->numseq[k+1]==3&&ct2->numseq[l-1]==4)||
                            (ct2->numseq[k+1]==4&&ct2->numseq[l-1]==3)||
                            (ct2->numseq[k-1]==3&&ct2->numseq[l+1]==4)||
                            (ct2->numseq[k-1]==4&&ct2->numseq[l+1]==3))
                          en2m = en2;
                        else en2m = DYNALIGN_INFINITY;
                      }
                      else en2m = en2;
                      helicalstack2=true;
                    }
												
                    // also allow single base pair insertions into one sequence only
                    if (singleinsert&&c==i+2&&d==j-2&&j-2>i+2&&
                        inc[ct1->numseq[i+1]][ct1->numseq[j-1]]&&
                        inc[ct1->numseq[i+2]][ct1->numseq[j-2]]&&i+1!=N&&
                        (j-1!=N)&&(j-2!=N)

                        && (!force ||
                            (!(fce1[jref(i+1,j-1,N)][iref(i+1,j-1,N)] & SINGLE) &&
                             !(fce1[jref(i+1,j-1,N)][iref(i+1,j-1,N)] & NOPAIR) &&
                             !mod1[i+1] &&
                             !mod1[j-1]
                             )
                            )
                        ) {

                      en1 = min(en1,erg1(i,j,i+1,j-1,ct1, data)+
                                erg1(i+1,j-1,i+2,j-2,ct1,data));
                    }
                  } else {
                    //internal loop

                    en2 = einternal2[e-k][l-f];
                    //erg2(k,l,e,f,ct2,data,fce2[jref(k,e,N2)][iref(k,e,N2)],fce2[jref(f,l,N2)][iref(f,l,N2)]);
                    //erg2(k,l,e,f,ct2,data,0,0);
                    if (force && modification) {
                      en2m = en2;
                    }
											
                    //also allow single base pair insertions into one sequence only
                    if (singleinsert&&e==k+2&&f==l-2&&l-2>k+2&&c==i+1&&d==j-1&&
                        inc[ct2->numseq[k+1]][ct2->numseq[l-1]]&&
                        inc[ct2->numseq[k+2]][ct2->numseq[l-2]]&&
                        inc[ct1->numseq[i+1]][ct1->numseq[j-1]]&&
                        i!=N&&(j-1!=N)&&k!=N2&&(l-1!=N2)&&
                        (k+1!=N2)&&(l-2!=N2)
                    
                        && (!force ||
                            (!(fce2[jref(k+1,l-1,N2)][iref(k+1,l-1,N2)] & SINGLE) &&
                             !(fce2[jref(k+1,l-1,N2)][iref(k+1,l-1,N2)] & NOPAIR) &&
                             !mod2[k+1] &&
                             !mod2[l-1]
                             )
                            )
                        ) {
                  
                      en2 = min(en2,erg1(k,l,k+1,l-1,ct2, data)+
                                erg1(k+1,l-1,k+2,l-2,ct2,data));
                    }
                  }
              
                  if (force) {
                    if (modification) {
                      if (helicalstack1||helicalstack2) {
                        en3 = min(en3,en1+en2+vmod->f(c,d,e,f)/*vmod[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                        en3m = min(en3m,en1m+en2m+vmod->f(c,d,e,f)/*vmod[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                      }
                      else {
                        en3 = min(en3,en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                        en3m = min(en3m,en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                      }
                  
                    } else {
                      en3 = min(en3,en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                    }
                  } else {
                    en3 = min(en3,en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                  }
                }
              }
            }
          }
        }

        if (force && modification) {
          /*if (j<=N)*/ vmod->f(i,j,k,l)/*vmod[j][i][a][b]*/=min(en3m,v->f(i,j,k,l)/*v[j][i][a][b]*/);
          //else vmod[j][i-j+N][a][b]=min(en3m,v[j][i-j+N][a][b]);
        }
    
        /*if (j<=N)*/ v->f(i,j,k,l)/*v[j][i][a][b]*/=min(en3,v->f(i,j,k,l)/*v[j][i][a][b]*/);
        //else v[j][i-j+N][a][b]=min(en3,v[j][i-j+N][a][b]);
      }
    else {
      v->f(i,j,k,l)/*v[j][iref(i,j,N)][a][b]*/ = DYNALIGN_INFINITY;
      //a bp wasn't allowed between either i and j or k and l
    }


						
  /////////////////////////////////////////////////////////
  //Fill w, the best energy for a fragment in a multiloop:

  //save a temporaray value for w in en1, a register short int
						

  //consider the possibilities for adding nucleotides to existing fragments
  //		i	j	k	l
						
  //2		0	0	0	1
  //3		0	0	1	0
  //4		0	0	1	1z
  //5		0	1	0	0
  //6		0	1	0	1z
  //7		0	1	1	0
  //8		0	1	1	1
  //9		1	0	0	0
  //10	1	0	0	1
  //11	1	0	1	0z
  //12	1	0	1	1
  //13	1	1	0	0z
  //14	1	1	0	1
  //15	1	1	1	0
  //16	1	1	1	1z

  //remember than j<=N is included fragment and j>N is excluded fragment
  en1=DYNALIGN_INFINITY;

  //case 6 0	1	0	1
  if (j-1 != N && l-1 != N2 && jldecrement &&
      (!force ||
       (!dbl1[j] && !dbl2[l])
       )
      ) {
    en1 = w->f(i,j-1,k,l-1)/*w[j-1][iref(i,j-1,N)][a][b]*/+2*data->eparam[6];

    //case 16 1	1	1	1
    if (j-1 > i+1 && i != N && k != N2 &&
        ikincrement && jldecrement &&
        (!force ||
         (!dbl1[i] && !dbl2[k])
         )
        ) {
      en1 = min(en1,w->f(i+1,j-1,k+1,l-1)/*w[j-1][iref(i+1,j-1,N)][a][b]*/+ 4*data->eparam[6]);
    }
  }

  //case 11	1	0	1	0
  if (i != N && k != N2 && ikincrement &&
      (!force ||
       (!dbl1[i] && !dbl2[k])
       )
      ) {
    en1 = min(en1,w->f(i+1,j,k+1,l)/*w[j][iref(i+1,j,N)][a][b]*/+2*data->eparam[6]);
  }
  if (k>=lowend[i+1]&&k<=highend[i+1]/*a>=1*/) {
    //case 9		1	0	0	0
    if (i != N &&
        (!force ||
         (!dbl1[i])
         )
        ) {
      en1 = min(en1,w->f(i+1,j,k,l)/*w[j][iref(i+1,j,N)][a-1][b]*/+data->eparam[6]+gap);

      //case 14	1	1	0	1
      if ((j-1>i+1)&&(j-1!=N)) {
        if (l-1 != N2 && jldecrement &&
            (!force ||
             (!dbl1[j]&&!dbl2[l])
             )
            ) {
          en1 = min(en1,w->f(i+1,j-1,k,l-1)/*w[j-1][iref(i+1,j-1,N)][a-1][b]*/+3*data->eparam[6]+gap);
        }
        if (l <= highend[j-1]&&l>=lowend[j-1]/*b+1<(2*maxseparation+2)*/ &&
            (!force ||
             (!dbl1[j])
             )
            ) {

          //case 13 1	1	0	0
          en1 = min(en1,w->f(i+1,j-1,k,l)/*w[j-1][iref(i+1,j-1,N)][a-1][b+1]*/+2*data->eparam[6]+2*gap);
        }
      }
    }
      
    if (l-1 >= lowend[j]/*b>=1*/ && i != N && l-1 != N2 &&
        (!force ||
         (!dbl1[i] && !dbl2[l])
         )
        ) {
      //case 10 1	0	0	1
      en1 = min(en1,w->f(i+1,j,k,l-1)/*w[j][iref(i+1,j,N)][a-1][b-1]*/+2*data->eparam[6]+2*gap);
    }
  }
  if (j-1 != N &&
      (l <= highend[j-1] && l>=lowend[j-1])) {
    //case 5		0	1	0	0
    if (!force ||
        (!dbl1[j])
        ) {
      en1 = min(en1,w->f(i,j-1,k,l)/*w[j-1][iref(i,j-1,N)][a][b+1]*/+data->eparam[6]+gap);
    }
    //case 15	1	1	1	0
    if (k!=N2) {
      if (j-1 > i+1 && i != N && ikincrement &&
          (!force ||
           (!dbl1[i] && !dbl1[j] && !dbl2[k])
           )
          ) {
        en1 = min(en1,w->f(i+1,j-1,k+1,l)/*w[j-1][iref(i+1,j-1,N)][a][b+1]*/+3*data->eparam[6]+gap);
      }
      if (k+1 <= highend[i]/*a+1<(2*maxseparation+2)*/ &&
          (!force ||
           (!dbl1[j] && !dbl2[k])
           )
          ) {
        //case 7		0	1	1	0
        en1=min(en1,w->f(i,j-1,k+1,l)/*w[j-1][iref(i,j-1,N)][a+1][b+1]*/+2*data->eparam[6]+2*gap);

      }
    }
  }
  if (k+1 <= highend[i]/*a+1<(2*maxseparation+2)*/ && k != N2 &&
      (!force ||
       (!dbl2[k])
       )
      ) {
    //case 3		0	0	1	0
    en1 = min(en1,w->f(i,j,k+1,l)/*w[j][iref(i,j,N)][a+1][b]*/+data->eparam[6]+gap);

    if (l-1 != N2 &&
        (!force ||
         (!dbl2[l])
         )
        ) {
      //case 8		0	1	1	1
      if (j-1 != N && jldecrement &&
          (!force ||
           (!dbl1[j])
           )
          ) {
        en1 = min(en1,w->f(i,j-1,k+1,l-1)/*w[j-1][iref(i,j-1,N)][a+1][b]*/+3*data->eparam[6]+gap);
      }					
      if (l-1 >= lowend[j]/*b>=1*/) {
        //case 4		0	0	1	1
        en1 = min(en1,w->f(i,j,k+1,l-1)/*w[j][iref(i,j,N)][a+1][b-1]*/+2*data->eparam[6]+2*gap);
      }
    }
  }
  if (l-1 >= lowend[j]/*(b>=1)*/ && l-1 != N2 &&
      (!force ||
       (!dbl2[l])
       )
      ) {
    //case 2		0	0	0	1
    en1=min(en1,w->f(i,j,k,l-1)/*w[j][iref(i,j,N)][a][b-1]*/+data->eparam[6]+gap);

    //case12
    //12	1	0	1	1
    if (i != N && k != N2 && ikincrement &&
        (!force ||
         (!dbl1[i] && !dbl2[k])
         )
        ) {
      en1 = min(en1,w->f(i+1,j,k+1,l-1)/*w[j][iref(i+1,j,N)][a][b-1]*/+3*data->eparam[6]+gap);
    }
  }
  
  //Consider the case where none of the four nucs (i,j,k,l) are paired
  //Consider whether any of the 4 nucleotides are stacked on a helix
  //There are 16 cases:
  //consider them in this order (0 unstacked, 1 stacked)
  //		i	j	k	l
  //1		0	0	0	0
  //2		0	0	0	1
  //3		0	0	1	0
  //4		0	0	1	1
  //5		0	1	0	0
  //6		0	1	0	1
  //7		0	1	1	0
  //8		0	1	1	1
  //9		1	0	0	0
  //10	1	0	0	1
  //11	1	0	1	0
  //12	1	0	1	1
  //13	1	1	0	0
  //14	1	1	0	1
  //15	1	1	1	0
  //16	1	1	1	1


  //note:
  //	a = k-i+maxsep;
  //	b = l-j+maxsep;

  //	so when addressing i+1 => address a-1 to keep k unchanged 
  //	and simlarly: j-1 => b+1 to keep l unchanged

  idangle=edangle5(i+1,j,i,ct1,data);
  jdangle=edangle3(j-1,i,j,ct1,data);
  kdangle=edangle5(k+1,l,k,ct2,data);
  ldangle=edangle3(l-1,k,l,ct2,data);
  ijdangle=edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data);
  kldangle=edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data);

  //case 1 - nothing stacked:
  en1=min(en1,v->f(i,j,k,l)/*v[j][iref(i,j,N)][a][b]*/+2*data->eparam[10]+penalty(i,j,ct1,data)+penalty(k,l,ct2,data));
						
  //case 6 - j and l stacked
  if ((j-1!=N)&&(l!=N2+1)&&jldecrement) {
    en1 = min(en1,v->f(i,j-1,k,l-1)/*v[j-1][iref(i,j-1,N)][a][b]*/+2*data->eparam[6]+jdangle
              +ldangle+2*data->eparam[10]+penalty(i,j-1,ct1,data)+penalty(k,l-1,ct2,data));
  }
  //case 11 - i and k stacked
  if ((i!=N)&&(k!=N2)&&ikincrement) {
    en1 = min(en1,v->f(i+1,j,k+1,l)/*v[j][iref(i+1,j,N)][a][b]*/+2*data->eparam[6]+kdangle+
              idangle+2*data->eparam[10]+penalty(i+1,j,ct1,data)+penalty(k+1,l,ct2,data));

    //case 16 - i, j, k, and l stacked
    if (j-1>i+1&&(j-1!=N)&&(l!=N2+1)&&ikincrement&&jldecrement) {
      en1=min(en1,v->f(i+1,j-1,k+1,l-1)/*v[j-1][iref(i+1,j-1,N)][a][b]*/+4*data->eparam[6]
              +kldangle
              +ijdangle
              +2*data->eparam[10]
              +penalty(i+1,j-1,ct1,data)+penalty(k+1,l-1,ct2,data));
    }
  }

  if (l-1>=lowend[j]/*(b-1>=0)*/&&(l!=N2+1)) {
    //case 2 - l stacked
    en1=min(en1,v->f(i,j,k,l-1)/*v[j][iref(i,j,N)][a][b-1]*/+data->eparam[6]+ldangle+2*data->eparam[10]+gap+
            penalty(i,j,ct1,data)+penalty(k,l-1,ct2,data));
						
    if (k+1<=highend[i]/*(a+1<2*maxsep+2)*/&&(k!=N2)) {
      //case 4 - l and k stacked
      en1=min(en1,v->f(i,j,k+1,l-1)/*v[j][iref(i,j,N)][a+1][b-1]*/+2*data->eparam[6]
              +kldangle+2*data->eparam[10]+2*gap+
              penalty(i,j,ct1,data)+penalty(k+1,l-1,ct2,data));
    }
  }

  if (k+1<=highend[i]/*(a+1<2*maxsep+2)*/&&(k!=N2)) {

							

    //case 3 - k stacked
    en1 = min(en1,v->f(i,j,k+1,l)/*v[j][iref(i,j,N)][a+1][b]*/+data->eparam[6]+kdangle+2*data->eparam[10]+gap
              +penalty(i,j,ct1,data)+penalty(k+1,l,ct2,data));

    if (j!=N+1) {
      if (l<=highend[j-1]&&l>=lowend[j-1]) {
        //case 7 - j and k stacked
        en1 = min(en1,v->f(i,j-1,k+1,l)/*v[j-1][iref(i,j-1,N)][a+1][b+1]*/+2*data->eparam[6]+jdangle
                  +kdangle+2*data->eparam[10]+2*gap+
                  penalty(i,j-1,ct1,data)+penalty(k+1,l,ct2,data));

      }

      //case 8 - j, k, and l stacked:
      if (l!=N2+1&&jldecrement) {
        en1 = min(en1,v->f(i,j-1,k+1,l-1)/*v[j-1][iref(i,j-1,N)][a+1][b]*/+3*data->eparam[6]+jdangle
                  +kldangle+2*data->eparam[10]+gap
                  +penalty(i,j-1,ct1,data)+penalty(k+1,l-1,ct2,data));
      }

    }
							

  }

  if (l<=highend[j-1]&&l>=lowend[j-1]&&(j!=N+1)) {
    //case 5 - j stacked
    en1 = min(en1,v->f(i,j-1,k,l)/*v[j-1][iref(i,j-1,N)][a][b+1]*/+data->eparam[6]+jdangle+2*data->eparam[10]+gap
              +penalty(i,j-1,ct1,data)+penalty(k,l,ct2,data));
					
    //case 15 - i,j, and k stacked
    if ((j-1>i+1)&&(i!=N)) {
      if (k!=N2&&ikincrement) {
        en1 = min(en1,v->f(i+1,j-1,k+1,l)/*v[j-1][iref(i+1,j-1,N)][a][b+1]*/+3*data->eparam[6]+kdangle
                  +ijdangle
                  +2*data->eparam[10]+gap
                  +penalty(i+1,j-1,ct1,data)+penalty(k+1,l,ct2,data));
      }

								 
      //case 13 - i and j stacked

								
      if (k>=lowend[i+1]&&k<=highend[i+1]) {
        en1 = min(en1,v->f(i+1,j-1,k,l)/*v[j-1][iref(i+1,j-1,N)][a-1][b+1]*/+2*data->eparam[6]+ijdangle
                  +2*data->eparam[10]+2*gap
                  +penalty(i+1,j-1,ct1,data)+penalty(k,l,ct2,data));
      }

    }
						
  }

  if (k>=lowend[i+1]&&k<=highend[i+1]&&(i!=N)) {
    //case 9 - i alone is stacked
    en1 = min(en1,v->f(i+1,j,k,l)/*v[j][iref(i+1,j,N)][a-1][b]*/+data->eparam[6]+idangle+2*data->eparam[10]+gap
              +penalty(i+1,j,ct1,data)+penalty(k,l,ct2,data));


    //case 14 - i, j, and l stacked
    if ((j-1>i+1)&&(j!=N+1)&&(l!=N2+1)&&jldecrement) {
      en1=min(en1,v->f(i+1,j-1,k,l-1)/*v[j-1][iref(i+1,j-1,N)][a-1][b]*/+3*data->eparam[6]+ldangle
              +ijdangle
              +2*data->eparam[10]+gap
              +penalty(i+1,j-1,ct1,data)+penalty(k,l-1,ct2,data));
    }
  }


  if (l-1>=lowend[j]/*b-1>=0*/&&(i!=N)&&(l!=N2+1)) {
    //case 12 - i, k, and l stacked:
    if (k!=N2&&ikincrement) {
      en1=min(en1,v->f(i+1,j,k+1,l-1)/*v[j][iref(i+1,j,N)][a][b-1]*/+3*data->eparam[6]+kldangle
              +idangle+2*data->eparam[10]+gap
              +penalty(i+1,j,ct1,data)+penalty(k+1,l-1,ct2,data));
    }
    if (k>=lowend[i+1]&&k<=highend[i+1]) {
      //case 10 - l and i stacked
      en1=min(en1,v->f(i+1,j,k,l-1)/*v[j][iref(i+1,j,N)][a-1][b-1]*/+2*data->eparam[6]+ldangle
              +idangle+2*data->eparam[10]+2*gap+
              penalty(i+1,j,ct1,data)+penalty(k,l-1,ct2,data));
    }
  }
						
  //calculate the free energy of 2 fragments merged:
  for (c=i+minloop;c<j-minloop;++c) {
    //if (c>N) {
    //	kp = Ndiff;
    //}
    //else kp=0;
    startd=max(k+minloop,lowend[c]);
    endd=min(l-minloop,highend[c]);
    for (d=startd/*max(k+minloop,c-maxsep-kp)*/;
         d<=endd/*d<l-minloop&&d<=c+maxsep-kp*/;++d) {
				
      //e = d-c+maxsep+kp;
								
      if ((c!=N)&&(d!=N2)&&d+1>=lowend[c+1]&&d+1<=highend[c+1]) {
        if (c<N&&d<N2) {
          en1 = min(en1,w->f(i,c,k,d)+w->f(c+1,j,d+1,l)/*w[c][iref(i,c,N)][a][e]+w[j][iref(c+1,j,N)][e][b]*/);									
        }
        else {
          if (d>N2&&c>N) {
            en1 = min(en1,w->f(i,c,k,d)+w->f(c+1,j,d+1,l)/*w[c][iref(i,c,N)][a][e]+w[j-N][iref(c+1-N,j-N,N)][e][b]*/);
          }
        }


      }

    }
  }
					
  w->f(i,j,k,l)=en1;
  	
} //end of dynalignstep


//traceback a single conserved structure (either interior or exterior fragments)
//return an int that indicates whether an error occurred
int  dyntrace(short i, short j, short a, short b, structure *ct1, structure *ct2,
              short structnum, short *alignment,
              dynalignarray *w, varray *v, wendarray *w3, wendarray *w5, 
              short *lowend, short *highend,datatable *data, short gap, dynalignarray *vmod,
              bool local, bool startopen) {

  short int (*edangle5)(int i, int j, int ip, structure *ct, datatable *data) = edangle5noforce;
  short int (*edangle3)(int i, int j, int ip, structure *ct, datatable *data) = edangle3noforce;

  short k,l,c,d,en1,en2,e,f,kp,ip,jp,dp,lp;
  dynalignstackclass stack;
  bool open,closed;
  integersize en3;
  register bool found;
  int constantclosure;
  bool ikincrement,jldecrement,ikincrement2,jldecrement2;
  int error = 0;
	
	
  //modification tracks whether chemical modification considerations apply
  bool modification,helicalstack;

  modification = ct1->nmod > 0 || ct2->nmod > 0;

  short int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
                       {0,1,0,1,0,0},{0,0,0,0,0,0}};

  //put the whole fragment on the stack
  //a and b originally referred to positions in the large arrays -
  //now a and b refer to the nucleotide positions
  if (startopen) stack.push(j,b,1,0,w5->f(j,b),startopen); 
  else stack.push(i,j,a,b,v->f(i,j,a,b),startopen);

  while (stack.pull(&i,&j, &k, &l /*&a,&b*/, &en3, &open)) {
    //cout << "pulled "<<i<< " "<<j<<" "<<a<<" "<<b<<"\n"<<flush;

		
    if (!open) {
      if (modification) {
        if (en3==v->f(i,j,k,l)/*v[j][iref(i,j,ct1->numofbases)][a][b]*/) closed = true;
        else if (en3==vmod->f(i,j,k,l)/*vmod[j][iref(i,j,ct1->numofbases)][a][b]*/) closed = true;
        else closed = false;

      }
      else {

        if (en3==v->f(i,j,k,l)) closed = true;
        else closed = false;

      }

    }

    if (open) {
      //This is from w3 or w5:
      if (k==ct1->numofbases) {
        ip =i;
        kp = j;
        a = i;
        b = j;


		

        //ap = j;
        //we are dealing with a fragment ending in ct1->numofbases, i.e. w3
		while (ip<ct1->numofbases) {								
			  found = false;
			//First check if there is no structure in the 3' ends.
			//This stopping rule is needed because the normal whittle away at the ends might
				//not work with the lowend/highend scheme for allowing aligned nucleotides.
			if (local) {
				  //local alignment calculation
				if(en3==0) {
					found = true;
					ip=ct1->numofbases;

				}
					
			}
			else {
				 //global alignment calculation
				if (en3==gap*(abs((ct1->numofbases-ip)-(ct2->numofbases-kp)))) {
					found =true;
					ip=ct1->numofbases;

				}
					
			}


        

          //check if it is safe to increment i and k
          ikincrement =
            kp+1 >= lowend[ip+1] &&
            kp+1 <= highend[ip+1];

          //kp = ip + ap -maxsep;
          if (ikincrement&&!found) {
            if (en3 == w3->f(ip+1,kp+1)) {//adding one nuc on both that is unpaired and unstacked
              found = true;	
              en3 = w3->f(ip+1,kp+1);
              ip++;
              kp++;//

            }
          }
					
          if (kp>=lowend[ip+1]&&kp<=highend[ip+1]&&!found){
            if (en3 == w3->f(ip+1,kp)+gap) {
              found = true;
              en3 = w3->f(ip+1,kp);
              ip++;
              

            }

          }


          if (kp+1<=highend[ip]/*highlimit(ip,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found) {
            if (en3==w3->f(ip,kp+1)+gap){
              found = true;
              en3 = w3->f(ip,kp+1);
              kp++;

            }
          }

          if (local) {

            if (kp>=lowend[ip+1]&&kp<=highend[ip+1]/*lowlimit(ip+1,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found){
              if (en3 == w3->f(ip+1,kp)) {
                found = true;
                en3 = w3->f(ip+1,kp);
                ip++;
								

              }

            }


            if (kp+1<=highend[ip]/*highlimit(ip,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found) {
              if (en3==w3->f(ip,kp+1)){
                found = true;
                en3 = w3->f(ip,kp+1);
                kp++;

              }
            }

          }

			
			
          for (jp=ct1->numofbases;jp>=ip+minloop;jp--) {
            dp = max(lowend[jp]/*lowlimit(jp,maxsep,ct1->numofbases,ct2->numofbases)*/,kp+minloop+1);
            for (lp=min(ct2->numofbases,highend[jp]/*highlimit(jp,maxsep,ct1->numofbases,ct2->numofbases)*/);lp>=dp&&!found;lp--) {
              //bp = lp-jp+maxsep;	

              //check whether mine[i][a] is split so that the lowest free energy is
              //an exterior fragment from 1 to j and a helix from j+1 to i; 

              //check all possible alignments (in index a(k) and b(l))
					
              //must consider whether an unpaired nucleotide is stacked
              //is stacked onto each of the nucleotides in a helix
              //There are 16 cases:
              //consider them in this order (0 unstacked, 1 stacked)
              //		jp+1	ip	lp+1	k
              //1		0		0	0		0
              //2		0		0	0		1
              //3		0		0	1		0
              //4		0		0	1		1
              //5		0		1	0		0
              //6		0		1	0		1
              //7		0		1	1		0
              //8		0		1	1		1
              //9		1		0	0		0
              //10		1		0	0		1
              //11		1		0	1		0
              //12		1		0	1		1
              //13		1		1	0		0
              //14		1		1	0		1
              //15		1		1	1		0
              //16		1		1	1		1

              //note that for these exterior loops:
              //j<i
              //l<k
							
              if (lp+1>=lowend[jp+1]&&lp+1<=highend[jp+1]/*lowlimit(jp+1,maxsep,ct1->numofbases,ct2->numofbases)*/) {

                jldecrement =
                  lp-1 <= highend[jp-1]/*highlimit(jp-1, maxsep, ct1->numofbases, ct2->numofbases)*/ &&
                  lp-1 >= lowend[jp-1]/*lowlimit(jp-1, maxsep, ct1->numofbases, ct2->numofbases)*/;

                //no stacking
                //case 1:
                if (en3==w3->f(jp+1,lp+1)+v->f(ip,jp,kp,lp)+
                    penalty(ip,jp,ct1,data)+penalty(kp,lp,ct2,data)) {
								
                  found = true;
                  stack.push(ip,jp,kp,lp,v->f(ip,jp,kp,lp));
                  en3=w3->f(jp+1,lp+1);
                  ip=jp+1;
                  kp=lp+1;;

                }

                //case 6:
                if (!found&&ikincrement) if (en3==w3->f(jp+1,lp+1)+v->f(ip+1,jp,kp+1,lp)+edangle5(ip+1,jp,ip,ct1,data)+
                                             edangle5(kp+1,lp,kp,ct2,data)+
                                             penalty(ip+1,jp,ct1,data)+penalty(kp+1,lp,ct2,data)){
				
                  found = true;
                  stack.push(ip+1,jp,kp+1,lp,v->f(ip+1,jp,kp+1,lp));
                  en3 = w3->f(jp+1,lp+1);
                  ip=jp+1;
                  kp=lp+1;


                }

                //case 11	1	0	1	0
                if (!found&&jldecrement) if(en3==w3->f(jp+1,lp+1)+v->f(ip,jp-1,kp,lp-1)+edangle3(jp-1,ip,jp,ct1,data)+
                                            edangle3(lp-1,kp,lp,ct2,data)+penalty(ip,jp-1,ct1,data)+penalty(lp-1,kp,ct2,data)) {

                  found = true;
                  stack.push(ip,jp-1,kp,lp-1,v->f(ip,jp-1,kp,lp-1));
                  en3 = w3->f(jp+1,lp+1);
                  ip = jp+1;
                  kp=lp+1;

                }


                //case 16	1	1	1	1
                if (!found&&ikincrement&&jldecrement) if(en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp+1,lp-1)/*v[jp-1][ip+1][ap][bp]*/+edangle3(jp-1,ip+1,jp,ct1,data)+
                                                         edangle5(ip+1,jp-1,ip,ct1,data)+edangle3(lp-1,kp+1,lp,ct2,data)+
                                                         edangle5(kp+1,lp-1,kp,ct2,data)+penalty(ip+1,jp-1,ct1,data)+penalty(lp-1,kp+1,ct2,data)) {

                  found = true;
                  stack.push(ip+1,jp-1,kp+1,lp-1,v->f(ip+1,jp-1,kp+1,lp-1));
                  en3 = w3->f(jp+1,lp+1);
                  ip = jp+1;
                  kp=lp+1;

                }

					
                if (kp+1<=highend[ip]/*highlimit(ip,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found) {
                  //case 2:
                  if (en3==w3->f(jp+1,lp+1)+v->f(ip,jp,kp+1,lp)+edangle5(kp+1,lp,kp,ct2,data)+
                      penalty(ip,jp,ct1,data)+penalty(kp+1,lp,ct2,data)+gap){

                    found = true;
                    stack.push(ip,jp,kp+1,lp,v->f(ip,jp,kp+1,lp));
                    en3=w3->f(jp+1,lp+1);
                    ip = jp+1;
                    kp=lp+1;

                  }

                  //case 12	1	0	1	1
                  if (!found&&jldecrement) if(en3==w3->f(jp+1,lp+1)+v->f(ip,jp-1,kp+1,lp-1)+edangle3(jp-1,ip,jp,ct1,data)+
                                              edangle3(lp-1,kp,lp,ct2,data) + edangle5(kp+1,lp-1,kp,ct2,data)+
                                              penalty(ip,jp-1,ct1,data)+penalty(kp+1,lp-1,ct2,data)+gap){

                    found = true;
                    stack.push(ip,jp-1,kp+1,lp-1,v->f(ip,jp-1,kp+1,lp-1));
                    en3 = w3->f(jp+1,lp+1);
                    ip = jp+1;
                    kp=lp+1;


                  }

                  if (lp-1>=lowend[jp]/*lowlimit(jp,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found) {
                    //case 4
                    if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip,jp,kp+1,lp-1)/*v[jp][ip][ap+1][bp-1]*/+edangle5(kp+1,lp-1,kp,ct2,data)+
                        edangle3(lp-1,kp+1,lp,ct2,data)+penalty(ip,jp,ct1,data)+penalty(kp+1,lp-1,ct2,data)+2*gap) {
		
                      found = true;
                      stack.push(ip,jp,kp+1/*ap+1*/,lp-1/*bp-1*/,v->f(ip,jp,kp+1,lp-1)/*v[jp][ip][ap+1][bp-1]*/);
                      en3 = w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                      ip = jp+1;
                      kp=lp+1;//ap = bp;

                    }
                  }
	
                }

					
                if (lp-1>=lowend[jp]/*lowlimit(jp,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found) {
                  //case 3:
                  if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip,jp,kp,lp-1)/*v[jp][ip][ap][bp-1]*/+edangle3(lp-1,kp,lp,ct2,data)+
                      penalty(ip,jp,ct1,data)+penalty(lp-1,kp,ct2,data)+gap) {

                    found = true;
                    stack.push(ip,jp,kp/*ap*/,lp-1/*bp-1*/,v->f(ip,jp,kp,lp-1)/*v[jp][ip][ap][bp-1]*/);
                    en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                    ip = jp+1;
                    kp=lp+1;//ap = bp;

                  }

                  //case 8:		0	1	1	1
                  if (!found&&ikincrement) if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip+1,jp,kp+1,lp-1)/*v[jp][ip+1][ap][bp-1]*/+edangle5(ip+1,jp,ip,ct1,data)+
                                               edangle3(lp-1,kp+1,lp,ct2,data)+edangle5(kp+1,lp-1,kp,ct2,data)+
                                               penalty(ip+1,jp,ct1,data)+penalty(lp-1,kp+1,ct2,data)+gap) {

                    found = true;
                    stack.push(ip+1,jp,kp+1/*ap*/,lp-1/*bp-1*/,v->f(ip+1,jp,kp+1,lp-1)/*v[jp][ip+1][ap][bp-1]*/);
                    en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                    ip = jp+1;
                    kp=lp+1;//ap = bp;


                  }

                  if (kp>=lowend[ip+1]&&kp<=highend[ip+1]&&!found) {
                    //case 7		0	1	1	0
                    if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip+1,jp,kp,lp-1)/*v[jp][ip+1][ap-1][bp-1]*/+edangle3(lp-1,kp,lp,ct2,data)+
                        edangle5(ip+1,jp,ip,ct1,data)+penalty(ip+1,jp,ct1,data)+penalty(lp-1,kp,ct2,data)+2*gap) {

                      found = true;
                      stack.push(ip+1,jp,kp/*ap-1*/,lp-1/*bp-1*/,v->f(ip+1,jp,kp,lp-1)/*v[jp][ip+1][ap-1][bp-1]*/);
                      en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                      ip = jp+1;
                      kp=lp+1;//ap = bp;

                    }
                  }


                }


                if (kp>=lowend[ip+1]&&kp<=highend[ip+1]&&!found) {
                  //case5: 0 1 0 0
                  if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip+1,jp,kp,lp)/*v[jp][ip+1][ap-1][bp]*/+edangle5(ip+1,jp,ip,ct1,data)+
                      penalty(ip+1,jp,ct1,data)+penalty(kp,lp,ct2,data)+gap) {

                    found = true;
                    stack.push(ip+1,jp,kp/*ap-1*/,lp/*bp*/,v->f(ip+1,jp,kp,lp)/*v[jp][ip+1][ap-1][bp]*/);
                    en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                    ip = jp+1;
                    kp=lp+1;//ap = bp;

                  }

                  //case 15	1	1	1	0
                  if (!found&&jldecrement) if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp,lp-1)/*v[jp-1][ip+1][ap-1][bp]*/ + edangle3(jp-1,ip+1,jp,ct1,data)+
                                               edangle5(ip+1,jp-1,ip,ct1,data)+edangle3(lp-1,kp,lp,ct2,data)+
                                               penalty(ip+1,jp-1,ct1,data)+penalty(lp-1,kp,ct2,data)+gap) {

                    found = true;
                    stack.push(ip+1,jp-1,kp/*ap-1*/,lp-1/*bp*/,v->f(ip+1,jp-1,kp,lp-1)/*v[jp-1][ip+1][ap-1][bp]*/);
                    en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                    ip = jp+1;
                    kp=lp+1;//ap = bp;
	
                  }

                  if (lp<=highend[jp-1]&&lp>=lowend[jp-1]/*highlimit(jp-1,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found) {
                    //case 13	1	1	0	0
                    if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp,lp)/*v[jp-1][ip+1][ap-1][bp+1]*/+edangle3(jp-1,ip+1,jp,ct1,data)+
                        edangle5(ip+1,jp-1,ip,ct1,data)+
                        penalty(ip+1,jp-1,ct1,data)+penalty(kp,lp,ct2,data)+2*gap) {


                      found = true;
                      stack.push(ip+1,jp-1,/*ap-1*/kp,lp/*bp+1*/,v->f(ip+1,jp-1,kp,lp)/*v[jp-1][ip+1][ap-1][bp+1]*/);
                      en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                      ip = jp+1;
                      kp=lp+1;//ap = bp;

                    }

                  }


                }
                if (lp<=highend[jp-1]&&lp>=lowend[jp-1]/*highlimit(jp-1,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found) {
                  //case 9		1	0	0	0
                  if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip,jp-1,kp,lp)/*v[jp-1][ip][ap][bp+1]*/+edangle3(jp-1,ip,jp,ct1,data)+
                      penalty(ip,jp-1,ct1,data)+penalty(kp,lp,ct2,data)+gap) {

                    found = true;
                    stack.push(ip,jp-1,kp/*ap*/,lp/*bp+1*/,v->f(ip,jp-1,kp,lp)/*v[jp-1][ip][ap][bp+1]*/);
                    en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                    ip = jp+1;
                    kp=lp+1;//ap = bp;

                  }

					
                  //case 14	1	1	0	1
                  if (!found&&ikincrement) if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp+1,lp)/*v[jp-1][ip+1][ap][bp+1]*/+edangle3(jp-1,ip+1,jp,ct1,data)+
                                               edangle5(ip+1,jp-1,ip,ct1,data)+edangle5(kp+1,lp,kp,ct2,data)+
                                               penalty(ip+1,jp-1,ct1,data)+penalty(kp+1,lp,ct2,data)+gap) {

                    found = true;
                    stack.push(ip+1,jp-1,kp+1/*ap*/,lp/*bp+1*/,v->f(ip+1,jp-1,kp+1,lp)/*v[jp-1][ip+1][ap][bp+1]*/);
                    en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                    ip=jp+1;
                    kp=lp+1;//ap = bp;

                  }


                  if (kp+1<=highend[ip]/*highlimit(ip,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found) {
                    //case 10	1	0	0	1
                    if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip,jp-1,kp+1,lp)/*v[jp-1][ip][ap+1][bp+1]*/+edangle3(jp-1,ip,jp,ct1,data)+
                        edangle5(kp+1,lp,kp,ct2,data)+penalty(ip,jp-1,ct1,data)+penalty(kp+1,lp,ct2,data)+2*gap) {


                      found = true;
                      stack.push(ip,jp-1,kp+1/*ap+1*/,lp/*bp+1*/,v->f(ip,jp-1,kp+1,lp)/*v[jp-1][ip][ap+1][bp+1]*/);
                      en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                      ip=jp+1;
                      kp=lp+1;//ap = bp;

                    }

                  }


					
                }


              }
            }
          }
          if (!found) {
		    //A traceback error occurred
            //cerr << "Traceback error at w3\n";
		    error = 14;
            break;

          }

        }
									

      }
      else {
        //dealing with w5

        //a = j;
        k = j;
        while(i>0) {
          found = false;
          //k = a+i-maxsep;




          //although this is really a decrement...
          ikincrement =
            k-1 <= highend[i-1]/*highlimit(i-1, maxsep, ct1->numofbases, ct2->numofbases)*/ &&
            k-1 >= lowend[i-1]/*lowlimit(i-1, maxsep, ct1->numofbases, ct2->numofbases)*/;

		  //check if there is no structure in the 5' direction.
		  //This is important because the usual whittle away might not work with the highlimit/lowlimit tracking
			//of what can be aligned...
		  if (local) {
			
			  if (en3==0) {

				  //If energy is zero, no need to traceback farther
				  found = true;
				  i=0;
			  }
			
		  } 
		  else {
									
			
					
			  if (en3 == gap* (abs(i - k))) {

				  //If energy is accounted with gaps, then no need for more traceback
				  found = true;
				  i=0;

			  }
			
		  }

		
          if (ikincrement&&!found) { 
            if (en3==w5->f(i-1,k-1)/*w5[i-1][a]*/) {
              found = true;
              en3 = w5->f(i-1,k-1)/*w5[i-1][a]*/;
              i--;
              k--;//

            }
          }

		  if (!found) {
			  if (k<=highend[i-1]&&k>=lowend[i-1]) { 
				if (en3==w5->f(i-1,k)/*w5[i-1][a+1]*/+gap){
				  found = true;
				  en3 = w5->f(i-1,k);/*w5[i-1][a+1]*/
				  i--;
				  //a++;
	

				}

			  }
		  }

					

          if (k-1>=lowend[i]&&k-1<=highend[i]&&!found) {
		
            if(en3==w5->f(i,k-1)/*w5[i][a-1]*/+gap) {
              found = true;
              en3 = w5->f(i,k-1)/*w5[i][a-1]*/;
              k--;//a--;

            }

          }

          else if (local) {
            if (k<=highend[i-1]&&k>=lowend[i-1]&&!found) { 
              if (en3==w5->f(i-1,k)){
                found = true;
                en3 = w5->f(i-1,k);
                i--;
								


              }

            }

					

            if (k-1>=lowend[i]&&k-1<=highend[i]&&!found) {
		
              if(en3==w5->f(i,k-1)) {
                found = true;
                en3 = w5->f(i,k-1);
                k--;

              }

            }

          }

          for (j=0;j+minloop<i&&!found;j++) {
            for (l=max(0,lowend[j]/*lowlimit(j,maxsep,ct1->numofbases,ct2->numofbases)/*j-maxsep*/);l<=highend[j]/*highlimit(j,maxsep,ct1->numofbases,ct2->numofbases)/*j+maxsep*/&&l<=ct2->numofbases&&l+minloop<k&&!found;l++) {
              //b = l-j+maxsep;	

              //check whether w5[i][a] is split so that the lowest free energy is
              //an exterior fragment from 1 to j and a helix from j+1 to i; 

              //check all possible alignments (in index a(k) and b(l))
					
						

              //no stacking
              //if(en3== w5[j][b]+v[i][j+1][b][a]) {
              //	found = true;
							
              //	stack.push(j+1,i,b,a,v[i][j+1][b][a]);
              //	en3 = w5[j][b];
              //	i = j;
              //	a = b;


              //}

              //check whether w5[i][a] is split so that the lowest free energy is
              //an exterior fragment from 1 to j and a helix from j+1 to i; 

              //check all possible alignments (in index a(k) and b(l))
					
              //must consider whether an unpaired nucleotide is stacked
              //is stacked onto each of the nucleotides in a helix
              //There are 16 cases:
              //consider them in this order (0 unstacked, 1 stacked)
              //		j+1	i	l+1	k
              //1		0	0	0	0
              //2		0	0	0	1
              //3		0	0	1	0
              //4		0	0	1	1
              //5		0	1	0	0
              //6		0	1	0	1
              //7		0	1	1	0
              //8		0	1	1	1
              //9		1	0	0	0
              //10		1	0	0	1
              //11		1	0	1	0
              //12		1	0	1	1
              //13		1	1	0	0
              //14		1	1	0	1
              //15		1	1	1	0
              //16		1	1	1	1

              //note that for these exterior loops:
              //j<i
              //l<k

              //although this is really an increment...
              jldecrement =
                l+1 >= lowend[j+1]/*lowlimit(j+1, maxsep, ct1->numofbases, ct2->numofbases)*/ &&
                l+1 <= highend[j+1]/*highlimit(j+1, maxsep, ct1->numofbases, ct2->numofbases)*/;
              //although this is really an increment...
              jldecrement2 =
                l+2 >= lowend[j+2]/*lowlimit(j+2, maxsep, ct1->numofbases, ct2->numofbases)*/ &&
                l+2 <= highend[j+2]/*highlimit(j+2, maxsep, ct1->numofbases, ct2->numofbases)*/;


              //no stacking
              //case 1:
              if (jldecrement) {
                if (en3 == w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i,l+1,k)/*v[i][j+1][b][a]*/+penalty(i,j+1,ct1,data)+penalty(k,l+1,ct2,data)) {
                  found = true;
                  en3 = w5->f(j,l)/*w5[j][b]*/;
                  stack.push(j+1,i,l+1/*b*/,k/*a*/,v->f(j+1,i,l+1,k)/*v[i][j+1][b][a]*/);
                  i=j;
                  k=l;
                  //a=b;

                }
              }
              //case 6:
					
              if (!found&&ikincrement&&jldecrement) {
                if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i-1,l+1,k-1)/*v[i-1][j+1][b][a]*/+edangle3(i-1,j+1,i,ct1,data)+
                   edangle3(k-1,l+1,k,ct2,data)+penalty(i-1,j+1,ct1,data)+penalty(k-1,l+1,ct2,data)){

                  found = true;
                  en3 = w5->f(j,l)/*w5[j][b]*/;
                  stack.push(j+1,i-1,l+1/*b*/,k-1/*a*/,v->f(j+1,i-1,l+1,k-1)/*v[i-1][j+1][b][a]*/);
                  i = j;
                  k=l;
                  //a=b;
                }
              }

              //case 11	1	0	1	0
              if (!found&&jldecrement2) {
                if (en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i,l+2,k)/*v[i][j+2][b][a]*/+edangle5(j+2,i,j+1,ct1,data)+
                    edangle5(l+2,k,l+1,ct2,data)+penalty(i,j+2,ct1,data)+penalty(l+2,k,ct2,data)) {
								
                  stack.push(j+2,i,l+2/*b*/,k/*a*/,v->f(j+2,i,l+2,k)/*v[i][j+2][b][a]*/);
                  found = true;
                  en3=w5->f(j,l)/*w5[j][b]*/;
                  i = j;
                  k=l;
                  //a=b;

					
                }
              }


              //case 16	1	1	1	1
              if (!found&&ikincrement&&jldecrement2) {
                if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i-1,l+2,k-1)/*v[i-1][j+2][b][a]*/+edangle5(j+2,i-1,j+1,ct1,data)+
                   edangle3(i-1,j+2,i,ct1,data)+edangle5(l+2,k-1,l+1,ct2,data)+
                   edangle3(k-1,l+2,k,ct2,data)+penalty(i-1,j+2,ct1,data)+penalty(k-1,l+2,ct2,data)) {

                  stack.push(j+2,i-1,l+2/*b*/,k-1/*a*/,v->f(j+2,i-1,l+2,k-1)/*v[i-1][j+2][b][a]*/);
                  found = true;
                  en3=w5->f(j,l)/*w5[j][b]*/;
                  i = j;
                  k=l;
                  //a=b;

					

                }
              }

					
					
              if (!found&&k-1>=lowend[i]) {
                //case 2: //gap added
                if (jldecrement) {
                  if (en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i,l+1,k-1)/*v[i][j+1][b][a-1]*/+edangle3(k-1,l+1,k,ct2,data)
                      +penalty(i,j+1,ct1,data)+penalty(k-1,l+1,ct2,data)+gap) {
									
                    stack.push(j+1,i,l+1/*b*/,k-1/*a-1*/,v->f(j+1,i,l+1,k-1)/*v[i][j+1][b][a-1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;

						

                  }
                }

                //case 12	1	0	1	1
                if(jldecrement2&&!found) {
									
                  if (en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i,l+2,k-1)/*v[i][j+2][b][a-1]*/+edangle5(j+2,i,j+1,ct1,data)+
                      edangle5(l+2,k-1,l+1,ct2,data) + edangle3(k-1,l+2,k,ct2,data)+
                      penalty(i,j+2,ct1,data)+penalty(k-1,l+2,ct2,data)+gap) {

                    stack.push(j+2,i,l+2/*b*/,k-1/*a-1*/,v->f(j+2,i,l+2,k-1)/*v[i][j+2][b][a-1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;

						

                  }
                }

                if (!found&&l+2<=highend[j+1]/*highlimit(j+1,maxsep,ct1->numofbases,ct2->numofbases)*/&&(l+2>=lowend[j+1]/*lowlimit(j+1,maxsep,ct1->numofbases,ct2->numofbases)*/)) {
                  //case 4
                  if	(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i,l+2,k-1)/*v[i][j+1][b+1][a-1]*/+edangle3(k-1,l+2,k,ct2,data)+
                       edangle5(l+2,k-1,l+1,ct2,data)+penalty(i,j+1,ct1,data)+penalty(k-1,l+2,ct2,data)+2*gap) {
										
                    stack.push(j+1,i,l+2/*b+1*/,k-1/*a-1*/,v->f(j+1,i,l+2,k-1)/*v[i][j+1][b+1][a-1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;

						

                  }

                }

              }

					
              if (l+2<=highend[j+1]/*highlimit(j+1,maxsep,ct1->numofbases,ct2->numofbases)*/&&l+2>=lowend[j+1]/*lowlimit(j+1,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found) {
                //case 3:
                if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i,l+2,k)/*v[i][j+1][b+1][a]*/+edangle5(l+2,k,l+1,ct2,data)+
                   penalty(i,j+1,ct1,data)+penalty(l+2,k,ct2,data)+gap) {
									
                  stack.push(j+1,i,l+2/*b+1*/,k/*a*/,v->f(j+1,i,l+2,k)/*v[i][j+1][b+1][a]*/);
                  found = true;
                  en3=w5->f(j,l)/*w5[j][b]*/;
                  i = j;
                  k=l;
                  //a=b;
	


                }

                //case 8:		0	1	1	1
                if (ikincrement&&!found) {
                  if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i-1,l+2,k-1)/*v[i-1][j+1][b+1][a]*/+edangle3(i-1,j+1,i,ct1,data)+
                     edangle5(l+2,k-1,l+1,ct2,data)+edangle3(k-1,l+2,k,ct2,data)+
									
                     penalty(i-1,j+1,ct1,data)+penalty(k-1,l+2,ct2,data)+gap) {
                    stack.push(j+1,i-1,l+2/*b+1*/,k-1/*a*/,v->f(j+1,i-1,l+2,k-1)/*v[i-1][j+1][b+1][a]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;

						

                  }
                }

                if (!found&&k<=highend[i-1]&&k>=lowend[i-1]) {
                  //case 7		0	1	1	0
                  if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i-1,l+2,k)/*v[i-1][j+1][b+1][a+1]*/+edangle5(l+2,k,l+1,ct2,data)+
                     edangle3(i-1,j+1,i,ct1,data)+penalty(i-1,j+1,ct1,data)+penalty(l+2,k,ct2,data)+2*gap) {
										
                    stack.push(j+1,i-1,l+2/*b+1*/,k/*a+1*/,v->f(j+1,i-1,l+2,k)/*v[i-1][j+1][b+1][a+1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;

							

                  }
                }
              }

		
			if (k<=highend[i-1]&&k>=lowend[i-1]&&!found) {
                //case5:
                if (jldecrement) {
                  if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i-1,l+1,k)/*v[i-1][j+1][b][a+1]*/+edangle3(i-1,j+1,i,ct1,data)+
                     penalty(i-1,j+1,ct1,data)+penalty(k,l+1,ct2,data)+gap) {
									
                    stack.push(j+1,i-1,l+1/*b*/,k/*a+1*/,v->f(j+1,i-1,l+1,k)/*v[i-1][j+1][b][a+1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;

                  }

                }

                //case 15	1	1	1	0
                if (!found&&jldecrement2) {
                  if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i-1,l+2,k)/*v[i-1][j+2][b][a+1]*/ + edangle5(j+2,i-1,j+1,ct1,data)+
                     edangle3(i-1,j+2,i,ct1,data)+edangle5(l+2,k,l+1,ct2,data)+
                     penalty(i-1,j+2,ct1,data)+penalty(l+2,k,ct2,data)+gap) {

                    stack.push(j+2,i-1,l+2/*b*/,k/*a+1*/,v->f(j+2,i-1,l+2,k)/*v[i-1][j+2][b][a+1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;

						

                  }
                }

                if (!found&&l+1>=lowend[j+2]&&l+1<=highend[j+2]/*lowlimit(j+2,maxsep,ct1->numofbases,ct2->numofbases)*/) {
                  //case 13	1	1	0	0
                  if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i-1,l+1,k)/*v[i-1][j+2][b-1][a+1]*/+edangle5(j+2,i-1,j+1,ct1,data)+
                     edangle3(i-1,j+2,i,ct1,data)+penalty(i-1,j+2,ct1,data)+penalty(k,l+1,ct2,data)+2*gap) {
										
                    stack.push(j+2,i-1,l+1/*b-1*/,k/*a+1*/,v->f(j+2,i-1,l+1,k)/*v[i-1][j+2][b-1][a+1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;

							


                  }
					
                }

              }
              if (l+1>=lowend[j+2]&&l+1<=highend[j+2]/*lowlimit(j+2,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found) {
                //case 9		1	0	0	0
                if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i,l+1,k)/*v[i][j+2][b-1][a]*/+edangle5(j+2,i,j+1,ct1,data)+
                   penalty(i,j+2,ct1,data)+penalty(k,l+1,ct2,data)+gap) {
									
                  stack.push(j+2,i,l+1/*b-1*/,k/*a*/,v->f(j+2,i,l+1,k)/*v[i][j+2][b-1][a]*/);
                  found = true;
                  en3=w5->f(j,l)/*w5[j][b]*/;
                  i = j;
                  k=l;
                  //a=b;

						


                }
					
                //case 14	1	1	0	1
                else if (ikincrement) {
                  if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i-1,l+1,k-1)/*v[i-1][j+2][b-1][a]*/+edangle5(j+2,i-1,j+1,ct1,data)+
                     edangle3(i-1,j+2,i,ct1,data)+edangle3(k-1,l+1,k,ct2,data)+
                     penalty(i-1,j+2,ct1,data)+penalty(k-1,l+1,ct2,data)+gap) {

                    stack.push(j+2,i-1,l+1/*b-1*/,k-1/*a*/,v->f(j+2,i-1,l+1,k-1)/*v[i-1][j+2][b-1][a]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;
	
						

                  }
                }

                if (!found&&k-1>=lowend[i]/*lowlimit(i,maxsep,ct1->numofbases,ct2->numofbases)*/) {
                  //case 10	1	0	0	1
                  if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i,l+1,k-1)/*v[i][j+2][b-1][a-1]*/+edangle5(j+2,i,j+1,ct1,data)+
                     edangle3(k-1,l+1,k,ct2,data)+penalty(i,j+2,ct1,data)+penalty(k-1,l+1,ct2,data)+2*gap) {

                    stack.push(j+2,i,l+1/*b-1*/,k-1/*a-1*/,v->f(j+2,i,l+1,k-1)/*v[i][j+2][b-1][a-1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;

							

                  }
                }
					
              }
					

            }

		  }
          if (!found) {
				//a traceback error occurred
            ///cerr << "Traceback error at w5!\n";
			  error = 14;
            break;

          }

        }
      }

    }
    else if (closed) {

      //check if it is safe to increment i and k
      ikincrement =
        k+1 >= lowend[i+1]/*lowlimit(i+1, maxsep, ct1->numofbases, ct2->numofbases)*/ &&
        k+1 <= highend[i+1]/*highlimit(i+1, maxsep, ct1->numofbases, ct2->numofbases)*/;

      ikincrement2 =
        k+2 >= lowend[i+2]/*lowlimit(i+2, maxsep, ct1->numofbases, ct2->numofbases)*/ &&
        k+2 <= highend[i+2]/*highlimit(i+2, maxsep, ct1->numofbases, ct2->numofbases)*/;

      //check to see if it is safe to decrement j and l
      jldecrement =
        l-1 <= highend[j-1]/*highlimit(j-1, maxsep, ct1->numofbases, ct2->numofbases)*/ &&
        l-1 >= lowend[j-1]/*lowlimit(j-1, maxsep, ct1->numofbases, ct2->numofbases)*/;

      jldecrement2 =
        l-2 <= highend[j-2]/*highlimit(j-2, maxsep, ct1->numofbases, ct2->numofbases)*/ &&
        l-2 >= lowend[j-2]/*lowlimit(j-2, maxsep, ct1->numofbases, ct2->numofbases)*/;


      //i and j are paired and aligned to a and b
			
      if (j<=ct1->numofbases) {
        ct1->basepr[ct1->numofstructures][i]=j;
        ct1->basepr[ct1->numofstructures][j]=i;
        ct2->basepr[ct2->numofstructures][k/*a+i-maxsep*/]=l/*b+j-maxsep*/;
        ct2->basepr[ct2->numofstructures][l/*b+j-maxsep*/]=k/*a+i-maxsep*/;
        alignment[i]=k/*a+i-maxsep*/;
        alignment[j]=l/*b+j-maxsep*/;
      }
      else {
        //j (and l) are > N (and N2)
        ct1->basepr[ct1->numofstructures][i]=j-ct1->numofbases;
        ct1->basepr[ct1->numofstructures][j-ct1->numofbases]=i;
        ct2->basepr[ct2->numofstructures][k]=l-ct2->numofbases;
        ct2->basepr[ct2->numofstructures][l-ct2->numofbases]=k;
        alignment[i]=k;
        alignment[j-ct1->numofbases]=l-ct2->numofbases;

      }

      //now find the next pair:

      if (en3!=erg3(i,j,ct1,data,0)+erg3(k,l,ct2,data,0)+gap*abs(j-i-l+k)){
        //the fragment does not close hairpins, therefore, the internal
        // fragment does need to be characterized

        //first check for internal loop:
        found= false;

        //Now consider internal loops/one stacked and other not
        for (c=i+1;c<=i+maxloop&&c<j&&!found&&c<=ct1->numofbases;c++) {
          for (d=j-1;d>=j-maxloop&&d>c;d--) {
            //if (d>ct1->numofbases) ap = ct1->numofbases-ct2->numofbases;
            //else ap =0;
            for (e=max(k+1,lowend[c]/*lowlimit(c,maxsep,ct1->numofbases,ct2->numofbases)*/);
                 e<=k+maxloop&&e<l&&e<=highend[c]/*highlimit(c,maxsep,ct1->numofbases,ct2->numofbases)*/&&e<=ct2->numofbases&&!found;e++) {
							
              for (f=min(l-1,highend[d]/*highlimit(d,maxsep,ct1->numofbases,ct2->numofbases)*/);f>=l-maxloop&&f>e&&f>=lowend[d]/*lowlimit(d,maxsep,ct1->numofbases,ct2->numofbases)*/
                     /*d-maxsep-ap*/&&!found;f--) {

                if (c-i+j-d-2<maxinternal&&e-k+l-f-2<maxinternal&&((d>ct1->numofbases)&&(f>ct2->numofbases))||((d<=ct1->numofbases)&&(f<=ct2->numofbases))) {						
                  helicalstack=false;		
                  if (c==i+1&&d==j-1&&(i!=ct1->numofbases)&&(j-1!=ct1->numofbases)&&
                      ((d>ct1->numofbases&&f>ct2->numofbases)||(d<=ct1->numofbases&&f<=ct2->numofbases))) {
                    //seq1 is helical stack
                    en1 = erg1(i,j,i+1,j-1,ct1, data);
                    if (modification) helicalstack=true;
									

                  }
                  else {
                    //internal loop

                    en1 = erg2(i,j,c,d,ct1,data,0,0);

                  }

                  if (e==k+1&&f==l-1&&k!=ct2->numofbases&&(l-1!=ct2->numofbases)) {
                    //seq2 is helical stack
                    en2 = erg1(k,l,k+1,l-1,ct2, data);
                    if (modification) helicalstack=true;
                  }
                  else {
                    //internal loop
                    en2 = erg2(k,l,e,f,ct2,data,0,0);

                  }
									
                  if (helicalstack) {
                    if(en3== en1+en2+vmod->f(c,d,e,f)/*vmod[d][iref(c,d,ct1->numofbases)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l))){
                      found = true;
                      stack.push(c,d,e/*-c+maxsep*/,f/*-d+maxsep+ap*/,vmod->f(c,d,e,f)/*vmod[d][iref(c,d,ct1->numofbases)][e-c+maxsep][ f-d+maxsep+ap]*/);
	
                    }	


                    else if (c==i+2&&d==j-2&&j-2>i+2&&
                             inc[ct1->numseq[i+1]][ct1->numseq[j-1]]&&
                             inc[ct1->numseq[i+2]][ct1->numseq[j-2]]&&
                             e==k+1&&f==l-1&&
                             inc[ct2->numseq[k+1]][ct2->numseq[l-1]]
                             &&i!=ct1->numofbases&&i+1!=ct1->numofbases&&
                             (j-1!=ct1->numofbases)&&(j-2!=ct1->numofbases)
                             &&k!=ct2->numofbases&&(l-1!=ct2->numofbases)) {

                      en1 = erg1(i,j,i+1,j-1,ct1,data)+
                        erg1(i+1,j-1,i+2,j-2,ct1,data);

                      if (en3==en1+en2+vmod->f(c,d,e,f)/*vmod[d][iref(c,d,ct1->numofbases)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l))) {
                        //base pairs with single base pair insertion in ct1
                        //ct1->basepr[ct1->numofstructures][i+1]=j-1;
                        //ct1->basepr[ct1->numofstructures][j-1]=i+1;
                        if (j<=ct1->numofbases) {
                          ct1->basepr[ct1->numofstructures][i+1]=j-1;
                          ct1->basepr[ct1->numofstructures][j-1]=i+1;
												
												
                        }
                        else {
                          //j (and l) are > N (and N2)
                          ct1->basepr[ct1->numofstructures][i+1]=j-1-ct1->numofbases;
                          ct1->basepr[ct1->numofstructures][j-1-ct1->numofbases]=i+1;
												
												

                        }
                        stack.push(c,d,e/*-c+maxsep*/,f/*-d+maxsep+ap*/,v->f(c,d,e,f)/*v[d][iref(c,d,ct1->numofbases)][e-c+maxsep][f-d+maxsep+ap]*/);
                        found = true;


                      }

											


                    }
                    else if (e==k+2&&f==l-2&&l-2>k+2&&c==i+1&&d==j-1&&
                             inc[ct2->numseq[k+1]][ct2->numseq[l-1]]&&
                             inc[ct2->numseq[k+2]][ct2->numseq[l-2]]&&
                             inc[ct1->numseq[i+1]][ct1->numseq[j-1]]&&
                             i!=ct1->numofbases&&(j-1!=ct1->numofbases)&&k!=ct2->numofbases&&(l-1!=ct2->numofbases)&&
                             (k+1!=ct2->numofbases)&&(l-2!=ct2->numofbases)) {

                      en2= erg1(k,l,k+1,l-1,ct2,data)+
                        erg1(k+1,l-1,k+2,l-2,ct2,data);

                      if (en3==en1+en2+vmod->f(c,d,e,f)/*vmod[d][iref(c,d,ct1->numofbases)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l))) {
                        //base pairs with single base pair insertion in ct1
                        //ct2->basepr[ct2->numofstructures][k+1]=l-1;
                        //ct2->basepr[ct2->numofstructures][l-1]=k+1;
                        if (l<=ct2->numofbases) {
												
                          ct2->basepr[ct2->numofstructures][k+1]=l-1;
                          ct2->basepr[ct2->numofstructures][l-1]=k+1;
												
                        }
                        else {
                          //j (and l) are > N (and N2)
												
                          ct2->basepr[ct2->numofstructures][k+1]=l-1-ct2->numofbases;
                          ct2->basepr[ct2->numofstructures][l-1-ct2->numofbases]=k+1;
												

                        }
                        stack.push(c,d,e/*-c+maxsep*/,f/*-d+maxsep+ap*/,v->f(c,d,e,f)/*v[d][iref(c,d,ct1->numofbases)][e-c+maxsep][f-d+maxsep+ap]*/);
										
                        found = true;

                      }

										
                    }
									

                  }
                  else {
                    if(en3== en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,ct1->numofbases)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l))){
                      found = true;
                      stack.push(c,d,e/*-c+maxsep*/,f/*-d+maxsep+ap*/,v->f(c,d,e,f)/*v[d][iref(c,d,ct1->numofbases)][e-c+maxsep][ f-d+maxsep+ap]*/);
									


                    }
								
                    else if (c==i+2&&d==j-2&&j-2>i+2&&
                             inc[ct1->numseq[i+1]][ct1->numseq[j-1]]&&
                             inc[ct1->numseq[i+2]][ct1->numseq[j-2]]&&
                             e==k+1&&f==l-1&&
                             inc[ct2->numseq[k+1]][ct2->numseq[l-1]]
                             &&i!=ct1->numofbases&&i+1!=ct1->numofbases&&
                             (j-1!=ct1->numofbases)&&(j-2!=ct1->numofbases)
                             &&k!=ct2->numofbases&&(l-1!=ct2->numofbases)) {

                      en1 = erg1(i,j,i+1,j-1,ct1,data)+
                        erg1(i+1,j-1,i+2,j-2,ct1,data);

                      if (en3==en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,ct1->numofbases)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l))) {
                        //base pairs with single base pair insertion in ct1
                        //ct1->basepr[ct1->numofstructures][i+1]=j-1;
                        //ct1->basepr[ct1->numofstructures][j-1]=i+1;
                        if (j<=ct1->numofbases) {
                          ct1->basepr[ct1->numofstructures][i+1]=j-1;
                          ct1->basepr[ct1->numofstructures][j-1]=i+1;
												
												
                        }
                        else {
                          //j (and l) are > N (and N2)
                          ct1->basepr[ct1->numofstructures][i+1]=j-1-ct1->numofbases;
                          ct1->basepr[ct1->numofstructures][j-1-ct1->numofbases]=i+1;
												
												

                        }
                        stack.push(c,d,e/*-c+maxsep*/,f/*-d+maxsep+ap*/,v->f(c,d,e,f)/*v[d][iref(c,d,ct1->numofbases)][e-c+maxsep][f-d+maxsep+ap]*/);
                        found = true;


                      }

											


                    }
                    else if (e==k+2&&f==l-2&&l-2>k+2&&c==i+1&&d==j-1&&
                             inc[ct2->numseq[k+1]][ct2->numseq[l-1]]&&
                             inc[ct2->numseq[k+2]][ct2->numseq[l-2]]&&
                             inc[ct1->numseq[i+1]][ct1->numseq[j-1]]&&
                             i!=ct1->numofbases&&(j-1!=ct1->numofbases)&&k!=ct2->numofbases&&(l-1!=ct2->numofbases)&&
                             (k+1!=ct2->numofbases)&&(l-2!=ct2->numofbases)) {

                      en2= erg1(k,l,k+1,l-1,ct2,data)+
                        erg1(k+1,l-1,k+2,l-2,ct2,data);

                      if (en3==en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,ct1->numofbases)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l))) {
                        //base pairs with single base pair insertion in ct1
                        //ct2->basepr[ct2->numofstructures][k+1]=l-1;
                        //ct2->basepr[ct2->numofstructures][l-1]=k+1;
                        if (l<=ct2->numofbases) {
												
                          ct2->basepr[ct2->numofstructures][k+1]=l-1;
                          ct2->basepr[ct2->numofstructures][l-1]=k+1;
												
                        }
                        else {
                          //j (and l) are > N (and N2)
												
                          ct2->basepr[ct2->numofstructures][k+1]=l-1-ct2->numofbases;
                          ct2->basepr[ct2->numofstructures][l-1-ct2->numofbases]=k+1;
												

                        }
                        stack.push(c,d,e/*-c+maxsep*/,f/*-d+maxsep+ap*/,v->f(c,d,e,f)/*v[d][iref(c,d,ct1->numofbases)][e-c+maxsep][f-d+maxsep+ap]*/);
										
                        found = true;

                      }

										
                    }
                  }
                }

              }
            }
          }
        }
        constantclosure = penalty(i,j,ct1,data)+penalty(k,l,ct2,data);
        if (j>ct1->numofbases&&!found) {
          //consider an exterior loop closed by i-j and k-l
          //There are 16 cases:
          //consider them in this order (0 unstacked, 1 stacked)
          //		i	j	k	l
          //1		0	0	0	0
          //2		0	0	0	1
          //3		0	0	1	0
          //4		0	0	1	1
          //5		0	1	0	0
          //6		0	1	0	1
          //7		0	1	1	0
          //8		0	1	1	1
          //9		1	0	0	0
          //10		1	0	0	1
          //11		1	0	1	0
          //12		1	0	1	1
          //13		1	1	0	0
          //14		1	1	0	1
          //15		1	1	1	0
          //16		1	1	1	1

		  //allow an exception for jldecrement2 if this reaches the end of the sequence
		  if (j-2==ct1->numofbases) jldecrement2=true;

          //case 1, no stacks
          if (jldecrement) {
			  if (ikincrement) {
				  if (en3==w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/+
						w3->f(i+1,k+1)/*w3[i+1][a]*/+constantclosure) {
						stack.push(j-1-ct1->numofbases,l-1-ct2->numofbases/*b*/,1,0,w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/,true);
						stack.push(i+1,k+1/*a*/,ct1->numofbases,0,w3->f(i+1,k+1)/*w3[i+1][a]*/,true);
						found = true;
					}

			  }
			  else if (k==ct2->numofbases&&local) { 
				if (en3==w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/+
					+constantclosure) {
					stack.push(j-1-ct1->numofbases,l-1-ct2->numofbases/*b*/,1,0,w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/,true);
              
					found = true;
				}

            }
			  else if (k==ct2->numofbases) {
				  if (en3==w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/+
						+constantclosure+gap*(ct1->numofbases-i)) {
						stack.push(j-1-ct1->numofbases,l-1-ct2->numofbases/*b*/,1,0,w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/,true);
              
						found = true;
					}

            }

          }

          //case 16, all stack
          if ((i<ct1->numofbases)&&(j>ct1->numofbases+1)&&!found&&jldecrement2&&ikincrement2&&k+2<=ct2->numofbases+1&&j-2-ct1->numofbases>=0&&l-2-ct2->numofbases>=0) { 
            if (en3==w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/+w3->f(i+2,k+2)/*w3[i+2][a]*/+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)
                +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+constantclosure) {
							
              stack.push(j-2-ct1->numofbases,l-2-ct2->numofbases/*b*/,1,0,w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/,true);
              stack.push(i+2,k+2/*a*/,ct1->numofbases,0,w3->f(i+2,k+2)/*w3[i+2][a]*/,true);
              found = true;
							
            }

          }
          //case 6, j and l stacked
          if (j>ct1->numofbases+1&&!found&&jldecrement2&&j-2-ct1->numofbases>=0&&l-2-ct2->numofbases>=0) {
			  if (ikincrement) {
				  if (en3==w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/+w3->f(i+1,k+1)/*w3[i+1][a]*/+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+constantclosure) {

					stack.push(j-2-ct1->numofbases,l-2-ct2->numofbases/*b*/,1,0,w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/,true);
					stack.push(i+1,k+1/*a*/,ct1->numofbases,0,w3->f(i+1,k+1)/*w3[i+1][a]*/,true);
					found = true;
				  }

            }
			  else if (k==ct2->numofbases&&local) {
				  if (en3==w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+constantclosure) {

					stack.push(j-2-ct1->numofbases,l-2-ct2->numofbases/*b*/,1,0,w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/,true);
              
					found = true;
				  }

            }
			  else if (k==ct2->numofbases) {
				  if (en3==w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+constantclosure+gap*(ct1->numofbases-i)) {

					stack.push(j-2-ct1->numofbases,l-2-ct2->numofbases/*b*/,1,0,w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/,true);
              
					found = true;	
				  }

            }
		  }

          //case 11, i and k stacked
          if ((i<ct1->numofbases)&&!found&&jldecrement&&ikincrement2&&k+2<=ct2->numofbases+1) {
            if (en3==w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/
                +w3->f(i+2,k+2)/*w3[i+2][a]*/+edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+constantclosure) {

              stack.push(j-1-ct1->numofbases,l-1-ct2->numofbases/*b*/,1,0,w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/,true);
              stack.push(i+2,k+2/*a*/,ct1->numofbases,0,w3->f(i+2,k+2)/*w3[i+2][a]*/,true);
              found = true;	


            }
          }

		  if (l-2-ct2->numofbases>=0) if (l-2-ct2->numofbases>=lowend[j-1-ct1->numofbases]&&l-2-ct2->numofbases<=highend[j-1-ct1->numofbases]/*lowlimit(j-1-ct1->numofbases,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found) {
            //case 2, l stacked
			  if (ikincrement) {
				  if (en3==w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/+
					w3->f(i+1,k+1)/*w3[i+1][a]*/+edangle5(l,k,l-1,ct2,data)+gap+constantclosure) {
					stack.push(j-1-ct1->numofbases,l-2-ct2->numofbases/*b-1*/,1,0,w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/,true);
					stack.push(i+1,k+1/*a*/,ct1->numofbases,0,w3->f(i+1,k+1)/*w3[i+1][a]*/,true);
					found = true;

				}
			  }
			  else if (k==ct2->numofbases&&local) {
				  if (en3==w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/+
					+edangle5(l,k,l-1,ct2,data)+gap+constantclosure) {
					stack.push(j-1-ct1->numofbases,l-2-ct2->numofbases/*b-1*/,1,0,w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/,true);
              
					found = true;
				  }

			  }
			  else if (k==ct2->numofbases){ 
				if (en3==w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/+
					+edangle5(l,k,l-1,ct2,data)+gap+constantclosure+gap*(ct1->numofbases-i)) {
					stack.push(j-1-ct1->numofbases,l-2-ct2->numofbases/*b-1*/,1,0,w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/,true);
              
					found = true;
				}

			  }

            //case 12, i, k, and l stacked
            if ((i<ct1->numofbases)&&!found&&ikincrement2&&k+2<=ct2->numofbases+1) {
							
              if (en3==w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/+
                  w3->f(i+2,k+2)/*w3[i+2][a]*/+edangle5(l,k,l-1,ct2,data)+edangle3(i,j,i+1,ct1,data)
                  +edangle3(k,l,k+1,ct2,data)+gap+constantclosure) {

                stack.push(j-1-ct1->numofbases,l-2-ct2->numofbases/*b-1*/,1,0,w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/,true);
                stack.push(i+2,k+2/*a*/,ct1->numofbases,0,w3->f(i+2,k+2)/*w3[i+2][a]*/,true);
                found = true;

              }
            }
            if (k+2<=highend[i+1]&&k+2>=lowend[i+1]&&k+2<=ct2->numofbases+1/*highlimit(i+1,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found) {
              //case 4, l and k stacked
              if (en3==w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/+
                  w3->f(i+1,k+2)/*w3[i+1][a+1]*/+edangle5(l,k,l-1,ct2,data)+
                  edangle3(k,l,k+1,ct2,data)+2*gap+constantclosure){

                stack.push(j-1-ct1->numofbases,l-2-ct2->numofbases/*b-1*/,1,0,w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/,true);
                stack.push(i+1,k+2/*a+1*/,ct1->numofbases,0,w3->f(i+1,k+2)/*w3[i+1][a+1]*/,true);
                found = true;

              }
            }

            if (k+1>=lowend[i+2]&&(i<ct1->numofbases)/*lowlimit(i+2,maxsep,ct1->numofbases,ct2->numofbases)/*a>0*/&&!found) {
              //case 10, i and l stacked
				if (k+1<=highend[i+2]) {
					if (en3==w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/+
						w3->f(i+2,k+1)/*w3[i+2][a-1]*/+edangle5(l,k,l-1,ct2,data)+
						edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure){

					  stack.push(j-1-ct1->numofbases,l-2-ct2->numofbases/*b-1*/,1,0,w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/,true);
					  stack.push(i+2,k+1/*a-1*/,ct1->numofbases,0,w3->f(i+2,k+1)/*w3[i+2][a-1]*/,true);
					  found = true;


					}
				}
				else if (k==ct2->numofbases&&local) {
					if (en3==w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/+
						+edangle5(l,k,l-1,ct2,data)+
						edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure){

					  stack.push(j-1-ct1->numofbases,l-2-ct2->numofbases/*b-1*/,1,0,w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/,true);
	                  
					  found = true;


					}
				}
				else if (k==ct2->numofbases) {
					if (en3==w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/+
						+edangle5(l,k,l-1,ct2,data)+
						edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure+gap*(ct1->numofbases-i-1)){

					  stack.push(j-1-ct1->numofbases,l-2-ct2->numofbases/*b-1*/,1,0,w5->f(j-1-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-1-ct1->numofbases][b-1]*/,true);
	                  
					  found = true;


					}
				}
              
            }
									

          }//

          if (k+2<=highend[i+1]&&k+2>=lowend[i+1]&&k+2<=ct2->numofbases+1/*highlimit(i+1,maxsep,ct1->numofbases,ct2->numofbases)/*a<2*maxsep+2*/&&!found) {
            //case 3, k alone stacked
			if (jldecrement) {
				if (en3==w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/+
					w3->f(i+1,k+2)/*w3[i+1][a+1]*/+edangle3(k,l,k+1,ct2,data)+gap+constantclosure){
								
				  stack.push(j-1-ct1->numofbases,l-1-ct2->numofbases/*b*/,1,0,w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/,true);
				  stack.push(i+1,k+2/*a+1*/,ct1->numofbases,0,w3->f(i+1,k+2)/*w3[i+1][a+1]*/,true);
				  found = true;

				}
		    }

            //case 8, j, k, and l stacked
            if ((j>ct1->numofbases+1)&&!found&&jldecrement2&&k+2<=ct2->numofbases+1&&l-2-ct2->numofbases>=0) {
              if (en3==w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/+
                  w3->f(i+1,k+2)/*w3[i+1][a+1]*/+edangle3(k,l,k+1,ct2,data)
                  +edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+gap+constantclosure) {
					
                stack.push(j-2-ct1->numofbases,l-2-ct2->numofbases/*b*/,1,0,w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/,true);
                stack.push(i+1,k+2/*a+1*/,ct1->numofbases,0,w3->f(i+1,k+2)/*w3[i+1][a+1]*/,true);
                found = true;

              }
            }
									
            if (j-2-ct1->numofbases>=0) if (l-1-ct2->numofbases<=highend[j-2-ct1->numofbases]&&l-1-ct2->numofbases>=lowend[j-2-ct1->numofbases]/*highlimit(j-2-ct1->numofbases,maxsep,ct1->numofbases,ct2->numofbases)/*b<2*maxsep+2*/&&!found) {
              //case 7, j and k stacked
              if ((j>ct1->numofbases+1)&&k+2<=ct2->numofbases+1) {
								
                if (en3==w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/+
                    w3->f(i+1,k+2)/*w3[i+1][a+1]*/+edangle3(k,l,k+1,ct2,data)+
                    edangle5(j,i,j-1,ct1,data)+2*gap+constantclosure){

                  stack.push(j-2-ct1->numofbases,l-1-ct2->numofbases/*b+1*/,1,0,w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/,true);
                  stack.push(i+1,k+2/*a+1*/,ct1->numofbases,0,w3->f(i+1,k+2)/*w3[i+1][a+1]*/,true);
                  found = true;	

                }
              }
            }
          }

		  if (j-2-ct1->numofbases>=0) if (l-1-ct2->numofbases<=highend[j-2-ct1->numofbases]&&l-1-ct2->numofbases>=lowend[j-2-ct1->numofbases]/*highlimit(j-2-ct1->numofbases,maxsep,ct1->numofbases,ct2->numofbases)/*b<2*maxsep+2*/&&!found) {
            //case 5, j stacked
            if ((ikincrement)&&(j>ct1->numofbases+1)) {
              if (en3==w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/+
                  w3->f(i+1,k+1)/*w3[i+1][a]*/+edangle5(j,i,j-1,ct1,data)+gap+constantclosure) {

                stack.push(j-2-ct1->numofbases,l-1-ct2->numofbases/*b+1*/,1,0,w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/,true);
                stack.push(i+1,k+1/*a*/,ct1->numofbases,0,w3->f(i+1,k+1)/*w3[i+1][a]*/,true);
                found = true;	

              }
            }
			else if ((k==ct2->numofbases&&local)&&(j>ct1->numofbases+1)) {
              if (en3==w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/+
                  +edangle5(j,i,j-1,ct1,data)+gap+constantclosure) {

                stack.push(j-2-ct1->numofbases,l-1-ct2->numofbases/*b+1*/,1,0,w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/,true);
                
                found = true;	

              }
            }
			else if ((k==ct2->numofbases)&&(j>ct1->numofbases+1)) {
              if (en3==w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/+
                  +edangle5(j,i,j-1,ct1,data)+gap+constantclosure+gap*(ct1->numofbases-i)) {

                stack.push(j-2-ct1->numofbases,l-1-ct2->numofbases/*b+1*/,1,0,w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/,true);
                
                found = true;	

              }
            }

            //case 15, i, j, and k stacked
            if ((i<ct1->numofbases)&&(j>ct1->numofbases+1)&&!found&&ikincrement2&&k+2<=ct2->numofbases+1) {
							
              if (en3==w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/+
                  w3->f(i+2,k+2)/*w3[i+2][a]*/+edangle5(j,i,j-1,ct1,data)
                  +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+gap+constantclosure) {

                stack.push(j-2-ct1->numofbases,l-1-ct2->numofbases/*b+1*/,1,0,w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/,true);
                stack.push(i+2,k+2/*a*/,ct1->numofbases,0,w3->f(i+2,k+2)/*w3[i+2][a]*/,true);
                found = true;


              }
							
            }
            if (k+1>=lowend[i+2]/*lowlimit(i+2,maxsep,ct1->numofbases,ct2->numofbases)/*a>0*/&&!found) {
              //case 13, j and i stacked
              if ((i<ct1->numofbases)&&(j>ct1->numofbases+1)) {
				  if (k+1<=highend[i+2]) {
					  if (en3==w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/+
						w3->f(i+2,k+1)/*w3[i+2][a-1]*/+edangle5(j,i,j-1,ct1,data)+
						edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure){

						stack.push(j-2-ct1->numofbases,l-1-ct2->numofbases/*b+1*/,1,0,w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/,true);
						stack.push(i+2,k+1/*a-1*/,ct1->numofbases,0,w3->f(i+2,k+1)/*w3[i+2][a-1]*/,true);
						found = true;
					  }

                }
				  else if (k==ct2->numofbases&&local) {
					  if (en3==w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/
						+edangle5(j,i,j-1,ct1,data)+
						edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure){

						stack.push(j-2-ct1->numofbases,l-1-ct2->numofbases/*b+1*/,1,0,w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/,true);
                  
						found = true;
					  }

                }
				  else if (k==ct2->numofbases) {
					  if (en3==w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/
						+edangle5(j,i,j-1,ct1,data)+
						edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure+gap*(ct1->numofbases-i-1)){

						stack.push(j-2-ct1->numofbases,l-1-ct2->numofbases/*b+1*/,1,0,w5->f(j-2-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-2-ct1->numofbases][b+1]*/,true);
                  
						found = true;
					  }

                }

              }
            }
          }

          if (k+1>=lowend[i+2]/*lowlimit(i+2,maxsep,ct1->numofbases,ct2->numofbases)/*a>0*/&&!found) {
            //case 9, i stacked
            if ((i<ct1->numofbases)&&jldecrement) {
							
				if (k+1<=highend[i+2]) {
					if (en3==w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/+
						w3->f(i+2,k+1)/*w3[i+2][a-1]*/+edangle3(i,j,i+1,ct1,data)+gap+constantclosure) {
							
						stack.push(j-1-ct1->numofbases,l-1-ct2->numofbases/*b*/,1,0,w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/,true);
						stack.push(i+2,k+1/*a-1*/,ct1->numofbases,0,w3->f(i+2,k+1)/*w3[i+2][a-1]*/,true);
						found = true;

					}
				}
				else if (k==ct2->numofbases&&local) {
					if (en3==w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/
						+edangle3(i,j,i+1,ct1,data)+gap+constantclosure) {
							
						stack.push(j-1-ct1->numofbases,l-1-ct2->numofbases/*b*/,1,0,w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/,true);
                
						found = true;
					}

				}
				else if (k==ct2->numofbases) {
					if (en3==w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/
						+edangle3(i,j,i+1,ct1,data)+gap+constantclosure+gap*(ct1->numofbases-i-1)) {
							
						stack.push(j-1-ct1->numofbases,l-1-ct2->numofbases/*b*/,1,0,w5->f(j-1-ct1->numofbases,l-1-ct2->numofbases)/*w5[j-1-ct1->numofbases][b]*/,true);
                
						found = true;
					}

				}

            }
            //case 14, i, j, and l stacked
			if ((i<ct1->numofbases)&&(j>ct1->numofbases+1)&&!found&&jldecrement2&&j-2-ct1->numofbases>=0&&l-2-ct2->numofbases>=0) {
				if (k+1<=highend[i+2]) { 
					if (en3==w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/+
						w3->f(i+2,k+1)/*w3[i+2][a-1]*/+edangle3(i,j,i+1,ct1,data)
						+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+gap+constantclosure) {

					stack.push(j-2-ct1->numofbases,l-2-ct2->numofbases/*b*/,1,0,w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/,true);
					stack.push(i+2,k+1/*a-1*/,ct1->numofbases,0,w3->f(i+2,k+1)/*w3[i+2][a-1]*/,true);
					found = true;

					}
				}
				else if (k==ct2->numofbases&&local) {
					if (en3==w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/+
						edangle3(i,j,i+1,ct1,data)
						+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+gap+constantclosure) {

						stack.push(j-2-ct1->numofbases,l-2-ct2->numofbases/*b*/,1,0,w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/,true);
                
						found = true;

					}
				}
				else if (k==ct2->numofbases) {
					if (en3==w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/+
						edangle3(i,j,i+1,ct1,data)
						 +edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+gap+constantclosure+gap*(ct1->numofbases-i-1)) {

						stack.push(j-2-ct1->numofbases,l-2-ct2->numofbases/*b*/,1,0,w5->f(j-2-ct1->numofbases,l-2-ct2->numofbases)/*w5[j-2-ct1->numofbases][b]*/,true);
                
						found = true;
					}

				}
            }
          }

		  //reset jl decrement2 for remainder of code... i.e. undo the exception above
		  if (j-2==ct1->numofbases) jldecrement2 = l-2 <= highend[j-2] && l-2 >= lowend[j-2];


		}
        //Now consider multiloop:
        //junction closed by i-j pair aligned with a and b
        //calculate the free energy of 2 fragments merged:
        if (i!=ct1->numofbases&&j!=ct1->numofbases+1&&k!=ct2->numofbases&&l!=ct2->numofbases+1) {
          for (c=i+minloop+1;c<j-minloop&&!found;c++) {
            //if (c>ct1->numofbases) kp=ct1->numofbases-ct2->numofbases;
            //else kp = 0;
            for (d=max(k+minloop+1,lowend[c]/*lowlimit(c,maxsep,ct1->numofbases,ct2->numofbases)/*c-maxsep-kp*/);
                 d<l-minloop&&d<=highend[c]/*highlimit(c,maxsep,ct1->numofbases,ct2->numofbases)/*(c+maxsep-kp*/&&!found;d++) {
              //e = d-c+maxsep+kp;
						
              if ((((c<ct1->numofbases)&&d<ct2->numofbases)||((c>ct1->numofbases)&&(d>ct2->numofbases)))
				  &&(c!=ct1->numofbases)&&(d!=ct2->numofbases)&&d+1>=lowend[c+1]/*lowlimit(c+1,maxsep,ct1->numofbases,ct2->numofbases)*/
				  &&d+1<=highend[c+1]&&d+1>=lowend[c+1]) {
								

                //must consider whether an unpaired nucleotide is stacked
                //is stacked onto each of the nucleotides in a helix
                //There are 16 cases:
                //consider rthem in this order (0 unstacked, 1 stacked)
                //		i	j	k	l
                //1		0	0	0	0
                //2		0	0	0	1
                //3		0	0	1	0
                //4		0	0	1	1
                //5		0	1	0	0
                //6		0	1	0	1
                //7		0	1	1	0
                //8		0	1	1	1
                //9		1	0	0	0
                //10	1	0	0	1
                //11	1	0	1	0
                //12	1	0	1	1
                //13	1	1	0	0
                //14	1	1	0	1
                //15	1	1	1	0
                //16	1	1	1	1

                //case 1 - no stacks
                if (ikincrement&&jldecrement) {
                  if(en3== w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->numofbases)][a][e]*/
                     +w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                     +constantclosure) {
	
                    stack.push(i+1,c,k+1,d,w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->numofbases)][a][e]*/);
                    stack.push(ideref(c+1,j-1,ct1->numofbases)/*ideref(c+1,j-1,ct1->numofbases)*/,jderef(c+1,j-1,ct1->numofbases)/*jderef(c+1,j-1,ct1->numofbases)*/,ideref(d+1,l-1,ct2->numofbases)/*e*/,jderef(d+1,l-1,ct2->numofbases)/*b*/,w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b]*/);
                    found = true;
                  }

                }
									
                //case 16 - all four stacked
                if (!found&&ikincrement2&&jldecrement2) {
                  if((en3 ==w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->numofbases)][a][e]*/
                      +w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                      +4*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
                      +edangle3(k,l,k+1,ct2,data)
                      +edangle5(j,i,j-1,ct1,data)
                      +edangle3(i,j,i+1,ct1,data)+constantclosure)&&
                     ((i+1)!=ct1->numofbases)&&((j-1)!=ct1->numofbases)&&((k+1)!=ct2->numofbases)&&((l-1)!=ct2->numofbases)) {
	
                    found = true;
                    stack.push(i+2,c,k+2/*a*/,d/*e*/,w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->numofbases)][a][e]*/);
                    stack.push(ideref(c+1,j-2,ct1->numofbases),jderef(c+1,j-2,ct1->numofbases),ideref(d+1,l-2,ct2->numofbases)/*e*/,jderef(d+1,l-2,ct2->numofbases)/*b*/,w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b]*/);

                  }
                }
									
                //case 6 - j and l stacked:
                if (ikincrement&&jldecrement2&&!found) {
                  if((en3 == w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->numofbases)][a][e]*/
                      +w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                      +2*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
                      +edangle5(l,k,l-1,ct2,data)+constantclosure)&&((j-1)!=ct1->numofbases)&&((l-1)!=ct2->numofbases)){

                    found = true;
                    stack.push(i+1,c,k+1/*a*/,d/*e*/,w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->numofbases)][a][e]*/);
                    stack.push(ideref(c+1,j-2,ct1->numofbases),jderef(c+1,j-2,ct1->numofbases),ideref(d+1,l-2,ct2->numofbases)/*e*/,jderef(d+1,l-2,ct2->numofbases)/*b*/,w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b]*/);
                  }
                }
									
                //case 11 - i and k stacked
                if (ikincrement2&&jldecrement&&!found) {
                  if((en3 ==w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->numofbases)][a][e]*/
                      +w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                      +2*data->eparam[6]+edangle3(k,l,k+1,ct2,data)+
                      edangle3(i,j,i+1,ct1,data)+constantclosure)&&((i+1)!=ct1->numofbases)&&((k+1)!=ct2->numofbases)){
						
                    found = true;
                    stack.push(i+2,c,k+2/*a*/,d/*e*/,w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->numofbases)][a][e]*/);
                    stack.push(ideref(c+1,j-1,ct1->numofbases),jderef(c+1,j-1,ct1->numofbases),ideref(d+1,l-1,ct2->numofbases)/*e*/,jderef(d+1,l-1,ct2->numofbases)/*b*/,w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b]*/);
                  }
                }
									
                if (l-2>=lowend[j-1]&&l-2<=highend[j-1]/*lowlimit(j-1,maxsep,ct1->numofbases,ct2->numofbases)/*(b-1>=0*/&&!found&&((l-1)!=ct2->numofbases)) {
                  //case 2 - stack on l
                  if (ikincrement) {
                    if(en3== w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->numofbases)][a][e]*/
                       +w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                       +data->eparam[6]+edangle5(l,k,l-1,ct2,data)+gap+constantclosure){

                      found = true;
                      stack.push(i+1,c,k+1/*a*/,d/*e*/,w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->numofbases)][a][e]*/);
                      stack.push(ideref(c+1,j-1,ct1->numofbases),jderef(c+1,j-1,ct1->numofbases),ideref(d+1,l-2,ct2->numofbases)/*e*/,jderef(d+1,l-2,ct2->numofbases)/*b-1*/,w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b-1]*/);

                    }
                  }
									
                  //case 12 - i, k, and l stacked
                  if (!found&&ikincrement2) {
                    if((en3 ==w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->numofbases)][a][e]*/
                        +w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                        +3*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
                        +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+gap
                        +constantclosure)&&((i+1)!=ct1->numofbases)){

		
                      found = true;
                      stack.push(i+2,c,k+2/*a*/,d/*e*/,w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->numofbases)][a][e]*/);
                      stack.push(ideref(c+1,j-1,ct1->numofbases),jderef(c+1,j-1,ct1->numofbases),ideref(d+1,l-2,ct2->numofbases)/*e*/,jderef(d+1,l-2,ct2->numofbases)/*b-1*/,w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b-1]*/);


                    }
                  }
									
                  if (k+1>=lowend[i+2]&&k+1<=highend[i+2]/*lowlimit(i+2,maxsep,ct1->numofbases,ct2->numofbases)/*a-1>0*/&&!found&&(i+1)!=ct1->numofbases) {
                    //case 10 - l and i stacked
                    if(en3 ==w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->numofbases)][a-1][e]*/
                       +w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                       +2*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
                       +edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure){


                      found = true;
                      stack.push(i+2,c,k+1/*a-1*/,d/*e*/,w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->numofbases)][a-1][e]*/);
                      stack.push(ideref(c+1,j-1,ct1->numofbases),jderef(c+1,j-1,ct1->numofbases),ideref(d+1,l-2,ct2->numofbases)/*e*/,jderef(d+1,l-2,ct2->numofbases)/*b-1*/,w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b-1]*/);


                    }
									
                  }

                  if (k+2<=highend[i+1]&&k+2>=lowend[i+1]/*highlimit(i+1,maxsep,ct1->numofbases,ct2->numofbases)/*a+1<2*maxsep+2*/&&!found&&((k+1)!=ct2->numofbases)) {
                    //case 4 - k and l stacked
                    if(en3 ==w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->numofbases)][a+1][e]*/
                       +w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                       +2*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
                       +edangle3(k,l,k+1,ct2,data)+2*gap+constantclosure){

                      found = true;
                      stack.push(i+1,c,k+2/*a+1*/,d/*e*/,w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->numofbases)][a+1][e]*/);
                      stack.push(ideref(c+1,j-1,ct1->numofbases),jderef(c+1,j-1,ct1->numofbases),ideref(d+1,l-2,ct2->numofbases)/*e*/,jderef(d+1,l-2,ct2->numofbases)/*b-1*/,w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b-1]*/);


                    }

											
											
										
                  }
										

										


                }
                if (k+1>=lowend[i+2]&&k+1<=highend[i+2]/*lowlimit(i+2,maxsep,ct1->numofbases,ct2->numofbases)/*a-1>0*/&&!found&&(i+1!=ct1->numofbases)) {
                  //case 14 - i, j, and l stacked
                  if (jldecrement2) {
                    if((en3 ==w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->numofbases)][a-1][e]*/
                        +w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                        +3*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
                        +edangle5(j,i,j-1,ct1,data)
                        +edangle3(i,j,i+1,ct1,data)+gap+constantclosure)&&(j-1!=ct1->numofbases)){


                      found = true;
                      stack.push(i+2,c,k+1/*a-1*/,d/*e*/,w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->numofbases)][a-1][e]*/);
                      stack.push(ideref(c+1,j-2,ct1->numofbases),jderef(c+1,j-2,ct1->numofbases),ideref(d+1,l-2,ct2->numofbases)/*e*/,jderef(d+1,l-2,ct2->numofbases)/*b*/,w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b]*/);

                    }
                  }
									
                  //case 9 - i stacked
                  if (jldecrement&&!found) {
                    if(en3 ==w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->numofbases)][a-1][e]*/
                       +w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                       +data->eparam[6]+edangle3(i,j,i+1,ct1,data)+gap+constantclosure){

                      found = true;
                      stack.push(i+2,c,k+1/*a-1*/,d/*e*/,w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->numofbases)][a-1][e]*/);
                      stack.push(ideref(c+1,j-1,ct1->numofbases),jderef(c+1,j-1,ct1->numofbases),ideref(d+1,l-1,ct2->numofbases)/*e*/,jderef(d+1,l-1,ct2->numofbases)/*b*/,w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b]*/);


                    }
                  }
									
                  if (!found&&l-1<=highend[j-2]&&l-1>=lowend[j-2]/*highlimit(j-2,maxsep,ct1->numofbases,ct2->numofbases)/*(b+1<2*maxsep+2)*/&&(j-1!=ct1->numofbases)) {
                    //case 13 - i and j stacked
                    if(en3 ==w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->numofbases)][a-1][e]*/
                       +w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                       +2*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
                       +edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure){


                      found = true;
                      stack.push(i+2,c,k+1/*a-1*/,d/*e*/,w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->numofbases)][a-1][e]*/);
                      stack.push(ideref(c+1,j-2,ct1->numofbases),jderef(c+1,j-2,ct1->numofbases),ideref(d+1,l-1,ct2->numofbases)/*e*/,jderef(d+1,l-1,ct2->numofbases)/*b+1*/,w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b+1]*/);

                    }
													
                  }

                }
									
                if (k+2<=highend[i+1]&&k+2>=lowend[i+1]/*highlimit(i+1,maxsep,ct1->numofbases,ct2->numofbases)/*a+1<2*maxsep+2*/&&!found&&(k+1!=ct2->numofbases)) {
                  //case 8 - j, k, and l stacked
                  if (jldecrement2) {
                    if((en3==w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->numofbases)][a+1][e]*/
                        +w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                        +3*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
                        +edangle5(l,k,l-1,ct2,data)+edangle3(k,l,k+1,ct2,data)+gap
                        +constantclosure)&&(j-1!=ct1->numofbases)){

                      found = true;
                      stack.push(i+1,c,k+2/*a+1*/,d/*e*/,w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->numofbases)][a+1][e]*/);
                      stack.push(ideref(c+1,j-2,ct1->numofbases),jderef(c+1,j-2,ct1->numofbases),ideref(d+1,l-2,ct2->numofbases)/*e*/,jderef(d+1,l-2,ct2->numofbases)/*b*/,w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b]*/);

                    }
                  }
												
                  //case 3 - stack on k 
                  if (jldecrement&&!found) {
                    if(en3 ==w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->numofbases)][a+1][e]*/+2*data->eparam[10]
                       +w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b]*/+2*data->eparam[5]
                       +data->eparam[6]+edangle3(k,l,k+1,ct2,data)+gap+constantclosure){

                      found = true;
                      stack.push(i+1,c,k+2/*a+1*/,d/*e*/,w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->numofbases)][a+1][e]*/);
                      stack.push(ideref(c+1,j-1,ct1->numofbases),jderef(c+1,j-1,ct1->numofbases),ideref(d+1,l-1,ct2->numofbases)/*e*/,jderef(d+1,l-1,ct2->numofbases)/*b*/,w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->numofbases)][iref(c+1,j-1,ct1->numofbases)][e][b]*/);

                    }
                  }
							
                  if (!found&&l-1<=highend[j-2]&&l-1>=lowend[j-2]/*highlimit(j-2,maxsep,ct1->numofbases,ct2->numofbases)/*b+1<2*maxsep+2*/) {
                    //case 7 - j and k stacked
                    if ((en3 ==w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->numofbases)][a+1][e]*/
                         +w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                         +2*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
                         +edangle3(k,l,k+1,ct2,data)+2*gap+constantclosure)&&(j-1!=ct1->numofbases)){

                      found = true;
                      stack.push(i+1,c,k+2/*a+1*/,d/*e*/,w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->numofbases)][a+1][e]*/);
                      stack.push(ideref(c+1,j-2,ct1->numofbases),jderef(c+1,j-2,ct1->numofbases),ideref(d+1,l-1,ct2->numofbases)/*e*/,jderef(d+1,l-1,ct2->numofbases)/*b+1*/,w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b+1]*/);

                    }
																				
                  }	
										

                }
                if (l-1<=highend[j-2]&&l-1>=lowend[j-2]/*highlimit(j-2,maxsep,ct1->numofbases,ct2->numofbases)/*b+1<2*maxsep+2*/&&!found&&(j-1!=ct1->numofbases)) {
                  //case 15 - i, j, and k stacked
                  if (ikincrement2) {
                    if((en3 ==w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->numofbases)][a][e]*/
                        +w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                        +3*data->eparam[6]+edangle3(k,l,k+1,ct2,data)
                        +edangle5(j,i,j-1,ct1,data)
                        +edangle3(i,j,i+1,ct1,data)+gap+constantclosure)&&(i+1!=ct1->numofbases)&&(k+1!=ct2->numofbases)){

	
                      found = true;
                      stack.push(i+2,c,k+2/*a*/,d/*e*/,w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->numofbases)][a][e]*/);
                      stack.push(ideref(c+1,j-2,ct1->numofbases),jderef(c+1,j-2,ct1->numofbases),ideref(d+1,l-1,ct2->numofbases)/*e*/,jderef(d+1,l-1,ct2->numofbases)/*b+1*/,w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b+1]*/);

                    }
                  }
								
                  //case 5 - j stacked
                  if (!found&&ikincrement) {
                    if(en3 ==w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->numofbases)][a][e]*/
                       +w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                       +data->eparam[6]+edangle5(j,i,j-1,ct1,data)+gap+constantclosure){

                      found = true;
                      stack.push(i+1,c,k+1/*a*/,d/*e*/,w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->numofbases)][a][e]*/);
                      stack.push(ideref(c+1,j-2,ct1->numofbases),jderef(c+1,j-2,ct1->numofbases),ideref(d+1,l-1,ct2->numofbases)/*e*/,jderef(d+1,l-1,ct2->numofbases)/*b+1*/,w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->numofbases)][iref(c+1,j-2,ct1->numofbases)][e][b+1]*/);

                    }
                  }
                }

              }			

									
							

            }
          }
        }

        if (!found) {
			//a traceback error occurred
          //cerr << "Traceback error!!\n";
			error = 14;
          return error;
        }

	
      }

			
			

    }
    else {//v[i][j][a][b]!=en3 so we need to search for w solutions
      found = false;

      //check if it is safe to increment i and k
      ikincrement =
        k+1 >= lowend[i+1]/*lowlimit(i+1, maxsep, ct1->numofbases, ct2->numofbases)*/ &&
        k+1 <= highend[i+1]/*highlimit(i+1, maxsep, ct1->numofbases, ct2->numofbases)*/;


      //check to see if it is safe to decrement j and l
      jldecrement =
        l-1 <= highend[j-1]/*highlimit(j-1, maxsep, ct1->numofbases, ct2->numofbases)*/ &&
        l-1 >= lowend[j-1]/*lowlimit(j-1, maxsep, ct1->numofbases, ct2->numofbases)*/;

			
      //case 6
      if (jldecrement) {
        if (en3==w->f(i,j-1,k,l-1)/*w[j-1][iref(i,j-1,ct1->numofbases)][a][b]*/+2*data->eparam[6]&&(j-1!=ct1->numofbases)&&(l-1!=ct2->numofbases)){
          found = true;
          stack.push(i,j-1,k/*a*/,l-1,w->f(i,j-1,k,l-1)/*w[j-1][iref(i,j-1,ct1->numofbases)][a][b]*/);
        }
      }
      //case 11	
      if (!found&&ikincrement){
        if (en3==w->f(i+1,j,k+1,l)/*w[j][iref(i+1,j,ct1->numofbases)][a][b]*/+2*data->eparam[6]&&i!=ct1->numofbases&&k!=ct2->numofbases) {

          found = true;
          stack.push(i+1,j,k+1,l,w->f(i+1,j,k+1,l));
        }
      }
      //case 16
      if (!found&&ikincrement&&jldecrement) {
        if (j-1>i+1&&i!=ct1->numofbases&&k!=ct2->numofbases) {
				
          if (en3 ==w->f(i+1,j-1,k+1,l-1)/*w[j-1][iref(i+1,j-1,ct1->numofbases)][a][b]*/+ 4*data->eparam[6]) {
            found = true;
            stack.push(i+1,j-1,k+1,l-1,w->f(i+1,j-1,k+1,l-1));
          }
        }
      }

						

      if (/*a>=1*/k>=lowend[i+1]&&k<=highend[i+1]/*lowlimit(i+1,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found) {
        //case 9
        if(en3==w->f(i+1,j,k,l)/*w[j][iref(i+1,j,ct1->numofbases)][a-1][b]*/+data->eparam[6]+gap&&(i!=ct1->numofbases)){
          found = true;
          stack.push(i+1,j,k,l,w->f(i+1,j,k,l));

        }

        //case 14
        else if((j-1>i+1)&&(j-1!=ct1->numofbases)&&jldecrement) {
					
          if(en3==w->f(i+1,j-1,k,l-1)/*w[j-1][iref(i+1,j-1,ct1->numofbases)][a-1][b]*/+3*data->eparam[6]+gap){
            found = true;
            stack.push(i+1,j-1,k,l-1,w->f(i+1,j-1,k,l-1));
          }

        }

        if (/*b+1<(2*maxsep+2)*/l<=highend[j-1]&&l>=lowend[j-1]/*highlimit(j-1,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found&&(j-1>i+1)&&(j-1!=ct1->numofbases)) {
          //case 13
          if(en3==w->f(i+1,j-1,k,l)/*w[j-1][iref(i+1,j-1,ct1->numofbases)][a-1][b+1]*/+2*data->eparam[6]+2*gap){
            found = true;
            stack.push(i+1,j-1,k,l,w->f(i+1,j-1,k,l));

          }

        }
        if (/*b>=1*/l-1>=lowend[j]/*lowlimit(j,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found&&(i!=ct1->numofbases)&&(l-1!=ct2->numofbases)) {
          //case 10
          if(en3==w->f(i+1,j,k,l-1)/*w[j][iref(i+1,j,ct1->numofbases)][a-1][b-1]*/+2*data->eparam[6]+2*gap){
            found = true;
            stack.push(i+1,j,k,l-1,w->f(i+1,j,k,l-1));

          }

        }
      }
      if (/*b+1<(2*maxsep+2)*/l<=highend[j-1]&&l>=lowend[j-1]/*highlimit(j-1,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found&&(j-1!=ct1->numofbases)) {
							
        //case 5
        if(en3==w->f(i,j-1,k,l)/*w[j-1][iref(i,j-1,ct1->numofbases)][a][b+1]*/+data->eparam[6]+gap){
          found = true;
          stack.push(i,j-1,k,l,w->f(i,j-1,k,l));

        }

				

        else if (k+1<=highend[i]/*highlimit(i,maxsep,ct1->numofbases,ct2->numofbases)/*a+1<(2*maxsep+2)*/&&(k!=ct2->numofbases)&&(j-1!=ct1->numofbases)) {
          //case 7
          if(en3==w->f(i,j-1,k+1,l)/*w[j-1][iref(i,j-1,ct1->numofbases)][a+1][b+1]*/+2*data->eparam[6]+2*gap){
            found = true;
            stack.push(i,j-1,k+1,l,w->f(i,j-1,k+1,l));

          }

        }

        //case 15
        else if((j-1>i+1)&&(k!=ct2->numofbases)&&(j-1!=ct1->numofbases)&&ikincrement) {
          if(en3==w->f(i+1,j-1,k+1,l)/*w[j-1][iref(i+1,j-1,ct1->numofbases)][a][b+1]*/+3*data->eparam[6]+gap) {
            found = true;
            stack.push(i+1,j-1,k+1,l,w->f(i+1,j-1,k+1,l));
          }

        }

      }
      if (/*a+1<(2*maxsep+2)*/k+1<=highend[i]/*highlimit(i,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found&&(k!=ct2->numofbases)) {
        //case 3
        if(en3==w->f(i,j,k+1,l)/*w[j][iref(i,j,ct1->numofbases)][a+1][b]*/+data->eparam[6]+gap){
          found = true;
          stack.push(i,j,k+1,l,w->f(i,j,k+1,l));
		
        }

        //case 8
        else if (jldecrement){
          if (en3==w->f(i,j-1,k+1,l-1)/*w[j-1][iref(i,j-1,ct1->numofbases)][a+1][b]*/+3*data->eparam[6]+gap&&(l-1!=ct2->numofbases)) {
            found = true;
            stack.push(i,j-1,k+1,l-1,w->f(i,j-1,k+1,l-1));
          }
        }
							
        if (!found&&/*b>=1*/l-1>=lowend[j]/*lowlimit(j,maxsep,ct1->numofbases,ct2->numofbases)*/&&(l-1!=ct2->numofbases)) {
          //case 4
          if(en3==w->f(i,j,k+1,l-1)/*w[j][iref(i,j,ct1->numofbases)][a+1][b-1]*/+2*data->eparam[6]+2*gap){
            found = true;
            stack.push(i,j,k+1,l-1,w->f(i,j,k+1,l-1));

          }
        }
      }
      if (/*b>=1*/l-1>=lowend[j]/*lowlimit(j,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found&&(l-1!=ct2->numofbases)) {
        //case 2
        if(en3==w->f(i,j,k,l-1)/*w[j][iref(i,j,ct1->numofbases)][a][b-1]*/+data->eparam[6]+gap) {
          found = true;
          stack.push(i,j,k,l-1,w->f(i,j,k,l-1));
        }

        //case12
        else if (ikincrement){
          if ((en3==w->f(i+1,j,k+1,l-1)/*w[j][iref(i+1,j,ct1->numofbases)][a][b-1]*/+3*data->eparam[6]+gap&&i!=ct1->numofbases)&&(k!=ct2->numofbases)){
            found = true;
            stack.push(i+1,j,k+1,l-1,w->f(i+1,j,k+1,l-1));
          }
        }
      }
      //Consider the case where none of the four nucs (i,j,k,l) are paired
      //Consider whether any of the 4 nucleotides are stacked on a helix
      //There are 16 cases:
      //consider them in this order (0 unstacked, 1 stacked)
      //		i	j	k	l
      //1		0	0	0	0
      //2		0	0	0	1
      //3		0	0	1	0
      //4		0	0	1	1
      //5		0	1	0	0
      //6		0	1	0	1
      //7		0	1	1	0
      //8		0	1	1	1
      //9		1	0	0	0
      //10	1	0	0	1
      //11	1	0	1	0
      //12	1	0	1	1
      //13	1	1	0	0
      //14	1	1	0	1
      //15	1	1	1	0
      //16	1	1	1	1
			
				
      //case 1 - nothing stacked:
      if(en3==v->f(i,j,k,l)/*v[j][iref(i,j,ct1->numofbases)][a][b]*/+2*data->eparam[10]+penalty(i,j,ct1,data)+penalty(k,l,ct2,data)&&!found) {
        found = true;
        stack.push(i,j,k,l,v->f(i,j,k,l));

      }

						
						
      //case 6 - j and l stacked
      if (!found&&jldecrement) {
        if(en3==v->f(i,j-1,k,l-1)/*v[j-1][iref(i,j-1,ct1->numofbases)][a][b]*/+2*data->eparam[6]+edangle3(j-1,i,j,ct1,data)
           +edangle3(l-1,k,l,ct2,data)+2*data->eparam[10]+penalty(i,j-1,ct1,data)+penalty(l-1,k,ct2,data)
           &&(j-1!=ct1->numofbases)&&(l!=ct2->numofbases+1)) {

          found = true;
          stack.push(i,j-1,k,l-1,v->f(i,j-1,k,l-1));
        }
      }

      //case 11 - i and k stacked
      if(!found&&ikincrement) {
        if ((i!=ct1->numofbases)&&(k!=ct2->numofbases)&&en3==v->f(i+1,j,k+1,l)/*v[j][iref(i+1,j,ct1->numofbases)][a][b]*/+2*data->eparam[6]+edangle5(k+1,l,k,ct2,data)+
            edangle5(i+1,j,i,ct1,data)+2*data->eparam[10]+penalty(i+1,j,ct1,data)+penalty(k+1,l,ct2,data)) {

          found = true;
          stack.push(i+1,j,k+1,l,v->f(i+1,j,k+1,l));
        }

      }

      //case 16 - i, j, k, and l stacked
      if(!found&&(j-1>i+1)&&(i!=ct1->numofbases)&&(k!=ct2->numofbases)&&(j-1!=ct1->numofbases)&&(l!=ct2->numofbases+1)&&ikincrement&&jldecrement) {
        if (en3==v->f(i+1,j-1,k+1,l-1)/*v[j-1][iref(i+1,j-1,ct1->numofbases)][a][b]*/+4*data->eparam[6]+edangle3(l-1,k+1,l,ct2,data)
            +edangle5(k+1,l-1,k,ct2,data)
            +edangle3(j-1,i+1,j,ct1,data)
            +edangle5(i+1,j-1,i,ct1,data)+2*data->eparam[10]+
            penalty(i+1,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)) {


          found = true;
          stack.push(i+1,j-1,k+1,l-1,v->f(i+1,j-1,k+1,l-1));

        }


      }



      if (!found&&/*b-1>=0*/l-1>=lowend[j]/*lowlimit(j,maxsep,ct1->numofbases,ct2->numofbases)*/&&(l!=ct2->numofbases+1)) {
        //case 2 - l stacked
        if(en3==v->f(i,j,k,l-1)/*v[j][iref(i,j,ct1->numofbases)][a][b-1]*/+data->eparam[6]+edangle3(l-1,k,l,ct2,data)+2*data->eparam[10]+gap+
           penalty(i,j,ct1,data)+penalty(l-1,k,ct2,data)) {


          found = true;
          stack.push(i,j,k,l-1,v->f(i,j,k,l-1));

        }
						
        else if (/*(a+1<2*maxsep+2)*/k+1<=highend[i]/*highlimit(i,maxsep,ct1->numofbases,ct2->numofbases)*/&&(k!=ct2->numofbases)) {
          //case 4 - l and k stacked
          if(en3==v->f(i,j,k+1,l-1)/*v[j][iref(i,j,ct1->numofbases)][a+1][b-1]*/+2*data->eparam[6]+edangle3(l-1,k+1,l,ct2,data)
             +edangle5(k+1,l-1,k,ct2,data)+2*data->eparam[10]+2*gap
             +penalty(i,j,ct1,data)+penalty(l-1,k+1,ct2,data)) {

            found = true;
            stack.push(i,j,k+1,l-1,v->f(i,j,k+1,l-1));

          }
        }

						

						
      }

      if (!found&&k+1<=highend[i]/*highlimit(i,maxsep,ct1->numofbases,ct2->numofbases)/*a+1<2*maxsep+2*/&&(k!=ct2->numofbases)) {

        //case 8 - j, k, and l stacked:
        if (jldecrement) {
          if(en3 ==v->f(i,j-1,k+1,l-1)/*v[j-1][iref(i,j-1,ct1->numofbases)][a+1][b]*/+3*data->eparam[6]+edangle3(j-1,i,j,ct1,data)
             +edangle3(l-1,k+1,l,ct2,data)+edangle5(k+1,l-1,k,ct2,data)+2*data->eparam[10]+gap
             +penalty(i,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)&&(l!=ct2->numofbases+1)&&(j-1!=ct1->numofbases)) {


            found = true;
            stack.push(i,j-1,k+1,l-1,v->f(i,j-1,k+1,l-1));
          }
        }

        //case 3 - k stacked
        if(!found&&en3==v->f(i,j,k+1,l)/*v[j][iref(i,j,ct1->numofbases)][a+1][b]*/+data->eparam[6]+edangle5(k+1,l,k,ct2,data)+2*data->eparam[10]+gap+
           penalty(i,j,ct1,data)+penalty(k+1,l,ct2,data)) {


          found = true;
          stack.push(i,j,k+1,l,v->f(i,j,k+1,l));

        }

        if (!found&&l<=highend[j-1]&&l>=lowend[j-1]/*highlimit(j-1,maxsep,ct1->numofbases,ct2->numofbases)/*b+1<2*maxsep+2*/&&(j!=ct1->numofbases+1)) {
          //case 7 - j and k stacked
          if(en3==v->f(i,j-1,k+1,l)/*v[j-1][iref(i,j-1,ct1->numofbases)][a+1][b+1]*/+2*data->eparam[6]+edangle3(j-1,i,j,ct1,data)
             +edangle5(k+1,l,k,ct2,data)+2*data->eparam[10]+2*gap+
             penalty(i,j-1,ct1,data)+penalty(k+1,l,ct2,data)) {


            found = true;
            stack.push(i,j-1,k+1,l,v->f(i,j-1,k+1,l));

          }

        }
							

      }

      if (l<=highend[j-1]&&l>=lowend[j-1]/*highlimit(j-1,maxsep,ct1->numofbases,ct2->numofbases)/*b+1<2*maxsep+2*/&&!found&&(j!=ct1->numofbases+1)) {
				
        //case 5 - j stacked
        if(en3==v->f(i,j-1,k,l)/*v[j-1][iref(i,j-1,ct1->numofbases)][a][b+1]*/+data->eparam[6]+edangle3(j-1,i,j,ct1,data)+2*data->eparam[10]+gap+
           penalty(i,j-1,ct1,data)+penalty(k,l,ct2,data)) {


          found = true;
          stack.push(i,j-1,k,l,v->f(i,j-1,k,l));

        }
					
        //case 15 - i,j, and k stacked
        else if (ikincrement) {
          if(en3==v->f(i+1,j-1,k+1,l)/*v[j-1][iref(i+1,j-1,ct1->numofbases)][a][b+1]*/+3*data->eparam[6]+edangle5(k+1,l,k,ct2,data)
             +edangle3(j-1,i+1,j,ct1,data)
             +edangle5(i+1,j-1,i,ct1,data)+2*data->eparam[10]+gap+
             penalty(i+1,j-1,ct1,data)+penalty(k+1,l,ct2,data)&&(i!=ct1->numofbases)) {

            found = true;
            stack.push(i+1,j-1,k+1,l,v->f(i+1,j-1,k+1,l));
          }

        }

        if (!found&&/*a>=1*/k>=lowend[i+1]&&k<=highend[i+1]/*lowlimit(i+1,maxsep,ct1->numofbases,ct2->numofbases)*/&&(j-1>i+1)&&(i!=ct1->numofbases)) {
          //case 13 - i and j stacked

          if(en3==v->f(i+1,j-1,k,l)/*v[j-1][iref(i+1,j-1,ct1->numofbases)][a-1][b+1]*/+2*data->eparam[6]+edangle3(j-1,i+1,j,ct1,data)
             +edangle5(i+1,j-1,i,ct1,data)+2*data->eparam[10]+2*gap+
             penalty(i+1,j-1,ct1,data)+penalty(k,l,ct2,data)) {


            found = true;
            stack.push(i+1,j-1,k,l,v->f(i+1,j-1,k,l));

          }
        }
						
      }

      if (!found&&k>=lowend[i+1]&&k<=highend[i+1]/*lowlimit(i+1,maxsep,ct1->numofbases,ct2->numofbases)/*a>=1*/&&(i!=ct1->numofbases)) {
        //case 9 - i alone is stacked
        if(en3==v->f(i+1,j,k,l)/*v[j][iref(i+1,j,ct1->numofbases)][a-1][b]*/+data->eparam[6]+edangle5(i+1,j,i,ct1,data)+2*data->eparam[10]+gap+
           penalty(i+1,j,ct1,data)+penalty(k,l,ct2,data)) {
          found = true;
          stack.push(i+1,j,k,l,v->f(i+1,j,k,l));


        }


        //case 14 - i, j, and l stacked
        else if(j-1>i+1&&(j!=ct1->numofbases+1)&&(l!=ct2->numofbases+1)&&jldecrement) {
          if(en3==v->f(i+1,j-1,k,l-1)/*v[j-1][iref(i+1,j-1,ct1->numofbases)][a-1][b]*/+3*data->eparam[6]+edangle3(l-1,k,l,ct2,data)
             +edangle3(j-1,i+1,j,ct1,data)
             +edangle5(i+1,j-1,i,ct1,data)+2*data->eparam[10]+gap+
             penalty(i+1,j-1,ct1,data)+penalty(k,l-1,ct2,data)) {
	
            found = true;
            stack.push(i+1,j-1,k,l-1,v->f(i+1,j-1,k,l-1));

          }


        }
      }


      if (!found&&l-1>=lowend[j]/*lowlimit(j,maxsep,ct1->numofbases,ct2->numofbases)/*b-1>=0*/&&(i!=ct1->numofbases)&&(l!=ct2->numofbases+1)) {
        //case 12 - i, k, and l stacked:
        if (ikincrement) {
          if(en3==v->f(i+1,j,k+1,l-1)/*v[j][iref(i+1,j,ct1->numofbases)][a][b-1]*/+3*data->eparam[6]+edangle3(l-1,k+1,l,ct2,data)
             +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+2*data->eparam[10]+gap+
             penalty(i+1,j,ct1,data)+penalty(l-1,k+1,ct2,data)&&(k!=ct2->numofbases)) {


            found = true;
            stack.push(i+1,j,k+1,l-1,v->f(i+1,j,k+1,l-1));

          }
        }
        if (!found&&k>=lowend[i+1]&&k<=highend[i+1]/*lowlimit(i+1,maxsep,ct1->numofbases,ct2->numofbases)/*a>=1*/) {
          //case 10 - l and i stacked
          if(en3==v->f(i+1,j,k,l-1)/*v[j][iref(i+1,j,ct1->numofbases)][a-1][b-1]*/+2*data->eparam[6]+edangle3(l-1,k,l,ct2,data)
             +edangle5(i+1,j,i,ct1,data)+2*data->eparam[10]+2*gap+
             penalty(i+1,j,ct1,data)+penalty(k,l-1,ct2,data)) {

            found = true;
            stack.push(i+1,j,k,l-1,v->f(i+1,j,k,l-1));

          }
        }
      }
						
      //calculate the free energy of 2 fragments merged:
      for (c=i+minloop;c<j-minloop&&!found;c++) {
        //if (c>ct1->numofbases) kp = ct1->numofbases-ct2->numofbases;
        //else kp = 0;
        for (d=max(k+minloop,/*c-maxsep-kp*/lowend[c]/*lowlimit(c,maxsep,ct1->numofbases,ct2->numofbases)*/);
             d<l-minloop&&d<=/*c+maxsep-kp*/highend[c]/*highlimit(c,maxsep,ct1->numofbases,ct2->numofbases)*/&&!found;d++) {
				
          //e = d-c+maxsep+kp;
          if ((c!=ct1->numofbases)&&(d!=ct2->numofbases)&&d+1>=lowend[c+1]&&d+1<=highend[c+1]) {			
            if (c<ct1->numofbases) {
              if (d<ct2->numofbases) {
				  if (d+1>=lowend[c+1]&&d+1<=highend[c+1]) {
					if(en3 ==w->f(i,c,k,d)/*w[c][i][a][e]*/+w->f(c+1,j,d+1,l)/*w[j][iref(c+1,j,ct1->numofbases)][e][b]*/) {
						found = true;
						stack.push(i,c,k,d,w->f(i,c,k,d));
						stack.push(c+1,j,d+1,l,w->f(c+1,j,d+1,l));

					}
				  }

              }
            }
            else {
              if (d>ct2->numofbases) {
                if(en3 ==w->f(i,c,k,d)/*w[c][iref(i,c,ct1->numofbases)][a][e]*/+
                   w->f(c+1-ct1->numofbases,j-ct1->numofbases,d+1-ct2->numofbases,l-ct2->numofbases)
                   /*w[j-ct1->numofbases][iref(c+1-ct1->numofbases,j-ct1->numofbases,ct1->numofbases)][e][b]*/) {
                  found = true;
                  stack.push(i,c,k,d,w->f(i,c,k,d));
                  stack.push(c+1-ct1->numofbases,j-ct1->numofbases,d+1-ct2->numofbases,l-ct2->numofbases,w->f(c+1-ct1->numofbases,j-ct1->numofbases,d+1-ct2->numofbases,l-ct2->numofbases));

                }
              }
            }
          }
        }
      }

      if (!found) {
		  //a traceback error occurred
        //cerr << "Traceback error!\n"; 
		  error = 14;
        return error;
      }
    }
  }
  //end of while (stack.pull(&i,&j, &k, &l /*&a,&b*/, &en3, &open))
  return error;
}


//Perform multiple trackbacks.
//Return an int that indicates whether an error occurred (0=no error).
int  dyntraceback(short maxtracebacks, short window, short awindow,
                  short percentsort,
                  varray *v, dynalignarray *w, wendarray *w3, wendarray *w5,
                  /*bool **pair,*/
                  structure *ct1, structure *ct2, short **alignment,
                  short *lowend, short *highend, short int gapincrease,
                  datatable *data, bool singleinsert,
                  integersize lowest, dynalignarray *vmod,
                  bool local) {

  short int i,j,k,l,gap,N,N2,index;
  int ip,jp,kp,cp;
  integersize en1;
  integersize crit;
  dynalignheap heap;
  int point;
  bool **mark1, **mark2, **marka;
  bool modification;
  int error = 0;
  int localerror;


  if (ct1->nmod>0||ct2->nmod>0) modification = true;
  else modification = false;

  short int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
                       {0,1,0,1,0,0},{0,0,0,0,0,0}};

  //declarations for traceback portion:
	
  if (lowest>=0) {
    //return with empty CT files
    ct1->numofstructures=1;
    ct2->numofstructures=1;
    //initialize the structures:
    for (index=1;index<=ct1->numofbases;index++) ct1->basepr[ct1->numofstructures][index] = 0;
    for (index=1;index<=ct2->numofbases;index++) ct2->basepr[ct2->numofstructures][index] = 0;

    //initialize the alignment:
    for (index=0;index<=ct1->numofbases;index++) alignment[ct1->numofstructures-1][index]=0;
    return 0;

  }

  //maxsep = maxseparation;
  gap = gapincrease;
  N = ct1->numofbases;
  N2 = ct2->numofbases;


  mark1 = new bool *[N+1];
  for (i=0;i<=N;i++) mark1[i] = new bool [i];
  mark2 = new bool *[N2+1];
  for (i=0;i<=N2;i++) mark2[i] = new bool [i];
  
  marka = new bool *[N+1];
  for (i=0;i<=N;i++) marka[i] = new bool [highend[i]-lowend[i]+1/*2*maxsep+2*/];

  for (i=0;i<=N;i++) {
    for (j=0;j<i;j++) {
      mark1[i][j] = false;
    }
  }
  for (i=0;i<=N2;i++) {
    for (j=0;j<i;j++) {
      mark2[i][j] = false;
    }
  }
  for (i=0;i<=N;i++) {
    for (j=0;j<highend[i]-lowend[i]+1/*2*maxsep+2*/;j++) {
      marka[i][j] = false;
    }
	//change the pointer for marka for faster indexing:
	marka[i]-=lowend[i];
  }

  



  //now traceback:

  //For suboptimal structure prediction, build a heap of possible traceback starts
  crit = DYNALIGN_INFINITY;
  for (i=1;i<=ct1->numofbases;i++) {
    for (/*a=0;a<2*maxsep+2;a++*/k=max(lowend[i]/*lowlimit(i,maxsep,ct1->numofbases,ct2->numofbases)*/,1);k<=highend[i]/*highlimit(i,maxsep,ct1->numofbases,ct2->numofbases)/*&&k<=N2*/;k++) {
		
      //k = a+N-maxsep;
      if (k<=N2) {
        if ((w5->f(i,k) /*w5[N][a]*/ +gap*abs(N-i-(N2-k)) /*N2 - N - a + maxsep*/ )<crit) {
          crit = w5->f(i,k)/*w5[N][a]*/+gap*abs(N-i-(N2-k))/*N2 - N - a + maxsep*/;	
        }
      }
		
    }
  }
  if (percentsort>0) crit = crit - (integersize)((float) crit*(((float) percentsort)/100.0 ));
  else crit = crit - percentsort;


  for (i=1;i<=N;i++) {
    for (j=i+minloop;j<=N;j++) {
      for (k=max(/*i-maxsep*/lowend[i]/*lowlimit(i,maxsep,ct1->numofbases,ct2->numofbases)*/,1);k<=(min(ct2->numofbases,/*i+maxsep*/highend[i]/*highlimit(i,maxsep,ct1->numofbases,ct2->numofbases)*/));k++) {
        for (l=max(/*j-maxsep*/lowend[j]/*lowlimit(j,maxsep,ct1->numofbases,ct2->numofbases)*/,k);l<=(min(ct2->numofbases,/*j+maxsep*/highend[j]/*highlimit(j,maxsep,ct1->numofbases,ct2->numofbases)*/));l++) {
					
          //use a and b when refering to the energy arrays
          //a = k-i+maxsep;
          //b = l-j+maxsep;	

          //if (/*(a+N2-N>=0)&&(a+N2-N<=2*maxsep)*/) {
          if (modification) {
            if (vmod->f(i,j,k,l)/*vmod[j][i][a][b]*/+vmod->f(j,i+N,l,k+N2)/*vmod[i+N][j-i][b][a]*/<=crit) 
              heap.push(i,j,k,l,vmod->f(i,j,k,l)+vmod->f(j,i+N,l,k+N2));	
          }
          else if (v->f(i,j,k,l)/*v[j][i][a][b]*/+v->f(j,i+N,l,k+N2)/*v[i+N][j-i][b][a]*/<=crit) 
            heap.push(i,j,k,l,v->f(i,j,k,l)+v->f(j,i+N,l,k+N2));
          //}

        }
      }
    }

  }

  //make a heap out of dynalign heap:
  for (ip=1;ip<heap.size;ip++) {
    jp = ip;
    kp = jp/2;
    while ((heap.peak(jp)<heap.peak(kp))&&kp>=0){
      heap.swap(jp,kp);
      jp = jp/2;
      kp = kp/2;
    }

  }
  //sort the dynalign heap with heap sort
  for (ip = heap.size-2;ip>=0;ip--) {
    heap.swap(ip+1,0);

    jp = 0;
    cp = 1;
    while (cp<=ip) {
      if (cp!=ip) {
        if (heap.peak(cp+1)<heap.peak(cp)) cp++;

      }
      if (heap.peak(cp)<heap.peak(jp)) {
        heap.swap(cp,jp);
        jp = cp;
        cp = 2*cp;
      }
      else cp = ip + 1;
    }

  }


  //call dyntrace to determine suboptimal structures:

  point = heap.size;
  ct1->numofstructures = 0;
  ct2->numofstructures = 0;
  while (point>0&&ct1->numofstructures<maxtracebacks) {
		
    point--;
		

    heap.read(point,&i,&j,&k,&l,&en1);	
    //k = a+i-maxsep;
    //l = b+j-maxsep;

    if (!mark1[j][i]||!mark2[l][k]||!marka[i][k/*reference(i, k, N,N2,maxsep)/*a*/]||!marka[j][l/*reference(j, l, N,N2,maxsep)/*b*/]) {
      ct1->numofstructures++;
      ct2->numofstructures++;
      ct1->checknumberofstructures();
      ct2->checknumberofstructures();
	  if(ct1->numofstructures!=1) {
		  strcpy(ct1->ctlabel[ct1->numofstructures],ct1->ctlabel[1]);
		
		  strcpy(ct2->ctlabel[ct1->numofstructures],ct2->ctlabel[1]);
	  }
      //initialize the structures:
      for (index=1;index<=ct1->numofbases;index++) ct1->basepr[ct1->numofstructures][index] = 0;
      for (index=1;index<=ct2->numofbases;index++) ct2->basepr[ct2->numofstructures][index] = 0;

      //initialize the alignment:
      for (index=0;index<=ct1->numofbases;index++) alignment[ct1->numofstructures-1][index]=0;

      //dyntrace(57,126,61,131,ct1,ct2,ct1->numofstructures,alignment[ct1->numofstructures-1],w,v,w3,w5,maxsep,data,gap,vmod);
      localerror = dyntrace(i,j,k,l,
               ct1,ct2,ct1->numofstructures,
               alignment[ct1->numofstructures-1],
               w,v,w3,w5,lowend,highend,data,gap,vmod,
               local);
	  if (localerror!=0) error = localerror;
      localerror = dyntrace(j,i+N,l,k+N2,
               ct1,ct2,ct1->numofstructures,
               alignment[ct1->numofstructures-1],
               w,v,w3,w5,lowend,highend,data,gap,vmod,
               local);
	  if (localerror!=0) error = localerror;
      ct1->energy[ct1->numofstructures] = en1;
      ct2->energy[ct1->numofstructures] = en1;

			
			
			
      for (i=1;i<N;i++) {
        if (ct1->basepr[ct1->numofstructures][i]>i) {
          if (!mark1[ct1->basepr[ct1->numofstructures][i]][i]) {
						
            for (j=i-window;j<=i+window;j++) {
              for (k=ct1->basepr[ct1->numofstructures][i]-window;
                   k<=ct1->basepr[ct1->numofstructures][i]+window;k++) {

                if (j<k&&j>0&&k<=ct1->numofbases) {
                  mark1[k][j] = true;
                }
              }
            }
          }
        }
      }

      for (i=1;i<N2;i++) {
        if (ct2->basepr[ct2->numofstructures][i]>i) {
          if (!mark2[ct2->basepr[ct2->numofstructures][i]][i]) {
						
            for (j=i-window;j<=i+window;j++) {
              for (k=ct2->basepr[ct2->numofstructures][i]-window;
                   k<=ct2->basepr[ct2->numofstructures][i]+window;k++) {

                if (j<k&&j>0&&k<=ct2->numofbases) {
                  mark2[k][j] = true;
                }
              }
            }
          }
        }
      }

      for (i=1;i<=N;i++) {
        if (alignment[ct1->numofstructures-1][i]) {
          if (!marka[i][alignment[ct1->numofstructures-1][i]]) {/*[reference(i, alignment[ct1->numofstructures-1][i], N,N2,maxsep)]/*alignment[ct1->numofstructures-1][i]-i+maxsep]*/
							
            for (j=i-awindow;j<=i+awindow;j++) {
			  for (k=alignment[ct1->numofstructures-1][i]-awindow;k<=alignment[ct1->numofstructures-1][i]+awindow;k++) {
			  //for (k=reference(i, alignment[ct1->numofstructures-1][i], N,N2,maxsep)/*alignment[ct1->numofstructures-1][i]-i+maxsep*/-awindow;
              //     k<=reference(i, alignment[ct1->numofstructures-1][i], N,N2,maxsep)/*alignment[ct1->numofstructures-1][i]-i+maxsep*/+awindow;
              //     k++) {
									
                if (j>0&&j<=ct1->numofbases&&k>=lowend[j]&&k<=highend[j]) {
                  marka[j][k] = true;
                }
              }
            }
          }
        }
      }
    }
  }

  for (i=0;i<=N;i++) delete[] mark1[i];
  delete[] mark1;
  for (i=0;i<=N2;i++) delete[] mark2[i];;
  delete[] mark2;
  for (i=0;i<=N;i++) {
	  //move back the pointer:
	  marka[i]+=lowend[i];
	  delete[] marka[i];
  }
  delete[] marka;

  //return the error code
  return error;
}

void alignout(short** align, const char *aout, structure *ct1, structure *ct2) {
  ofstream out;
  short i,j,last,next,index;
  char *line1,*line2,*line3;
	
  line1 = new char [ct1->numofbases+ct2->numofbases+100];
  line2 = new char [ct1->numofbases+ct2->numofbases+100];
  line3 = new char [ct1->numofbases+ct2->numofbases+100];

  out.open(aout);

  for (index=0;index<ct1->numofstructures;index++) {
    strcpy(line1,"");
    strcpy(line2,"");
    strcpy(line3,"");

		

    last = 0;
    for (i=1;i<=ct1->numofbases;i++) {
      if (last==ct2->numofbases) {

        //nothing more can be put down in line 2 for the second sequence

        line1[strlen(line1)+1]='\0';
        line1[strlen(line1)]=ct1->nucs[i];
        strcat(line2,"-");
        strcat(line3," ");

      }
      else if (align[index][i]>0) {
        //this nucleotide is aligned to something
        while (align[index][i]!=last+1) {
          //need to lay nucleotides down in second sequence
          strcat(line1,"-");
          last++;
          line2[strlen(line2)+1]='\0';
          line2[strlen(line2)]=ct2->nucs[last];
          strcat(line3," ");
			
        }
        //now put i in alignment with last+1
        line1[strlen(line1)+1]='\0';
        line1[strlen(line1)]=ct1->nucs[i];
        last++;
        line2[strlen(line2)+1]='\0';
        line2[strlen(line2)]=ct2->nucs[last];
        strcat(line3,"^");

      }
      else {//align[i]=0
        next=0;
        for (j=i+1;j<=ct1->numofbases&&next==0;j++) next = align[index][j]; 
        if (next==last+1) {
          //no nucs from seq2 can be put down next to seq1
          line1[strlen(line1)+1]='\0';
          line1[strlen(line1)]=ct1->nucs[i];
          strcat(line2,"-");
          strcat(line3," ");

        }
        else {
          line1[strlen(line1)+1]='\0';
          line1[strlen(line1)]=ct1->nucs[i];
          last++;
          line2[strlen(line2)+1]='\0';
          line2[strlen(line2)]=ct2->nucs[last];
          strcat(line3," ");

        }


      }


    }

    //put down any nucs from seq2 that did not go down:
    for (i=last+1;i<=ct2->numofbases;i++) {
      strcat(line1,"-");
      line2[strlen(line2)+1]='\0';
      line2[strlen(line2)]=ct2->nucs[i];
      strcat(line3," ");
    }

    out<<"Alignment #"<<index+1<<" Score= "<<ct1->energy[index+1]<<"\n";
    out <<line1<<"\n"<<line2<<"\n"<<line3<<"\n\n\n";
  }

  out.close();
  delete[] line1;
  delete[] line2;
  delete[] line3;
}

void parse(structure *ct, char *seq1, char *seq2) {
  short int *seqone, *seqtwo;
  short int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
                       {0,1,0,1,0,0},{0,0,0,0,0,0}};
  short int i,j;
  char temp[5];

  seqone = new short int [ct->numofbases+1];
  seqtwo = new short int [ct->numofbases+1];

  j = 1;
  for (i=0;i<(short) (strlen(seq1));i++) {
    strncpy(temp,seq1+i,1);
    strcpy(temp+1,"\0");
    if (strcmp(temp,"-")) {
      seqone[j] = tonumi(temp);
      strncpy(temp,seq2+i,1);
      strcpy(temp+1,"\0");
      seqtwo[j] = tonumi(temp);
      j++;
    }
  }

  for (i=1;i<=ct->numofbases;i++) {
    for (j=i+1;j<=ct->numofbases;j++) {
      if (inc[seqone[i]][seqone[j]]&&inc[seqtwo[i]][seqtwo[j]]) {
        ct->tem[j][i] = true;
      }
else {
        ct->tem[j][i] = false;
      }
    }
  }
  delete[] seqone;
  delete[] seqtwo;
}



void opendynalignsavefile(const char *filename, structure *ct1, structure *ct2,
                        varray *v, dynalignarray *w, dynalignarray *vmod,
                          wendarray *w3, wendarray *w5, datatable *data,
                          bool *singleinsert, short *maxsep, short *gap,
                          short *lowest, bool *local, bool **allowed_alignments, 
						  short *lowend, short *highend) {
  short i,j,k,l,a,b,c,d,I;
  short int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
                       {0,1,0,1,0,0},{0,0,0,0,0,0}};
  int temp;
  bool optimalonly,modification;

  ifstream sav(filename,ios::binary);

  //start with structure information
  read(&sav,&temp);//read the modification/optimalonly key - this has been noted previously
  //because it is needed for array allocation
  //The flag returns:
  //0 - suboptimal, no modifications
  //1 - suboptimal, with modifications
  //2 - optimal only, no modifications
  //3 - optimal only, with modifications
  if (temp==1||temp==3) modification = true;
  else modification = false;
  read(&sav,&(ct1->numofbases));
  read(&sav,&(ct2->numofbases));
  read(&sav,maxsep);
  ct1->allocate(ct1->numofbases);
  read(&sav,&(ct1->npair));
  for (i=0;i<=ct1->npair;i++) {
    read(&sav,&(ct1->pair[i][0]));
    read(&sav,&(ct1->pair[i][1]));
  }
  for (i=0;i<=ct1->numofbases;i++) {
		
    read(&sav,&(ct1->hnumber[i]));
    sav.read(&(ct1->nucs[i]),1);

  }

  for (i=0;i<=2*ct1->numofbases;i++) read(&sav,&(ct1->numseq[i]));
	
  read(&sav,&(ct1->ndbl));
  for (i=0;i<=ct1->ndbl;i++) read(&sav,&(ct1->dbl[i]));

  read(&sav,&(ct1->intermolecular));
  if (ct1->intermolecular) {
    for (i=0;i<3;i++) read(&sav,&(ct1->inter[i]));
		
  }

  read(&sav,&(ct1->nnopair));
  for (i=0;i<=ct1->nnopair;i++) read(&sav,&(ct1->nopair[i]));

  read(&sav,&(ct1->nmod));
  for (i=0;i<=ct1->nmod;i++) read(&sav,&(ct1->mod[i]));
	
  read(&sav,&(ct1->ngu));
  for (i=0;i<=ct1->ngu;i++) read(&sav,&(ct1->gu[i]));  
	
  read(&sav,ct1->ctlabel[1]);



  read(&sav,&(ct1->templated));
  if (ct1->templated) {	
    ct1->allocatetem();
    for (i=0;i<=ct1->numofbases;i++) {
      for (j=0;j<=i;j++) {
        read(&sav,&(ct1->tem[i][j]));	
      }
    }
  }
    
  ct2->allocate(ct2->numofbases);
  read(&sav,&(ct2->npair));
  for (i=0;i<=ct2->npair;i++) {
    read(&sav,&(ct2->pair[i][0]));
    read(&sav,&(ct2->pair[i][1]));
  }
  for (i=0;i<=ct2->numofbases;i++) {
    read(&sav,&(ct2->hnumber[i]));
    sav.read(&(ct2->nucs[i]),1);
  }

  for (i=0;i<=2*ct2->numofbases;i++) read(&sav,&(ct2->numseq[i]));
	
  read(&sav,&(ct2->ndbl));
  for (i=0;i<=ct2->ndbl;i++) read(&sav,&(ct2->dbl[i]));

  read(&sav,&(ct2->intermolecular));
  if (ct2->intermolecular) {
    for (i=0;i<3;i++) read(&sav,&(ct2->inter[i]));
		
  }

  read(&sav,&(ct2->nnopair));
  for (i=0;i<=ct2->nnopair;i++) read(&sav,&(ct2->nopair[i]));

  read(&sav,&(ct2->nmod));
  for (i=0;i<=ct2->nmod;i++) read(&sav,&(ct2->mod[i]));
	
  read(&sav,&(ct2->ngu));
  for (i=0;i<=ct2->ngu;i++) read(&sav,&(ct2->gu[i]));  
	
  read(&sav,ct2->ctlabel[1]);

  read(&sav,&(ct2->templated));
	
  if (ct2->templated) {
    ct2->allocatetem();
    for (i=0;i<=ct2->numofbases;i++) {
      for (j=0;j<=i;j++) read(&sav,&(ct2->tem[i][j]));	

    }

  }




  //now read the thermodynamic data:
  for (i=0;i<5;i++) read(&sav,&(data->poppen[i]));
  read(&sav,&(data->maxpen));
  for (i=0;i<11;i++) read(&sav,&(data->eparam[i]));
  for (i=0;i<31;i++) {
    read(&sav,&(data->inter[i]));
    read(&sav,&(data->bulge[i]));
    read(&sav,&(data->hairpin[i]));

  }

  for (i=0;i<6;i++) {
    for (j=0;j<6;j++) {
      for (k=0;k<6;k++) {
        for (l=0;l<3;l++) {
          read(&sav,&(data->dangle[i][j][k][l]));	
        }
        for (l=0;l<6;l++) {
          read(&sav,&(data->stack[i][j][k][l]));
          read(&sav,&(data->tstkh[i][j][k][l]));
          read(&sav,&(data->tstki[i][j][k][l]));
          read(&sav,&(data->coax[i][j][k][l]));
          read(&sav,&(data->tstackcoax[i][j][k][l]));
          read(&sav,&(data->coaxstack[i][j][k][l]));
          read(&sav,&(data->tstack[i][j][k][l]));
          read(&sav,&(data->tstkm[i][j][k][l]));
          read(&sav,&(data->tstki23[i][j][k][l]));
          read(&sav,&(data->tstki1n[i][j][k][l]));
          for (a=0;a<6;a++) {
            for (b=0;b<6;b++) {
              read(&sav,&(data->iloop11[i][j][k][l][a][b]));
              for (c=0;c<6;c++) {
                if (inc[i][j]&&inc[b][c]) read(&sav,&(data->iloop21[i][j][k][l][a][b][c]));
                for (d=0;d<6;d++) {
                  if (inc[i][k]&&inc[j][l])
                    read(&sav,&(data->iloop22[i][j][k][l][a][b][c][d]));
                }
              }
            }
          }
        }
      }
    }
  }
  read(&sav,&(data->numoftloops));
  for (i=0;i<=data->numoftloops;i++) {
    for (j=0;j<2;j++) read(&sav,&(data->tloop[i][j]));

  }
  read(&sav,&(data->numoftriloops));
  for (i=0;i<=data->numoftriloops;i++) {
    for (j=0;j<2;j++) read(&sav,&(data->triloop[i][j]));

  }
  read(&sav,&(data->numofhexaloops));
  for (i=0;i<=data->numofhexaloops;i++) {
    for (j=0;j<2;j++) read(&sav,&(data->hexaloop[i][j]));

  }
  read(&sav,&(data->auend));
  read(&sav,&(data->gubonus));
  read(&sav,&(data->cint));
  read(&sav,&(data->cslope));
  read(&sav,&(data->c3));
  read(&sav,&(data->efn2a));
  read(&sav,&(data->efn2b));
  read(&sav,&(data->efn2c));
  read(&sav,&(data->init));
  read(&sav,&(data->mlasym));
  read(&sav,&(data->strain));
  read(&sav,&(data->prelog));
  read(&sav,&(data->singlecbulge));
	
  //now read the folding data:
  optimalonly = temp==2 || temp==3;
 

  read(&sav,gap);
  read(&sav,lowest);
  read(&sav,singleinsert);

	
	//read allowed_alignments if maxsep < 0
	if (*maxsep < 0) {

		for (i=0;i<=ct1->numofbases;++i) {
			for (j=0;j<=ct2->numofbases;++j) {
				read(&sav,&(allowed_alignments[i][j]));
			}
		}
	}

	//now calculate the lowend and highend arrays that indicate the first and last index to which a nucleotide
		//from sequence 1 can be aligned.

	if (*maxsep>0) {
		//Using the traditional M parameter to constrain the alignment space
		for (i=0;i<2*ct1->numofbases;++i) {
			lowend[i]=lowlimit(i,*maxsep,ct1->numofbases,ct2->numofbases);
			highend[i]=highlimit(i,*maxsep,ct1->numofbases,ct2->numofbases);
		}
	}
	else {
		//allowed_alignments must be defined to constrain the alignment space
		//need to determine lowend and highend
		for (i=0;i<2*ct1->numofbases;++i) {
			lowend[i]=lowlimit(i,allowed_alignments,ct1->numofbases,ct2->numofbases);
			highend[i]=highlimit(i,allowed_alignments,ct1->numofbases,ct2->numofbases);
		}


	}

	//Allocate v, w, vmod, w3, w5
	v->allocate(ct1->numofbases,ct2->numofbases,lowend,highend,ct1->tem,optimalonly);
	w->allocate(ct1->numofbases,ct2->numofbases,lowend,highend,optimalonly);
	if (modification) vmod->allocate(ct1->numofbases,ct2->numofbases,lowend,highend,optimalonly);

	if (!optimalonly) w3->allocate(ct1->numofbases, ct2->numofbases, lowend,highend);
	w5->allocate(ct1->numofbases, ct2->numofbases, lowend,highend);


  for (i=0;i<=ct1->numofbases;++i) {

    if (temp==2||temp==3) I = ct1->numofbases;
    else I = i+ct1->numofbases-1;
					
    for (j=i;j<=I;++j) {
				
				
      for (k=lowend[i];k<=highend[i];k++) {
					
					
        for (l=lowend[j];l<=highend[j];l++) {
          if (j>ct1->numofbases) {
            b = i;
            a = j-ct1->numofbases;
          }
          else {
            b = j;
            a = i;
          }

          if (ct1->tem[b][a]) read(&sav,&(v->f(i,j,k,l)));
          read(&sav,&(w->f(i,j,k,l)));
          
          if (modification) read(&sav,&(vmod->f(i,j,k,l)));
          
						
        }
      }

    }
  }
			


		


  for (i=0;i<=ct1->numofbases;i++) {
    for (j=lowend[i];j<=highend[i];j++) {
      if (temp < 2) read(&sav,&(w3->f(i,j))); //do not attemptto read w3 if this is an optimal-only file
      read(&sav,&(w5->f(i,j)));

    }

  }
  if (temp < 2) for (j=0/*lowlimit(ct1->numofbases+1,*maxsep,ct1->numofbases,ct2->numofbases)*/;j<=ct2->numofbases+1/*highlimit(ct1->numofbases+1,*maxsep,ct1->numofbases,ct2->numofbases)*/;j++) read(&sav,&(w3->f(ct1->numofbases+1,j)));

  read(&sav,&temp);
  if (temp==0) *local=false;
  else *local=true;
	
  sav.close();
}


//Refold dynalign using save file information
//return an int that indicates whther an error occurred (0 = no error)
int refolddynalign(const char* filename, structure *ct1, structure *ct2,
                    short **alignment, short maxtracebacks,
                    short window, short awindow, short percentsort) {
  short i,ip,kp,ipp,kpp,cp;
  bool singleinsert;
	
  short maxsep,gap,lowest,index,*lowend,*highend;
  //integersize ****v,****w,****vmod,**w3,**w5;
  dynalignarray *w,*vmod;
  varray *v;
  wendarray *w3,*w5;
  datatable *data;
  short int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
                       {0,1,0,1,0,0},{0,0,0,0,0,0}};
  int modificationflag;
  bool optimalonly;
  bool local,found;
  bool **allowed_alignments;
  int error;
  

  

  data = new datatable();

  //open the save file to peek at the sizes needed to allocate arrays:
  ifstream sav(filename,ios::binary);

  read(&sav, &modificationflag); //note that this includes modification and
  //optimal only data

  if (modificationflag==2||modificationflag==3) optimalonly=true;
  else optimalonly=false;

  //start with structure information
  read(&sav,&(ct1->numofbases));
  read(&sav,&(ct2->numofbases));
  read(&sav,&maxsep);
  sav.close();

	
  if (maxsep<0) {
	//This means that the allowed alignments were constrained by the HMM alignment
	  //probs and not by the traditional M parameter.
	  //Therefore, space must be allocated for the allowed_alignment arrays.

	  allowed_alignments=new bool *[ct1->numofbases+1];
	  for (i=0;i<=ct1->numofbases;i++) allowed_alignments[i]=new bool [ct2->numofbases+1];

  }
  else allowed_alignments=NULL;
  //fill the low and highend arrays:
  //allocate space:
  lowend = new short [2*ct1->numofbases];
  highend = new short [2*ct1->numofbases];
	
  //DHM (11/30/06): move these allocations to where the file is read
  if (modificationflag==1||modificationflag==3) {

    vmod = new dynalignarray();


  }
  else vmod=NULL;

  v = new varray();
  w = new dynalignarray();


	

  if (!optimalonly) w3 = new wendarray();
  else w3= NULL;
  w5 = new wendarray();

	
	
  opendynalignsavefile(filename, ct1, ct2, v, w, vmod, w3, w5, data, &singleinsert,  &maxsep, &gap, &lowest,&local,allowed_alignments,lowend,highend);

  //double up the sequences as part of suboptimal traceback support:
  for (i=1;i<=ct1->numofbases;i++) {
    ct1->numseq[(ct1->numofbases)+i] = ct1->numseq[i];
  }
  for (i=1;i<ct2->numofbases;i++) ct2->numseq[ct2->numofbases+i]=ct2->numseq[i];

  if (!optimalonly) {
    error = dyntraceback(maxtracebacks,window,awindow,percentsort,
                 v,w,w3,w5,
                 ct1,ct2,alignment,lowend,highend,gap,data,singleinsert,lowest,vmod,local);
  } else {

    //need to preset the structure and alignment info to zeros
    for (index=1;index<=ct1->numofbases;index++) ct1->basepr[1][index]=0;
    for (index=1;index<=ct2->numofbases;index++) ct2->basepr[1][index]=0;

    for (index=0;index<=ct1->numofbases;index++) alignment[0][index]=0;

    ct1->numofstructures=1;
    ct2->numofstructures=1;

    //calculate the total energy 
    ct1->energy[1]=lowest;
    ct2->energy[1]=lowest;
    found = false;
    for (ip=1;ip<=ct1->numofbases&&!found;++ip) {

								
		
      cp = min(ct2->numofbases,highlimit(ip,maxsep,ct1->numofbases,ct2->numofbases));
      for (kp=max(1,lowlimit(ip,maxsep,ct1->numofbases,ct2->numofbases));kp<=cp&&!found;++kp) {
        if (local) {
          if (lowest==w5->f(ip,kp)) {
            ipp=ip;
            kpp=kp;
            found=true;
          }

        }
        else if (abs((ct1->numofbases-ip)-(ct2->numofbases-kp))*gap+w5->f(ip,kp)==lowest) {
          ipp=ip;
          kpp=kp;
          found=true;
        }
      }
    }

    error = dyntrace(1, ipp, 1, kpp, ct1, ct2, 1, alignment[0],
             w, v, w3, w5, lowend, highend,data, gap, vmod, local, true);
  }

	

  delete v;
  delete w;
  delete w3;
  delete w5;

  delete lowend;
  delete highend;

  if (maxsep<0) {
	  
		for (i=0;i<=ct1->numofbases;i++) delete[] allowed_alignments[i];
		delete[] allowed_alignments;
  }

  if (modificationflag==1||modificationflag==3) delete vmod;
		

	
  //delete[] pair[0];
  //delete[] pair[1];
  //delete[] pair;
  delete data;

  return error;
}




//Fold a single sequence to decide what pairs should be allowed in a subsequent dynalign calculation.
void templatefromfold(structure *ct, datatable *data, int singlefold_subopt_percent) {
  integersize *w5,*w3;
  bool *lfce,*mod;
  short crit,i,j;
  int vmin;
  arrayclass *w2,*wmb2;

  //allocate everything
  arrayclass w(ct->numofbases);
  arrayclass v(ct->numofbases);
  arrayclass wmb(ct->numofbases);	
  forceclass fce(ct->numofbases);


  lfce = new bool [2*ct->numofbases+1];
  mod = new bool [2*ct->numofbases+1];

  for (i=0;i<=2*ct->numofbases;i++) {
    lfce[i] = false;
    mod[i] = false;
  }
	
  w5 = new integersize [ct->numofbases+1];
  w3 = new integersize [ct->numofbases+2];
  

  for (i=0;i<=ct->numofbases;i++) {
    w5[i] = 0;
    w3[i] = 0;
    
  }
  w3[ct->numofbases+1] = 0;

  if (ct->intermolecular) {
    w2 = new arrayclass(ct->numofbases);
    wmb2 = new arrayclass(ct->numofbases);
		
    
		
  }
  else {
    w2 = NULL;
    wmb2 = NULL;
    
  }

  force(ct,&fce,lfce);
  vmin=DYNALIGN_INFINITY;

  //perform the fill steps:(i.e. fill arrays v and w.)
  fill(ct, v, w, wmb, fce, vmin,lfce, mod,w5, w3, false, data, w2, wmb2);

  //find the dots here:
	

	crit= ((short) (abs(vmin)*(float (singlefold_subopt_percent)/100.0)));

  crit = crit + vmin;

	
  i = 1;
  j = 2;
  while (i<(ct->numofbases)) {

    if ((v.f(i,j)+v.f(j,i+ct->numofbases))>crit) {
      //don't allow this pair
      ct->tem[j][i]=false;
   			
    }
    j++;
    if (j>ct->numofbases) {
      i++;
      j=i+1;
    }
  }
	
  delete[] lfce;
  delete[] mod;

  delete[] w5;
  delete[] w3;

  if (ct->intermolecular) {
    delete w2;
    delete wmb2;

    
  }
}


//load a dynalign saved file and determine all possible pairs
//in structures within percentdots of the lowest free energy structure.  These pairs are then
//the only allowed pairs for dyalign structure prediction for the passed structure.
int templatefromdsv(structure *cttemplate, const char *savefilename, float maxdsvchange, int maxpairs) 
{
	
	datatable *data;

	short i,j,k,l,crit;
	bool singleinsert,local;
	short maxsep,gap,lowest,*lowend,*highend;
	
	dynalignarray *w,*vmod;
	varray *v;
	wendarray *w3,*w5;

	bool optimalonly;
	
	short int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};
	int modificationflag;

	dynalignheap heap;

	bool **allowed_alignments;
	

	structure *ct1, *ct2;
	ct1 = new structure();
	ct2 = new structure();
		
	
	data = new datatable();

	

	//open the save file to peek at the sizes needed to allocate arrays:
	ifstream sav(savefilename,ios::binary);


	//read a flag that indicates whether chemical modification constraints were used
	read(&sav, &modificationflag);

	if (modificationflag==2||modificationflag==3) optimalonly=true;
	else optimalonly=false;

	//start with structure information needed for array allocation
	read(&sav,&(ct1->numofbases));
	read(&sav,&(ct2->numofbases));
	read(&sav,&maxsep);
	sav.close();

	if (maxsep<0) {
	//This means that the allowed alignments were constrained by the HMM alignment
	  //probs and not by the traditional M parameter.
	  //Therefore, space must be allocated for the allowed_alignment arrays.

	  allowed_alignments=new bool *[ct1->numofbases+1];
	  for (i=0;i<=ct1->numofbases;i++) allowed_alignments[i]=new bool [ct2->numofbases+1];

	}
	else allowed_alignments=NULL;
	
	//allocate the arrays
	if (modificationflag==1||modificationflag==3) {
		vmod = new dynalignarray();
	} else {
		vmod = NULL;
	}
  
	v = new varray();
	w = new dynalignarray();

	if (!optimalonly) w3 = new wendarray();
	else w3= NULL;
	w5 = new wendarray();


	lowend = new short [2*ct1->numofbases];
	highend = new short [2*ct1->numofbases];

	//open the save file
	// this function does not work on iced for some reason
	opendynalignsavefile(savefilename, ct1, ct2, v, w, vmod, w3, w5, data, &singleinsert,  &maxsep, &gap, &lowest,&local,allowed_alignments,lowend,highend);

	//double up the sequences as part of suboptimal traceback support:
	for (i=1;i<=ct1->numofbases;i++) {
		ct1->numseq[(ct1->numofbases)+i] = ct1->numseq[i];
	}
	for (i=1;i<ct2->numofbases;i++) ct2->numseq[ct2->numofbases+i]=ct2->numseq[i];
	
	

	crit = DYNALIGN_INFINITY;
	for (i=1;i<=ct1->numofbases;i++) {
		for (k=max(lowend[i],1);k<=min(highend[i],ct2->numofbases);k++) {
		
			
			
				if (local) { 
					if ((w5->f(i,k))<crit) {
						crit = w5->f(i,k);
					}
				}
				else { 
					if ((w5->f(i,k)+gap*abs(ct1->numofbases-i-(ct2->numofbases-k))  )<crit) {
						crit = w5->f(i,k)+gap*abs(ct1->numofbases-i-(ct2->numofbases-k));
					}
				}

		}
	}

//zane 07/20/2009: 
	int minimum = crit; // lowest free energy * 10, usually a negative number
	const float stepsize = 0.0002; //i.e. 0.02%, it is the stepsize of energy increment in percentage 
	unsigned int step = 1; // the number of steps
	int dsv_energy = (1-stepsize*step)*minimum; //It's gonna be the energy threshold defined by maxpairs
	unsigned int PairsCount=0; // count the number of pairs, with each of which the lowest energy structure has energy lower than dsv_energy
	while ( PairsCount < maxpairs ){
	  PairsCount=0;
	  for (i=1;i<=ct1->numofbases;i++)
	    for (j=i+minloop;j<=ct1->numofbases;j++) {
	      bool lower = false;
	      
	      for (k=max(lowend[i],1);k<=(min(ct2->numofbases,highend[i]));k++) {
		for (l=max(lowend[j],k);l<=(min(ct2->numofbases,highend[j]));l++)
		  if (v->f(i,j,k,l)+v->f(j,i+ct1->numofbases,l,k+ct2->numofbases)<dsv_energy){
		    lower = true;//the lowest energy structure( with i paired with j, k paired with l, i and j aligned to k and l respectively) is found to be lower than dsv_energy
		    break;//terminate l loop when this lowest energy structure is lower than dsv_energy
		  }		
		if(lower)
		  break; //terminate k loop when this lowest energy structure is lower than dsv_energy
	      }
	      
	      if(lower) // count this i-j pair in PairsCount
		++PairsCount;
	    }
	  //cout << "step" << step <<" : " << PairsCount << endl;
	  ++step;
	  dsv_energy = (1-stepsize*step)*minimum;
	  if (dsv_energy > 0)//if the maxpairs is so large that dsv_energy is increased to zero
	    break; //report the error code
	}
if ( dsv_energy > 0) maxdsvchange =99;
else {
	--step;
   float change = step * stepsize * 100;
   maxdsvchange = maxdsvchange < change ? change : maxdsvchange;
  } cout << "set actual maxdsvchange: " << maxdsvchange << endl;


	crit = (short) ((float) crit -  ((float) crit*(((float) maxdsvchange)/100.0 )));

	////
	for (i=1;i<=ct1->numofbases;i++) {
		for (j=i+minloop;j<=ct1->numofbases;j++) {

			cttemplate->tem[j][i]=false; //the default is to not allow the pair...
			for (k=max(lowend[i],1);k<=(min(ct2->numofbases,highend[i]));k++) {
				for (l=max(lowend[j],k);l<=(min(ct2->numofbases,highend[j]));l++) {
						if (v->f(i,j,k,l)+v->f(j,i+ct1->numofbases,l,k+ct2->numofbases)<=crit) 
							cttemplate->tem[j][i]=true;
				}
			}
		}
	}

	delete w3;
	delete w5;
	delete w;
	delete v;

	if (modificationflag==1||modificationflag==3) delete vmod;

	if (maxsep<0) {
	  
		for (i=0;i<=ct1->numofbases;i++) delete[] allowed_alignments[i];
		delete[] allowed_alignments;
	}

	delete data;
	delete ct1;
	delete ct2;
	return 0;
}

//Using a set of known pairs for ct that have been previously loaded, set those as the only allowed pairs for a Dynalign calculation.
void templatefromct(structure *ct)
{
	short i,j;
	
	for (i=1;i<=ct->numofbases;i++) {
		for (j=i+minloop;j<=ct->numofbases;j++) {

			if (ct->basepr[1][i]==j) ct->tem[j][i]=true; //allow this pair because it occurs in the structure 
			else ct->tem[j][i]=false; //the default is to not allow the pair...
			
		}
	}

	
}


//Read alignment constraints from disk
void readalignmentconstraints(const char *filename, short **forcealign, structure *ct1, structure *ct2) {
	ifstream in;
	int i,j;
	
	in.open(filename);
	

	in >> i;
	in >> j;
	while (i!=-1) {
		forcealign[0][i]=j;
		forcealign[1][j]=i;

	        in >> i;
	        in >> j;
	}
	
	in.close();


}
