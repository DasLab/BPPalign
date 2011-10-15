#include "stdafx.h"
#include "dpdoc.h"
#include "../src/algorithm.h"
#include "../src/defines.h"

using namespace std;

IMPLEMENT_DYNCREATE(DPDoc, CDocument)

DPDoc::DPDoc(void)
{
}

DPDoc::DPDoc(CString Filename, CRNAstructureApp *app, bool istextfile):PlotDoc(app)
{
	//Place all the energies in the array[j][i] array. (Inherited from PlotDoc.)
	//(Note that array is a double.)
	short vers;//short for saving the save file version	
	int i,j;
	PFPRECISION ene;

	low = INFINITE_ENERGY;
	high = 0;

	if (!istextfile) {
		register int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
		{0,1,0,1,0,0},{0,0,0,0,0,0}};
		
		arrayclass *w2,*wmb2;
		integersize *w5,*w3;
		int vmin;
		bool *lfce,*mod;
		datatable *data;

		data = new datatable;

		//peek at the length of the sequence and whether the folding is intermolecular to allocate arrays:
		ifstream sav(Filename.GetBuffer(10),ios::binary);
	
		read(&sav,&vers);
		if (vers!=safiversion) {
			//at a future time, add save file version checking here.
		}
		read(&sav,&(ct.numofbases));
		read(&sav,&(ct.intermolecular));
		sav.close();


		//allocate everything
		xstart = 1;
		xstop = ct.numofbases;
		ystart = 1;
		ystop = ct.numofbases;
		ct.allocate(ct.numofbases);
	
		arrayclass w(ct.numofbases);
		arrayclass v(ct.numofbases);
		arrayclass wmb(ct.numofbases);	
		forceclass fce(ct.numofbases);


		lfce = new bool [2*ct.numofbases+1];
		mod = new bool [2*ct.numofbases+1];
	
		w5 = new integersize [ct.numofbases+1];
		w3 = new integersize [ct.numofbases+2];

		if (ct.intermolecular) {
			w2 = new arrayclass(ct.numofbases);
			wmb2 = new arrayclass(ct.numofbases);
		
			for (i=0;i<3;i++) read(&sav,&(ct.inter[i]));
		
		}
		else {
			w2 = NULL;
			wmb2 = NULL;
		}

		array = new PFPRECISION *[ct.numofbases+1];
		for (i=1;i<=ct.numofbases;i++) {
			array[i] = new PFPRECISION [i+1];
		}
		for (i=1;i<=ct.numofbases;i++) {
			for (j=1;j<=i;j++) {
				array[i][j]=0;
			}
		}
	
		readsav(Filename.GetBuffer(10), &ct, w2, wmb2, w5, w3, lfce, mod, data,
			 &v, &w, &wmb, &fce, &vmin);

		for (i=1;i<=ct.numofbases;i++) {
			for (j=i+1;j<=ct.numofbases;j++) {

				array[j][i] = ((PFPRECISION) (v.f(i,j)+v.f(j,i+ct.numofbases)))/conversionfactor;
				if (array[j][i]<low) low = array[j][i];
				//if (array[j][i]>high) high = array[j][i];
			}
		}

		delete[] lfce;
		delete[] mod;
	



		delete[] w5;
		delete[] w3;


		if (ct.intermolecular) {
			delete w2;
			delete wmb2;
		}

		delete data;

	}
	else {//This is a plain text file
	
		ifstream *in;
		in = new ifstream();
		in->open(Filename.GetBuffer(0));
		
		*in >> ct.numofbases;

		//foo is a string used to temporarily store some info 
		char foo[maxfil];
		
		*in >>foo;//read the "i" header
		*in >>foo;//read the "j" header
		*in >>foo;//read the type of data

		//set message to the type of data
		message = foo;
		
		ct.allocate(ct.numofbases);
		
		for (i=1;i<=ct.numofbases;i++) ct.nucs[i]='X';
		
		
		//allocate everything
		xstart = 1;
		xstop = ct.numofbases;
		ystart = 1;
		ystop = ct.numofbases;

		array = new PFPRECISION *[ct.numofbases+1];
		for (i=1;i<=ct.numofbases;i++) {
			array[i] = new PFPRECISION [i+1];
		}
		for (i=1;i<=ct.numofbases;i++) {
			for (j=1;j<=i;j++) {
				array[i][j]=INFINITE_ENERGY;
			}
		}

		

		//read the file once to determine maximum i and j;
		*in >> i;
		*in >> j;
		*in >> ene;
		
		
		while(!in->eof()) {
			if (i<j) array[j][i]=ene;
			else array[i][j]=ene;
			if (array[j][i]<low) low = array[j][i];
			if (array[j][i]>high) high=array[j][i];
			*in >> i;
			*in >> j;
			*in >> ene;

			
			

		}

		in->close();
		delete in;

		
	}


	originallow = low;
	originalhigh = high;
   
	//select the ranges for each color in display
	colorranges();

	
	SetTitle(Filename);
	outside = "infinity";
}

DPDoc::~DPDoc(void)
{
	int i;
	for (i=1;i<=ct.numofbases;i++) {
		delete[] array[i];
	}
	delete[] array;
}
