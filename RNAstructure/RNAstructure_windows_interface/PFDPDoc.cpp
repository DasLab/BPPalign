#include "stdafx.h"
#include "pfdpdoc.h"
#include <math.h>

using namespace std;

IMPLEMENT_DYNCREATE(PFDPDoc, CDocument)

PFDPDoc::PFDPDoc():PlotDoc() {

}

PFDPDoc::PFDPDoc(CString Filename, CRNAstructureApp *app):PlotDoc(app)
{
	filename = Filename;
	PFPRECISION logint;
	int i,j;
	low = 1;
	high = -INFINITE_ENERGY;
	short vers;

	//allocate the ct file by reading the save file:
	ifstream sav(filename.GetBuffer(10),ios::binary);

	
	read(&sav,&(vers));//read the version of the save file
		//right now there is no infrastructure to indicate the wrong version is being read. 
		//This should be changed in the future...

	read(&sav,&(ct.numofbases));
	
	sav.close();
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
			array[i][j]=0;
		}
	}
	
	ct.allocate(ct.numofbases);

	w = new pfunctionclass(ct.numofbases);
	v = new pfunctionclass(ct.numofbases);
	wmb = new pfunctionclass(ct.numofbases);
	wmbl = new pfunctionclass(ct.numofbases);
	wcoax = new pfunctionclass(ct.numofbases);
	wl = new pfunctionclass(ct.numofbases);
	fce = new forceclass(ct.numofbases);

	w5 = new PFPRECISION [ct.numofbases+1];
	w3 = new PFPRECISION [ct.numofbases+2];

	lfce = new bool [2*ct.numofbases+1];
    mod = new bool [2*ct.numofbases+1];

	data = new pfdatatable();

	//load all the data from the pfsavefile:
	readpfsave(filename.GetBuffer(10), &ct, w5, w3,v, w, wmb,wl, wmbl, wcoax, fce,&scaling,mod,lfce,data);

	//fill array with the values for the plot:
	for (i=1;i<ct.numofbases;i++) {
		for (j=i+1;j<=ct.numofbases;j++) {

			if ((v->f(i,j)*v->f(j,i+ct.numofbases))>0) {
				logint = calculateprobability(i,j,v,w5,&ct,data,lfce,mod,scaling,fce);
				array[j][i] = -log10(logint);

				if (array[j][i]<low) low = array[j][i];
				if (array[j][i]>high) high = array[j][i];
			}
			else array[j][i] = 1e20;
		}
	}

	originallow = low;
	originalhigh = high;
   
	//select the ranges for each color in display
	colorranges();

	
	message = " -log10(BP Probability) ";
	SetTitle(Filename);

	outside = "infinity";

}

PFDPDoc::~PFDPDoc(void)
{
	int i;
	delete w;
	delete v;
	delete wmb;
	delete fce;
	delete[] w5;
	delete[] w3;
	for (i=1;i<=ct.numofbases;i++) {
		delete[] array[i];
	}
	delete[] array;
	delete[] lfce;
	delete[] mod;
	delete data;
	delete wmbl;
	delete wl;
	delete wcoax;


}
