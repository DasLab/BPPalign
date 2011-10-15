// DrawDoc.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "DrawDoc.h"
#include "../src/algorithm.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDrawDoc

IMPLEMENT_DYNCREATE(CDrawDoc, CCTDoc)

CDrawDoc::CDrawDoc() {

	

}

CDrawDoc::CDrawDoc(char *ctfilename, bool *Clockwise, CMainFrame *pMainFrame, CRNAstructureApp *App)
{
	CString name;

	clockwise = Clockwise;
	pmainframe = pMainFrame;

	name = ctfilename;
	openct(&ct,ctfilename);	
	SetTitle(name);

	//allocate space for the coordinates
	out = new coordinates(ct.numofbases);
	iscolorannotated = false;
	isSHAPEannotated = false;
	app = App;
	

}

void CDrawDoc::NewStructure(int number) {
	StructureNumber = number;
	int i;

	for (i=1;i<=ct.numofbases;i++) {
		if (ct.basepr[number][i]) {
			Calculate();
			nopair = false;
			return;
		}
	}
	nopair = true;
}

void CDrawDoc::Calculate() {
	//Calculate the coordinates for structure number StructureNumber

	//remove any pseudoknots, but keep the pairs in a list
	//findpseudo(&ct,StructureNumber,&npseudo,pseudo);

	PK = HasPseudo(&ct, StructureNumber); 	
	if (PK) {
		//use a circular drawing routine if there are pseudoknots
		placepk(&ct,out,height,width);

	}
	//use a standard diagram if there are no pseudoknots
	else place(StructureNumber,&ct,out,height,width);

	sortxy(out,!*clockwise,height,width);//counter specifies clockwise or counterclockwise


	xmax = out->x[1];
	ymax = out->y[1];

	for (int i=2;i<=(ct.numofbases);i++) {
   		if (out->x[i] > xmax) xmax = out->x[i];
		if (out->y[i] > ymax) ymax = out->y[i];
	}
	for (int i=10;i<=ct.numofbases;i+=10) {
		if (out->num[i/10][0] > xmax) {
			xmax = out->num[i/10][0];
		}
		if (out->num[i/10][1] > ymax) {
			ymax = out->num[i/10][1];
		}
	}

}

BOOL CDrawDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	return TRUE;
}



CDrawDoc::~CDrawDoc()
{

	delete out;
	if (iscolorannotated) delete[] color;

}


BEGIN_MESSAGE_MAP(CDrawDoc, CCTDoc)
	//{{AFX_MSG_MAP(CDrawDoc)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDrawDoc diagnostics

#ifdef _DEBUG
void CDrawDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CDrawDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CDrawDoc serialization

void CDrawDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}


void CDrawDoc::colorannotate() {
	//the structure will be color annotated

	color = new short [ct.numofbases+1];
	

}

void CDrawDoc::OutPutHelices(char *filename) {

	writehelixfile(filename,&ct,StructureNumber);

}
/////////////////////////////////////////////////////////////////////////////
// CDrawDoc commands
