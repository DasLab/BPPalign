// BimolPFView.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "BimolPFView.h"
#include "RNAstructure.h"
#include "temp_Dialog.h"

#include "../src/algorithm.h"

#include <direct.h>
#include "../src/pfunction.h"


// CBimolPFView

IMPLEMENT_DYNCREATE(CBimolPFView, CFormView)

CBimolPFView::CBimolPFView()
	: CFormView(CBimolPFView::IDD)
	, m_seq1(_T(""))
	, m_seq2(_T(""))
	, m_save(_T(""))
{
}

CBimolPFView::~CBimolPFView()
{
}

void CBimolPFView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_SEQUENCENAME, m_seq1);
	DDX_Text(pDX, IDC_SEQUENCE2NAME, m_seq2);
	DDX_Text(pDX, IDC_SAVENAME, m_save);
}

BEGIN_MESSAGE_MAP(CBimolPFView, CFormView)
	ON_BN_CLICKED(IDC_START, OnBnClickedStart)
	ON_BN_CLICKED(IDC_SEQUENCEBUTTON, OnBnClickedSequencebutton)
	ON_BN_CLICKED(IDC_SEQUENCE2BUTTON, OnBnClickedSequence2button)
	ON_BN_CLICKED(IDC_SAVEBUTTON, OnBnClickedSavebutton)
	ON_MESSAGE(ID_FOLDDONE, DoneFolding)
	ON_COMMAND(ID_TEMPERATURE, OnTemperature)
END_MESSAGE_MAP()


// CBimolPFView diagnostics

#ifdef _DEBUG
void CBimolPFView::AssertValid() const
{
	CFormView::AssertValid();
}

void CBimolPFView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif //_DEBUG

void CBimolPFView::OnUpdate(CView*, LPARAM, CObject*)
{

	ResizeParentToFit(FALSE);
	ResizeParentToFit();
	

}


// CBimolPFView message handlers

void CBimolPFView::OnBnClickedStart()
{
	int i;
	structure ct1, ct2;

	//Now do the work of the partition function calculation
	UpdateData(TRUE);
	
	if(m_seq1=="") {
		MessageBox("Please specify a sequence name.");
		return;
	}
	else if(m_seq2=="") {
		MessageBox("Please specify a sequence name.");
		return;
	}

	//check to see if the temperature has been changed.
	if (GetFoldDocument()->T<310||GetFoldDocument()->T>311) {

		//change the temperature from 310.15
		if (GetFoldDocument()->newtemp()==0) {
			//if newtemp returned zero, pass a warning to the user
			AfxMessageBox( "An enthalpy data file could not be found!\nTemperature of prediction will revert back to 37 degrees C.\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
				MB_OK|MB_ICONHAND);

		}

	}

	progress = new TProgressDialog();//parent,"Folding the RNA..."
	CRect *rect;
	rect = new CRect(10,40,210,120);
	
	if (GetFoldDocument()->ISRNA)
	progress->Create(NULL,"Folding the RNA...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);//(IDD_PROGRESSDIALOG);//
	else progress->Create(NULL,"Folding the DNA...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);
	delete rect;

	if((openseq (&ct1,m_seq1.GetBuffer(10)))==0||
		(openseq (&ct2,m_seq2.GetBuffer(10)))==0) {

		AfxMessageBox( "One of the sequences is unrecognizable.", 
			MB_OK|MB_ICONSTOP);

		return;

	}

	strcpy(GetFoldDocument()->ct.ctlabel[1],ct1.ctlabel[1]);
	//remove the new line at the end of ct.ctlabel[1]
	i = strlen(GetFoldDocument()->ct.ctlabel[1]);
	GetFoldDocument()->ct.ctlabel[1][i-1]='\0';

	strcat(GetFoldDocument()->ct.ctlabel[1],"_");
	strcat(GetFoldDocument()->ct.ctlabel[1],ct2.ctlabel[1]);

	//prepare ct3 with both sequences:
	GetFoldDocument()->ct.numofbases=ct1.numofbases+ct2.numofbases+3;
	GetFoldDocument()->ct.allocate(ct1.numofbases+ct2.numofbases+3);
	for (i=1;i<=ct1.numofbases;i++) {
		GetFoldDocument()->ct.numseq[i] = ct1.numseq[i];
		GetFoldDocument()->ct.nucs[i] = ct1.nucs[i];
		GetFoldDocument()->ct.hnumber[i] = ct1.hnumber[i];

	}
	
	for (i=1;i<=ct2.numofbases;i++) {
		GetFoldDocument()->ct.numseq[i+ct1.numofbases+3] = ct2.numseq[i];
		GetFoldDocument()->ct.nucs[i+ct1.numofbases+3] = ct2.nucs[i];
		GetFoldDocument()->ct.hnumber[i+ct1.numofbases+3] = ct2.hnumber[i];

	} 	
      
   
   GetFoldDocument()->ct.numseq[ct1.numofbases+1] = 5;
   GetFoldDocument()->ct.numseq[ct1.numofbases+2] = 5;
   GetFoldDocument()->ct.numseq[ct1.numofbases+3] = 5;

   GetFoldDocument()->ct.nucs[ct1.numofbases+1] = 'I';
   GetFoldDocument()->ct.nucs[ct1.numofbases+2] = 'I';
   GetFoldDocument()->ct.nucs[ct1.numofbases+3] = 'I';

   GetFoldDocument()->ct.hnumber[ct1.numofbases+1] = 0;
   GetFoldDocument()->ct.hnumber[ct1.numofbases+2] = 0;
   GetFoldDocument()->ct.hnumber[ct1.numofbases+3] = 0;


   GetFoldDocument()->ct.inter[0] = ct1.numofbases+1;
   GetFoldDocument()->ct.inter[1] = ct1.numofbases+2;
   GetFoldDocument()->ct.inter[2] = ct1.numofbases+3;

   GetFoldDocument()->ct.intermolecular = true;

   //Also copy information about nucleotides that must be single stranded
		//(These were entered as lowercase by the user.)

   for (i=1;i<=ct1.nnopair;i++) {
		GetFoldDocument()->ct.nnopair++;
		GetFoldDocument()->ct.nopair[GetFoldDocument()->ct.nnopair] = ct1.nopair[i];

   }
   for (i=1;i<=ct2.nnopair;i++) {
		GetFoldDocument()->ct.nnopair++;
		GetFoldDocument()->ct.nopair[GetFoldDocument()->ct.nnopair] = ct2.nopair[i]+ct1.numofbases+3;

   }


	pfobject.parent=this;
	
	pfobject.data=&GetFoldDocument()->data;
	pfobject.savefile = new char[m_save.GetLength()+1];
	strcpy(pfobject.savefile,m_save.GetBuffer(10));
	pfobject.progress = progress;
	pfobject.ct=&(GetFoldDocument()->ct);
	pfobject.Temperature=GetFoldDocument()->T;

	
	AfxBeginThread(PFProc,&pfobject);

	
	
}

void CBimolPFView::OnBnClickedSequencebutton()
{
	CFileDialog *filedialog;
		
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Sequence Files (*.seq)|*.seq||");

	
	filedialog->m_ofn.lpstrInitialDir=GetFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		
		m_seq1=(filedialog->GetPathName()).GetBuffer(30);
		
		
		
		UpdateData(FALSE);
		
		//now store the path in Startpath so that the program can start here next time:
		_getcwd(GetFoldDocument()->startpath,_MAX_PATH);

		
		

	}
	delete filedialog;
}

void CBimolPFView::OnBnClickedSequence2button()
{
	CFileDialog *filedialog;
	int i;
	char *seqname1,*seqname2;
		
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Sequence Files (*.seq)|*.seq||");

	
	filedialog->m_ofn.lpstrInitialDir=GetFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		
		m_seq2=(filedialog->GetPathName()).GetBuffer(30);
		
		if (m_seq1!="") {
			i = m_seq1.GetLength();
			seqname1= new char [(i)+1];
			i = m_seq2.GetLength();
			seqname2= new char [(i)+1];

			strcpy(seqname1,m_seq1.GetBuffer(1));
			strcpy(seqname2,m_seq2.GetBuffer(2));

			i = strlen(seqname1);
			while(i>=0){
			
				if (seqname1[i]=='.') break;
				i--;
			}
			if (i!=0) seqname1[i] = '\0';
			
			
			
			
			
			i = strlen(seqname2);
			while(i>=0){
			
				if (seqname2[i]=='.') break;
				i--;
			}
			if (i!=0) seqname2[i] = '\0';

			//get rid of the path from seqname2
			i = strlen(seqname2);
			while (i>=0) {
				if (seqname2[i]=='\\') break;
				i--;

			}
			strcpy(seqname2,seqname2+i+1);
			
			m_save+=seqname1;
			m_save+="_";
			m_save+=seqname2;
			m_save+=".pfs";

			delete[] seqname1;
			delete[] seqname2;
		}
		
		UpdateData(FALSE);
		
		//now store the path in Startpath so that the program can start here next time:
		_getcwd(GetFoldDocument()->startpath,_MAX_PATH);

		
		

	}
	delete filedialog;
}

void CBimolPFView::OnBnClickedSavebutton()
{
	//The user is specifying the save file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".pfs",m_save,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Partition Function Save Files (*.pfs)|*.pfs|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		m_save=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
}

CFoldDoc *CBimolPFView::GetFoldDocument() {
	
	return ((CFoldDoc*) GetDocument());	

}

LRESULT CBimolPFView::DoneFolding(WPARAM wParam, LPARAM lParam) {
	delete[] pfobject.savefile;

	progress->SendMessage (WM_CLOSE);

	((CRNAstructureApp*)GetFoldDocument()->pMainFrame)->BoxPlot(m_save);
	
	GetFoldDocument()->Frame->SendMessage(WM_CLOSE);

	
	return 0;
}

void CBimolPFView::OnTemperature()
{
	//Allow the user to specify a new temperature.
	CTemp_Dialog *temp;
	temp=new CTemp_Dialog(&(GetFoldDocument()->T));

	temp->DoModal();
	delete temp;
	


}
