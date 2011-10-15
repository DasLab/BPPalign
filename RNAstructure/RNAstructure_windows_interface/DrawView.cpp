// DrawView.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "DrawView.h"
#include "zoom.h"
#include "../src/algorithm.h"

#include <afxwin.h>
#include <Windows.h>


#include <iostream>
using namespace std;

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#define SCROLLCHANGE 5

/////////////////////////////////////////////////////////////////////////////
// CDrawView

IMPLEMENT_DYNCREATE(CDrawView, CScrollView)

CDrawView::CDrawView()
{

	resize = true;
	
}

CDrawView::~CDrawView()
{

	if (pDoc->iscolorannotated) {
		
		delete[] crColor;
		delete[] w5;
		delete v;
	}
	delete Font;

}

void CDrawView::NewSize(int Size) {
	
	pDoc = (CDrawDoc*) GetDocument();

	pDoc->FontSize = Size;
	NewFont();

}


void CDrawView::NewFont() {
	pDoc = (CDrawDoc*) GetDocument();
	LOGFONT pLogFont;
	
	pLogFont.lfHeight = pDoc->FontSize;
	pLogFont.lfWidth = 0;
	pLogFont.lfEscapement = 0;
	pLogFont.lfOrientation = 0;
	pLogFont.lfWeight = FW_NORMAL;
	pLogFont.lfItalic = 0;
	pLogFont.lfUnderline = 0;
	pLogFont.lfStrikeOut = 0;
	pLogFont.lfCharSet = ANSI_CHARSET;
	pLogFont.lfOutPrecision = OUT_DEFAULT_PRECIS;
	pLogFont.lfClipPrecision = CLIP_DEFAULT_PRECIS;
	pLogFont.lfQuality = PROOF_QUALITY;
	pLogFont.lfPitchAndFamily = VARIABLE_PITCH|FF_MODERN;
	strcpy(pLogFont.lfFaceName,"Courier New");

	Font= new CFont;
	Font->CreateFontIndirect(&pLogFont);

	pDoc->height = pDoc->FontSize;
	pDoc->width = pDoc->FontSize/2;

	

	//Font.CreatePointFont( pDoc->FontSize, "Courier New",  NULL );

	



}

void CDrawView::NewStructure(int number) {
	pDoc = (CDrawDoc*) GetDocument();
	
	pDoc->NewStructure(number);
	
	


}


BEGIN_MESSAGE_MAP(CDrawView, CScrollView)
	//{{AFX_MSG_MAP(CDrawView)
	ON_COMMAND(ID_DRAW_ZOOM, OnDrawZoom)
	ON_COMMAND(ID_DRAW_STRUCTURENUMBER, OnDrawStructurenumber)
	ON_WM_KEYDOWN()
	ON_COMMAND(ID_DRAW_CLOCKWISE, OnDrawClockwise)
	ON_COMMAND(ID_FILE_PRINT, OnFilePrint)
	ON_COMMAND(ID_FILE_PRINTPREVIEW, OnFilePrintPreview)
	ON_COMMAND(ID_CPY, OnCpy)
	//}}AFX_MSG_MAP
	ON_COMMAND(ID_DRAW_ADDCOLORANNOTATION, OnDrawAddcolorannotation)
	ON_COMMAND(ID_DRAW_SHOWCOLORANNOTATIONKEY, OnDrawShowcolorannotationkey)
	ON_COMMAND(ID_DRAW_EXPORTSTRUCTURETOTEXTFILE, OnDrawExportstructuretotextfile)
	ON_COMMAND(ID_DRAW_ADDSHAPEANNOTATION, OnDrawAddshapeannotation)
	ON_COMMAND(ID_DRAW_SHOWSHAPEANNOTATIONKEY, OnDrawShowshapeannotationkey)
	ON_COMMAND(ID_DRAW_EXPORTTODOT, &CDrawView::OnDrawExporttodot)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDrawView drawing

void CDrawView::OnInitialUpdate()
{
	CScrollView::OnInitialUpdate();
	CSize sizeTotal;


	pDoc = (CDrawDoc*) GetDocument();
	NewSize(16);
	

	NewStructure(1);
	

	
	

	zoom = 100;
	//set dummy scroll size for now:
	sizeTotal.cx = 10;
	sizeTotal.cy = 10;
	SetScrollSizes(MM_TEXT, sizeTotal);




}


	

	





void CDrawView::OnDraw(CDC* pDC)
{
	
	CString title;
	CSize  oldVExt, oldWExt;
	int prevMode,ht;
	CSize sizeTotal;
	char number[8];
	
	CRect rect;
	CPoint topleft,bottomright;
	int zoom2;

	if (pDoc->iscolorannotated) DetermineColor();
	if (pDoc->isSHAPEannotated)	DetermineSHAPEColor();

	if (pDoc->nopair) {
		itoa(pDoc->StructureNumber,number,10);
		title = "Structure #";
		title += number;
		title += "  ";
		title += pDoc->ct.ctlabel[pDoc->StructureNumber];
      	title += "  Structure contains no pair.";
      	pDC->TextOut(0, 0, title);
        return;
    }
	prevMode = pDC->SetMapMode(MM_ANISOTROPIC);
	
	CFont *oldFont = pDC->SelectObject(Font);

	itoa(pDoc->StructureNumber,number,10);
	title = "Structure #";
	title += number;
	title += " of ";
	itoa(pDoc->ct.numofstructures,number,10);
	title += number;
	title += "   ";
	title += pDoc->ct.ctlabel[pDoc->StructureNumber];

	//Add the title of the sequence and the current structure number	
	pDC->TextOut(0, 0, title);

	if (resize) {
		//zoom out so that the whole structure fits in the client area
		GetClientRect( rect );

		topleft = rect.TopLeft();
		bottomright = rect.BottomRight();
		ht = bottomright.y-topleft.y;
		if ((ht)<pDoc->ymax) {

			zoom = int (100.0*((float (ht)) / (float (pDoc->ymax))) );
			resize=false;
		}

		
		ht = bottomright.x-topleft.x;
		if ((ht)<pDoc->xmax) {

			zoom2 = int (100.0*((float (ht)) / (float (pDoc->xmax))) );
			if (zoom2<zoom) zoom = zoom2;
			resize=false;

		}

	}

	

	CSize PageSize(int (float(pDoc->xmax)*((float (zoom)/float (100)))),
   	int (float(pDoc->ymax)*((float (zoom)/float (100)))));
	CSize WindowSize(pDoc->xmax,pDoc->ymax);
	
	//set the scroll bar sizes:
	//sizeTotal.cx = pDoc->xmax*((long (zoom)/long (100)))+15;
	//sizeTotal.cy = pDoc->ymax*((long (zoom)/long (100)))+15;
	//SetScrollSizes(MM_TEXT, sizeTotal);
	sizeTotal.cx = PageSize.cx+15;
	sizeTotal.cy = PageSize.cy+15;
	SetScrollSizes(MM_TEXT, sizeTotal);



	oldVExt = pDC->SetViewportExt(PageSize.cx,PageSize.cy);

	oldWExt = pDC->SetWindowExt(WindowSize.cx,WindowSize.cy); 
	
	
	if (resize) {
		resize = false;

		//shrink the view if the structure is smaller 
		ResizeParentToFit();


		

	   	
	}

	
	
	draw(pDC);

	
	pDC->SelectObject(oldFont);

	pDC->SetMapMode(prevMode);

	pDC->SetViewportExt(oldVExt);


	pDC->SetWindowExt(oldWExt);


}


void CDrawView::draw(CDC* pDC) {

	CPen *Pen;
	short i;
	COLORREF startcolor;

	CString title;
	
	char number[8];

	//draw the backbone connection
	Pen = new CPen(PS_SOLID,1,RGB(0,0,0));//pick a skinny, black pen
   	CPen *Oldpen = pDC->SelectObject(Pen);
   	pDC->MoveTo(pDoc->out->x[1],pDoc->out->y[1]);
   	for (i=2;i<=pDoc->ct.numofbases;i++) {//draw in the backbone connection
      	if (pDoc->ct.numseq[i]!=5) {
         	if (pDoc->ct.numseq[i-1]==5) {
            	pDC->MoveTo(pDoc->out->x[i]+pDoc->width/2,pDoc->out->y[i]+pDoc->height/2);
            }
            else {
   				pDC->LineTo(pDoc->out->x[i]+pDoc->width/2,pDoc->out->y[i]+pDoc->height/2);
            }
         }
   	}
	
	//place the nucleotide labels only if this is not an intermolecular fold
    if (!pDoc->ct.intermolecular) {
   		for (i=10;i<=pDoc->ct.numofbases;i=i+10) {//draw number to base conection
			if (!(pDoc->out->num[i/10][0]==0&&pDoc->out->num[i/10][1]==0)) {
				//Filter out cases where the coordinates are 0,0.
				//This is a flag to indicate that the label conflicts with a nuc position, so it 
					//should not be placed.
   				pDC->MoveTo(pDoc->out->x[i]+pDoc->width/2,pDoc->out->y[i]+pDoc->height/2);
      			pDC->LineTo(pDoc->out->num[i/10][0]+pDoc->width/2,pDoc->out->num[i/10][1]+pDoc->height/2);
   
   				itoa(i,number,10);
    			title = number;
      			pDC->TextOut(pDoc->out->num[i/10][0],pDoc->out->num[i/10][1], title);
			}
      	}
      }


	delete Pen;
	Pen = new CPen(PS_SOLID,3,RGB(0,0,0));//now use a fat black pen
	pDC->SelectObject(Pen);

	//draw in the basepairs
	if (pDoc->PK) {
		//If there was a pseudoknot in the structure, this is a circle and it requires bezier lines

		for (i=1;i<=pDoc->ct.numofbases;i++) {
   			if ((pDoc->ct.basepr[pDoc->StructureNumber][i])>i) {//draw in pairs
      			
				int numBasesBetween = pDoc->ct.basepr[pDoc->StructureNumber][i]-i;
				if( numBasesBetween > pDoc->ct.numofbases / 2 ) numBasesBetween = pDoc->ct.numofbases - numBasesBetween;


				int x1 = pDoc->out->x[i]+pDoc->width/2;
				int x2 = pDoc->out->x[pDoc->ct.basepr[pDoc->StructureNumber][i]]+pDoc->width/2;
				int y1 = pDoc->out->y[i]+pDoc->height/2;
				int y2 = pDoc->out->y[pDoc->ct.basepr[pDoc->StructureNumber][i]]+pDoc->height/2;

				int centerx = (pDoc->out->x[1] + pDoc->out->x[pDoc->ct.numofbases/2])/2;
				int centery = (pDoc->out->y[1] + pDoc->out->y[pDoc->ct.numofbases/2])/2;

				double diameter = sqrt(pow((double)(pDoc->out->x[1] - pDoc->out->x[pDoc->ct.numofbases/2]),2)+
					pow((double)(pDoc->out->y[1] - pDoc->out->y[pDoc->ct.numofbases/2]),2));

				int midpointx = ( x1 + x2 ) / 2;
				int midpointy = ( y1 + y2 ) / 2;
				double gamma = 0.9;
				double centerThreshold = diameter / 8.0 * 5.0;


				double distance = sqrt(pow((double)(x1 - x2),2)+
					pow((double)(y1 - y2),2));
				double maxBendDistance = (double)(numBasesBetween) * 2.0 / 
					((double)pDoc->ct.numofbases) * (diameter/2) * gamma;
				double lineAngle = atan2( ((double)(centery - midpointy)) , ((double)(centerx - midpointx)) );
				double distX = maxBendDistance * cos( lineAngle ) * 2;
				double distY = maxBendDistance * sin( lineAngle ) * 2;


				int controlPointX = midpointx + (int) distX;
				int controlPointY = midpointy + (int) distY;


				if( distance >= centerThreshold ) {
				   controlPointX = centerx;
				   controlPointY = centery;
				}

				POINT points[4];
				points[0].x=(long) x1;
				points[0].y=(long) y1;
				points[1].x=(long) controlPointX;
				points[1].y=(long) controlPointY;
				points[2].x=(long) controlPointX;
				points[2].y=(long) controlPointY;
				points[3].x=x2;
				points[3].y=y2;

				




				pDC->PolyBezier(points,4);
				
				
				
      		}
   		}


	}
	else {
		//There were no pseudoknots, draw lines
		for (i=1;i<=pDoc->ct.numofbases;i++) {
   			if ((pDoc->ct.basepr[pDoc->StructureNumber][i])>i) {//draw in pairs
      			pDC->MoveTo(pDoc->out->x[i]+pDoc->width/2,pDoc->out->y[i]+pDoc->height/2);
         		pDC->LineTo(pDoc->out->x[pDoc->ct.basepr[pDoc->StructureNumber][i]]+pDoc->width/2,
         			pDoc->out->y[pDoc->ct.basepr[pDoc->StructureNumber][i]]+pDoc->height/2);
      		}
   		}
	}
	
	delete Pen;

	if (pDoc->iscolorannotated||pDoc->isSHAPEannotated) {
				 startcolor = SetTextColor(
					pDC->m_hDC ,           
					crColor[pDoc->color[1]]);

	}
	//Place each nucleotide:
	for (i=1;i<=pDoc->ct.numofbases;i++) {
		if (pDoc->ct.numseq[i]!=5) {
			if (pDoc->iscolorannotated||pDoc->isSHAPEannotated) {
				 SetTextColor(
					pDC->m_hDC ,           
					crColor[pDoc->color[i]]);

			}
			title = pDoc->ct.nucs[i];
			pDC->TextOut(pDoc->out->x[i],pDoc->out->y[i],title);
		}

	}  	



	//Restore the DC to former settings:
	pDC->SelectObject(Oldpen);
	
	if (pDoc->iscolorannotated) SetTextColor(pDC->m_hDC ,startcolor);
	
	//pDC->SetMapMode(oldMapMode);

}


/////////////////////////////////////////////////////////////////////////////
// CDrawView diagnostics

#ifdef _DEBUG
void CDrawView::AssertValid() const
{
	CScrollView::AssertValid();
}

void CDrawView::Dump(CDumpContext& dc) const
{
	CScrollView::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CDrawView message handlers

void CDrawView::OnDrawZoom() 
{
	//Open the Zoom Dialog
	CZoom *zoomdialog;

	zoomdialog = new CZoom(&zoom,this);

	zoomdialog->Create(IDD_ZOOM,this);
	
}

void CDrawView::OnDrawStructurenumber() 
{
	CNumber *numberdialog;

	numberdialog = new CNumber((CDrawDoc*) GetDocument(),this);

	numberdialog->Create(IDD_STRUCTURENUMBER,this);
	
}





void CDrawView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags) 
{
	short control;
	CPoint point = GetScrollPosition( );
	POINT pt;

	control = (GetKeyState(VK_CONTROL)&(-128));

      switch (nChar)
         {
         case VK_UP:
            if (control) Cm_UP();
            else {
				if (point.y!=0) {
					pt.x = point.x;
					pt.y = point.y - SCROLLCHANGE;
					ScrollToPosition(pt);
				}
            }
            break;
         case VK_DOWN:
            if (control) Cm_DOWN();
         	else {
            	//if (point.y!=0) {
					pt.x = point.x;
					pt.y = point.y + SCROLLCHANGE;
					ScrollToPosition(pt);
				//}
            }
            break;
         case VK_LEFT:
         	if (control) Cm_SMALLER();
            else {
            	if (point.x!=0) {
					pt.x = point.x - SCROLLCHANGE;
					pt.y = point.y ;
					ScrollToPosition(pt);
				}
            }
            break;
         case VK_RIGHT:
         	if (control) Cm_LARGER();
            else {
            	//if (point.y!=0) {
					pt.x = point.x + SCROLLCHANGE;
					pt.y = point.y ;
					ScrollToPosition(pt);
				//}
            }
            break;
         case VK_HOME:
            //if (point.y!=0) {
					pt.x = 0;
					pt.y = 0;
					ScrollToPosition(pt);
				//}
            break;
         case VK_END:
            if (point.y!=0) {
					pt.x = point.x;
					pt.y = point.y - SCROLLCHANGE;
					ScrollToPosition(pt);
				}
            break;
         case VK_PRIOR:
            if (point.y!=0) {
					pt.x = point.x;
					pt.y = point.y + 10*SCROLLCHANGE;
					ScrollToPosition(pt);
				}
            break;
         case VK_NEXT:
            if (point.y!=0) {
					pt.x = point.x;
					pt.y = point.y - 10*SCROLLCHANGE;
					ScrollToPosition(pt);
				}
            break;

         //case '2':
         	//if (Scroller) Scroller->VScroll(SB_LINEDOWN,0);

         }
	
	CScrollView::OnKeyDown(nChar, nRepCnt, nFlags);
}

void CDrawView::Cm_SMALLER() {

	if ((zoom - 10) >0) {
		zoom = zoom - 10;
		Invalidate();
	}

}


void CDrawView::Cm_LARGER() {
	
		zoom = zoom + 10;
		Invalidate();
	


}

void CDrawView::Cm_UP() {

	if (pDoc->StructureNumber<pDoc->ct.numofstructures) {
		pDoc->NewStructure(pDoc->StructureNumber+1);
		Invalidate();
	}

}

void CDrawView::Cm_DOWN() {

	if (pDoc->StructureNumber>1) {
		pDoc->NewStructure(pDoc->StructureNumber-1);
		Invalidate();
	}


}

void CDrawView::OnDrawClockwise() 
{

	CMenu* menu = pDoc->pmainframe->GetMenu( );
	
	*pDoc->clockwise = !*pDoc->clockwise;
	if (*pDoc->clockwise) menu->CheckMenuItem(ID_DRAW_CLOCKWISE,MF_CHECKED);
	else menu->CheckMenuItem(ID_DRAW_CLOCKWISE,MF_UNCHECKED); 


	pDoc->NewStructure(pDoc->StructureNumber);

	Invalidate();
	
}

void CDrawView::OnFilePrint() 
{
	CView::OnFilePrint();
	
}

void CDrawView::OnFilePrintPreview() 
{
	CView::OnFilePrintPreview();
	
}

BOOL CDrawView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CDrawView::OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo)
{
	//Scale to one page:
	pInfo->SetMaxPage(1);
}

void CDrawView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}

void CDrawView::OnPrint( CDC *pDC, CPrintInfo *pInfo ) {
	int pageheight, pagewidth;
	float scalev,scaleh,scale;


	if (pDoc->nopair) {
      	AfxMessageBox( "Structure Contains No Pairs.",MB_OK|MB_ICONSTOP);
        return;
    }

	//Scale the image to fit one page:
	pageheight = pDC->GetDeviceCaps(VERTRES);
	pagewidth = pDC->GetDeviceCaps(HORZRES);
	scalev = float(pageheight)/float(pDoc->ymax);
	scaleh = float(pagewidth)/float(pDoc->ymax);

	if(scalev<scaleh) scale = scalev;
	else scale = scaleh;

	
	
	
	//OnDraw( pDC );//true indicates printing
	CString title;
	CSize  oldVExt, oldWExt;
	int prevMode;
	CSize sizeTotal;
	char number[8];
	
	CRect rect;
	CPoint topleft,bottomright;
	



	
	prevMode = pDC->SetMapMode(MM_ISOTROPIC);
	
	CFont *oldFont = pDC->SelectObject(Font);

	




	oldWExt = pDC->SetWindowExt(pDoc->xmax+500,pDoc->ymax);
	oldVExt = pDC->SetViewportExt(pageheight,pagewidth);

	
	itoa(pDoc->StructureNumber,number,10);
	title = "Structure #";
	title += number;
	title += "  ";
	title += pDoc->ct.ctlabel[pDoc->StructureNumber];

	//Add the title of the sequence and the current structure number	
	pDC->TextOut(0, 0, title); 
	

	draw(pDC);

	
	pDC->SelectObject(oldFont);

	
	pDC->SetWindowExt(oldWExt);
	pDC->SetViewportExt(oldVExt);
	pDC->SetMapMode(prevMode);


	


}

void CDrawView::OnCpy() 
{
	
    CBitmap     cBmp;
    CClientDC   cWndDC(this);   // View is an hWnd, so we can use "this"
    CDC         cMemDC;         // Handle to a memory DC
    CRect     rect;             // For storing the size of the window
	CBrush blankbrush;
	CRgn region;

    cMemDC.CreateCompatibleDC(&cWndDC); // Create the memory DC.

    //GetWindowRect(rect);         // Get the size of the window
    // Here we are going to do a little drawing so we can see a line.
    //cWndDC.MoveTo(0,0);
    //cWndDC.LineTo(rect.Width(),rect.Height());
	
	
	//OnDraw(&cMemDC);
	//int foo = cMemDC.GetBkColor( ); 
	//cMemDC.SetBkColor( RGB(255,255,255));
	//blankbrush.CreateSolidBrush(RGB(255,255,255));
	//region.CreateRectRgn(0,0,pDoc->xmax+2*pDoc->width,pDoc->ymax+2*pDoc->height );
	//cMemDC.FillRgn(&region,&blankbrush);
    //cBmp.CreateCompatibleBitmap(&cWndDC, pDoc->xmax+2*pDoc->width,pDoc->ymax+2*pDoc->height );


	if (pDoc->iscolorannotated) cBmp.CreateCompatibleBitmap(&cWndDC,pDoc->xmax+2*pDoc->width, pDoc->ymax+2*pDoc->height);
	else cBmp.CreateBitmap(pDoc->xmax+2*pDoc->width, pDoc->ymax+2*pDoc->height, 1, 1, NULL );

    // Keep the old bitmap
	CBitmap* pOldBitmap = cMemDC.SelectObject(&cBmp); 
	
	cMemDC.SetBkColor( RGB(255,255,255));
	blankbrush.CreateSolidBrush(RGB(255,255,255));
	region.CreateRectRgn(0,0,pDoc->xmax+2*pDoc->width,pDoc->ymax+2*pDoc->height );
	cMemDC.FillRgn(&region,&blankbrush);
    //cMemDC.BitBlt(0, 0, pDoc->xmax+15,pDoc->ymax+15, &cWndDC, 0, 0,
	//	SRCCOPY);
	draw(&cMemDC);
//pDoc->xmax,pDoc->ymax

    // here are the actual clipboard functions.
    AfxGetApp()->m_pMainWnd->OpenClipboard() ;
    EmptyClipboard() ;
    HANDLE foo = SetClipboardData (CF_BITMAP, cBmp.GetSafeHandle() ) ;
/*	LPVOID lpMsgBuf;
FormatMessage( 
    FORMAT_MESSAGE_ALLOCATE_BUFFER | 
    FORMAT_MESSAGE_FROM_SYSTEM | 
    FORMAT_MESSAGE_IGNORE_INSERTS,
    NULL,
    GetLastError(),
    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
    (LPTSTR) &lpMsgBuf,
    0,
    NULL 
);
MessageBox( NULL, (LPCTSTR)lpMsgBuf, MB_OK | MB_ICONINFORMATION );
// Free the buffer.
LocalFree( lpMsgBuf );*/
 



    CloseClipboard () ;
        // next we select the old bitmap back into the memory DC
        // so that our bitmap is not deleted when cMemDC is destroyed.
        // Then we detach the bitmap handle from the cBmp object so that
        // the bitmap is not deleted when cBmp is destroyed.
    cMemDC.SelectObject(pOldBitmap);
    cBmp.Detach();

    return;


	
}

void CDrawView::SetColors() {
		crColor = new COLORREF[10];
		crColor[0]=RGB(0,0,0);
		crColor[1]=RGB(colorintensity,0,0);
		crColor[2]=RGB(colorintensity,colorintensity/2,0);
		crColor[3]=RGB(4*colorintensity/5,3*colorintensity/4,0);
		crColor[4]=RGB(colorintensity/2,colorintensity/2,0);
		crColor[5]=RGB(0,colorintensity,0);
		crColor[6]=RGB(0,colorintensity,colorintensity);
		crColor[7]=RGB(0,0,colorintensity);
		crColor[8]=RGB(colorintensity,0,colorintensity);
		crColor[9]=RGB(colorintensity*2/3,colorintensity*2/3,colorintensity*2/3);
}

void CDrawView::OnDrawAddcolorannotation()
{

	//Add color annotation to the structure based on a partition function calculation.

	short vers;

	//Open the file dialog:
	CFileDialog *filedialog;
	filedialog = new CFileDialog(TRUE,".ct",NULL,OFN_FILEMUSTEXIST,
		"Partition Function Save Files (*.pfs)|*.pfs|All Files|*.*||");
	//filedialog->m_ofn.lpstrInitialDir=startpath;
	if (filedialog->DoModal()==IDOK) {
		
	
		
		
		structure ct;
		pfunctionclass *w,*wmb,*wmbl,*wcoax,*wl;
		forceclass *fce;
		PFPRECISION *w3;
		bool *lfce,*mod;
		pfdatatable *data;
		

		//allocate the ct file by reading the save file:
		ifstream sav(filedialog->GetPathName().GetBuffer(10),ios::binary);

		read(&sav,&vers);
		if (vers!=pfsaveversion) {
			//the file version doesn't match the current version
			AfxMessageBox( "Error: This save file was created with a different version of RNAstructure.", 
			MB_OK|MB_ICONEXCLAMATION);
			
			delete filedialog;

			return;

		}
	
		read(&sav,&(ct.numofbases));

		if (ct.numofbases!=pDoc->ct.numofbases) {
			//the partition function save file has a different number of nucs than the currently drawn file
			AfxMessageBox( "Error: The partition function calculation is for a different length sequence.", 
			MB_OK|MB_ICONEXCLAMATION);
			
			delete filedialog;

			return;
			

		}

		


		sav.close();
		//allocate everything
		
		data = new pfdatatable;
	
		ct.allocate(ct.numofbases);

		w = new pfunctionclass(ct.numofbases);
		v = new pfunctionclass(ct.numofbases);
		wmb = new pfunctionclass(ct.numofbases);
		fce = new forceclass(ct.numofbases);
		wl = new pfunctionclass(ct.numofbases);
		wcoax = new pfunctionclass(ct.numofbases);
		wmbl = new pfunctionclass(ct.numofbases);

		w5 = new PFPRECISION [ct.numofbases+1];
		w3 = new PFPRECISION [ct.numofbases+2];

		lfce = new bool [2*ct.numofbases+1];
		mod = new bool [2*ct.numofbases+1];

		//load all the data from the pfsavefile:
		readpfsave(filedialog->GetPathName().GetBuffer(10), &ct, w5, w3,v, w, wmb,wl, wmbl, wcoax, fce,&scaling,mod,lfce,data);

		//fill array with the values for the plot:
		pDoc->colorannotate();
		pDoc->iscolorannotated=true;

		DetermineColor();

		SetColors();


		
		delete[] w3;
		delete w;
		
		delete wmb;
		delete fce;
		delete[] lfce;
		delete[] mod;
		delete data;
		delete wcoax;
		delete wmbl;
		delete wl;

	}
	delete filedialog;


	//recreate the font to bold so that color annotation is more visible
	LOGFONT pLogFont;
	
	pLogFont.lfHeight = pDoc->FontSize;
	pLogFont.lfWidth = 0;
	pLogFont.lfEscapement = 0;
	pLogFont.lfOrientation = 0;
	pLogFont.lfWeight = FW_BOLD;
	pLogFont.lfItalic = 0;
	pLogFont.lfUnderline = 0;
	pLogFont.lfStrikeOut = 0;
	pLogFont.lfCharSet = ANSI_CHARSET;
	pLogFont.lfOutPrecision = OUT_DEFAULT_PRECIS;
	pLogFont.lfClipPrecision = CLIP_DEFAULT_PRECIS;
	pLogFont.lfQuality = PROOF_QUALITY;
	pLogFont.lfPitchAndFamily = VARIABLE_PITCH|FF_MODERN;
	strcpy(pLogFont.lfFaceName,"Courier New");

	delete Font;
	Font = new CFont;
	Font->CreateFontIndirect(&pLogFont);

	
	Invalidate();
}

void CDrawView::OnDrawShowcolorannotationkey()
{
	pDoc->app->DisplayColorKeyWindow();
}

void CDrawView::OnDrawExportstructuretotextfile()
{
	//Output the current structure to a text file that can be read by XRNA
	//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".txt",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Helix Files (*.txt)|*.txt|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		
		pDoc->OutPutHelices((filedialog->GetPathName()).GetBuffer(0));		

	}
	delete filedialog;
}

void CDrawView::DetermineColor() {
int i;
		for (i=1;i<=pDoc->ct.numofbases;i++) {
			//stratify into 9 classes of color - hard-wired currently
			if (pDoc->ct.basepr[pDoc->StructureNumber][i]==0) pDoc->color[i]=0;
			else if (pDoc->ct.basepr[pDoc->StructureNumber][i]>i) {
				if ((v->f(i,pDoc->ct.basepr[pDoc->StructureNumber][i])*v->f(pDoc->ct.basepr[pDoc->StructureNumber][i],i+pDoc->ct.numofbases))/(w5[pDoc->ct.numofbases]*scaling*scaling)>=.99) {
					pDoc->color[i]=1; 
					pDoc->color[pDoc->ct.basepr[pDoc->StructureNumber][i]]=1;
				}
				else if ((v->f(i,pDoc->ct.basepr[pDoc->StructureNumber][i])*v->f(pDoc->ct.basepr[pDoc->StructureNumber][i],i+pDoc->ct.numofbases))/(w5[pDoc->ct.numofbases]*scaling*scaling)>=.95) {
					pDoc->color[i]=2;
					pDoc->color[pDoc->ct.basepr[pDoc->StructureNumber][i]]=2;
				}
				else if ((v->f(i,pDoc->ct.basepr[pDoc->StructureNumber][i])*v->f(pDoc->ct.basepr[pDoc->StructureNumber][i],i+pDoc->ct.numofbases))/(w5[pDoc->ct.numofbases]*scaling*scaling)>=.90) {
					pDoc->color[i]=3;
					pDoc->color[pDoc->ct.basepr[pDoc->StructureNumber][i]]=3;
				}
				else if ((v->f(i,pDoc->ct.basepr[pDoc->StructureNumber][i])*v->f(pDoc->ct.basepr[pDoc->StructureNumber][i],i+pDoc->ct.numofbases))/(w5[pDoc->ct.numofbases]*scaling*scaling)>=.80) {
					pDoc->color[i]=4;
					pDoc->color[pDoc->ct.basepr[pDoc->StructureNumber][i]]=4;
				}
				else if ((v->f(i,pDoc->ct.basepr[pDoc->StructureNumber][i])*v->f(pDoc->ct.basepr[pDoc->StructureNumber][i],i+pDoc->ct.numofbases))/(w5[pDoc->ct.numofbases]*scaling*scaling)>=.70) {
					pDoc->color[i]=5;
					pDoc->color[pDoc->ct.basepr[pDoc->StructureNumber][i]]=5;
				}
				else if ((v->f(i,pDoc->ct.basepr[pDoc->StructureNumber][i])*v->f(pDoc->ct.basepr[pDoc->StructureNumber][i],i+pDoc->ct.numofbases))/(w5[pDoc->ct.numofbases]*scaling*scaling)>=.60) {
					pDoc->color[i]=6;
					pDoc->color[pDoc->ct.basepr[pDoc->StructureNumber][i]]=6;
				}
				else if ((v->f(i,pDoc->ct.basepr[pDoc->StructureNumber][i])*v->f(pDoc->ct.basepr[pDoc->StructureNumber][i],i+pDoc->ct.numofbases))/(w5[pDoc->ct.numofbases]*scaling*scaling)>.50) {
					pDoc->color[i]=7;
					pDoc->color[pDoc->ct.basepr[pDoc->StructureNumber][i]]=7;
				}
				else {
					pDoc->color[i]=8;
					pDoc->color[pDoc->ct.basepr[pDoc->StructureNumber][i]]=8;
				}
			}
		}
}

//Function to determine drawing colors with SHAPE annotation
void CDrawView::DetermineSHAPEColor() {
	int i;
		for (i=1;i<=pDoc->ct.numofbases;i++) {
			//stratify into 7 classes of color - hard-wired currently
			if (pDoc->ct.SHAPE[i]<-900) pDoc->color[i]=9;
			else  {
				if (pDoc->ct.SHAPE[i]>=0.7) {
					pDoc->color[i]=1; 
					
				}
				else if (pDoc->ct.SHAPE[i]>=0.3) {
					pDoc->color[i]=2;
					
				}
				//else if (pDoc->ct.SHAPE[i]>=1.0) {
				//	pDoc->color[i]=3;
					
				//}
				
				//else if (pDoc->ct.SHAPE[i]>=0.5) {
				//	pDoc->color[i]=5;
					
				//}
				//else if (pDoc->ct.SHAPE[i]>=0.0) {
				//	pDoc->color[i]=6;
					
				//}
				else {//if (pDoc-ct->SHAPE[]>=2.0) {
					pDoc->color[i]=0;
					
				}
				
			}
		}
}

void CDrawView::OnDrawAddshapeannotation()
{
	//Add a color annotation based on SHAPE data
	CString InputFilename;


	//Get the input filename
	CFileDialog *filedialog;

	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"SHAPE Files (*.shape)|*.shape|All Files (*.*)|*.*||");

	
	if (filedialog->DoModal()==IDOK) {
		//Good result from file choice -- proceed.
		InputFilename=(filedialog->GetPathName()).GetBuffer(30);

		pDoc->ct.ReadSHAPE(InputFilename.GetBuffer(40),false);

		pDoc->colorannotate();
		pDoc->isSHAPEannotated=true;

		DetermineSHAPEColor();

		SetColors();

		/*crColor = new COLORREF[9];
		crColor[0]=RGB(0,0,0);
		crColor[1]=RGB(colorintensity,0,0);
		crColor[2]=RGB(colorintensity,colorintensity/2,0);
		crColor[3]=RGB(4*colorintensity/5,3*colorintensity/4,0);
		crColor[4]=RGB(colorintensity/2,colorintensity/2,0);
		crColor[5]=RGB(0,colorintensity,0);
		crColor[6]=RGB(0,colorintensity,colorintensity);
		crColor[7]=RGB(0,0,colorintensity);
		crColor[8]=RGB(colorintensity,0,colorintensity);*/

		//recreate the font to bold so that color annotation is more visible
		LOGFONT pLogFont;
	
		pLogFont.lfHeight = pDoc->FontSize;
		pLogFont.lfWidth = 0;
		pLogFont.lfEscapement = 0;
		pLogFont.lfOrientation = 0;
		pLogFont.lfWeight = FW_BOLD;
		pLogFont.lfItalic = 0;
		pLogFont.lfUnderline = 0;
		pLogFont.lfStrikeOut = 0;
		pLogFont.lfCharSet = ANSI_CHARSET;
		pLogFont.lfOutPrecision = OUT_DEFAULT_PRECIS;
		pLogFont.lfClipPrecision = CLIP_DEFAULT_PRECIS;
		pLogFont.lfQuality = PROOF_QUALITY;
		pLogFont.lfPitchAndFamily = VARIABLE_PITCH|FF_MODERN;
		strcpy(pLogFont.lfFaceName,"Courier New");

		delete Font;
		Font = new CFont;
		Font->CreateFontIndirect(&pLogFont);

	
		Invalidate();



	}
	delete filedialog;



}



void CDrawView::OnDrawShowshapeannotationkey()
{
	pDoc->app->DisplaySHAPEKeyWindow();
}

void CDrawView::OnDrawExporttodot()
{
	//Export to .bracket notation.

	CFileDialog *filedialog;


	filedialog = new CFileDialog(FALSE,".dbn",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Dot-bracket file (*.dbn)|*.dbn|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		pDoc->ct.writedotbracket(filedialog->GetPathName().GetBuffer(0));

	}
	delete filedialog;


}