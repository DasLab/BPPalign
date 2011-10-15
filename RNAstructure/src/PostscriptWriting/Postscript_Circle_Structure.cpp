/*
 * An implementation file for a class that writes a circular Postscript
 * structure from varied input sources.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#include "Postscript_Circle_Structure.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Circle_Structure::Postscript_Circle_Structure( string in,
							  string out ) :
  Postscript_Structure( in, out ) {

}

///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Circle_Structure::~Postscript_Circle_Structure() {}

///////////////////////////////////////////////////////////////////////////////
// Set the location of a legend for this structure.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Circle_Structure::placeLegend( Postscript_Legend* legend ) {

  legend->setLocation( 40, 595 );

}

///////////////////////////////////////////////////////////////////////////////
// Write the nucleotides out as the circular backbone.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Circle_Structure::writeBackbone( ofstream &out ) {

  // Write opening of loop that writes nucleotides on the circular path.
  out << "% Use repeat loop to write nucleotides on the circular path." << endl
      << "0 1 numBases {" << endl
      << tab << "clear" << endl << endl;

  // Write commands that move to the appropriate point in the circle to place
  // the next nucleotide.
  out << tab << "% Move to appropriate point, angle to write next nucleotide."
      << endl
      << tab << "center center moveto" << endl
      << tab << "currentBase angle mul rotate" << endl
      << tab << "quarterFont radius rmoveto" << endl << endl;

  // Determine the x,y coordinates of the nucleotide location.
  out << tab << "% Determine x and y coordinates of the nucleotide location."
      << endl
      << tab << "currentBase angle mul -1 mul rotate" << endl
      << tab << "currentpoint" << endl
      << tab << "/y exch cvi def" << endl
      << tab << "/x exch cvi def" << endl
      << tab << "currentBase angle mul rotate" << endl
      << tab << "quarterFont -1 mul 0 rmoveto" << endl << endl;

  // Save the nucleotide location and show the nucleotide.
  out << tab << "% Save the nucleotide location and show the nucleotide."
      << endl
      << tab << "/point [ x y ] def" << endl
      << tab << "basePoints currentBase point put" << endl << endl;

  // Write variables describing conditions where number labels are needed.
  out << tab << "% Set variables for conditions where number labels are found."
      << endl
      << tab << "/numCond1 currentBase 1 add 10 mod 0 eq def" << endl
      << tab << "/numCond2 currentBase 1 add 1 eq def" << endl
      << tab << "/numCond3 currentBase 1 add numBases eq def" << endl << endl;

  // If the structure is annotated, set the proper base pair color.
  if( isProbabilityAnnotated || isSHAPEAnnotated ) {
    out << tab << "/annotation annotations currentBase get def" << endl
	<< tab
	<< "annotation 0 get annotation 1 get annotation 2 get setrgbcolor"
	<< endl << endl;
  }

  // Write the conditional statement that places number labels if one or more
  // conditions are met.
  out << tab << "% If a condition is met for a number label, write a label."
      << endl
      << tab << "numCond1 numCond2 or numCond3 or {" << endl
      << tab << tab << "/Courier-Bold findfont fontSize scalefont setfont"
      << endl
      << tab << tab << "bases currentBase get show" << endl
      << tab << tab << "halfFont -1 mul labelSpace rmoveto" << endl
      << tab << tab << "/numString 7 string def" << endl
      << tab << tab << "currentBase 1 add numString cvs show" << endl
      << tab << tab << "0 labelSpace -1 mul rmoveto" << endl
      << tab << tab << "/Courier findfont fontSize scalefont setfont" << endl
      << tab << "}" << endl
      << tab << "{ bases currentBase get show } ifelse" << endl;

  // If the structure is annotated, set the base color back to black.
  if( isProbabilityAnnotated || isSHAPEAnnotated ) {
    out << tab << "0 setgray" << endl;
  }
  out << endl;

  // Write movement back to the center of the circle, rotation back to the
  // original orientation, and incrementing of current base.
  out << tab << "% Return to circle center, rotate to 0, increment base."
      << endl
      << tab << "0 radius -1 mul rmoveto" << endl
      << tab << "currentBase angle mul -1 mul rotate" << endl
      << tab << "/currentBase currentBase 1 add def" << endl;

  // Write closing of loop that writes nucleotides on the circular path.
  out << "} repeat" << endl << endl;

}

///////////////////////////////////////////////////////////////////////////////
// Write the descriptor (file name) of the plot.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Circle_Structure::writeDescriptor( ofstream &out ) {

  // Write the variables for the main description/comment string.
  out << "% Set variables for the main description string." << endl
      << "/horizontalKeyLocation 40 def" << endl
      << "/verticalDescriptionLocation 725 def" << endl
      << "/descriptionString ("
      << escapeAndTrim( strand->GetCommentString( structureNumber ).c_str(),
			TRUNC_END )
      << ") def"
      << endl << endl;

  // Write the separator for the plot writing section.
  out << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
      << "% Write out the circular plot." << endl
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl;

  // Write the structure description string.
  out << "% Write the structure description string." << endl
      << "horizontalKeyLocation verticalDescriptionLocation moveto" << endl
      << "descriptionString show" << endl << endl;

  // Write annotation data, if necessary.
  Postscript_Structure::writeBaseAnnotations( out );

}

///////////////////////////////////////////////////////////////////////////////
// Write circular structure Postscript that handles how a pair is written.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Circle_Structure::writePair( ofstream &out ) {

  // Write the determination of the number of bases between those paired.
  out << tab << "/between base2 base1 sub def"
      << endl << tab
      << "between numBases 2 div gt { /between numBases between sub def } if"
      << endl << endl;

  // Write determination of the midpoint coordinates between the paired bases.
  out << tab << "/midX x1 x2 add 2 div def" << endl
      << tab << "/midY y1 y2 add 2 div def"
      << endl << endl;

  // Write the arbitrary constant values that handle curvature.
  out << tab << "/gamma 0.9 def" << endl
      << tab << "/centerThresh radius 2 mul 8 div 5 mul def"
      << endl << endl;

  // Write the distance between the paired nucleotides.
  out << tab
      << "/ends x2 x1 sub x2 x1 sub mul y2 y1 sub y2 y1 sub mul add sqrt def"
      << endl;

  // Write calculation of the appropriate X and Y distance for curves.
  out << tab << "/distance between 2 mul numBases div radius mul gamma mul def"
      << endl
      << tab << "/lineAngle center midY sub center midX sub atan def" << endl
      << tab << "/distX lineAngle cos distance mul 2 mul def" << endl
      << tab << "/distY lineAngle sin distance mul 2 mul def"
      << endl << endl;

  // Write the control points.
  out << tab << "/controlX midX distX add def" << endl
      << tab << "/controlY midY distY add def"
      << endl << endl;

  // Write editing of the control points as an if statement (if necessary).
  out << tab << "ends centerThresh ge {" << endl
      << tab << tab << "/controlX center def" << endl
      << tab << tab << "/controlY center def" << endl
      << tab << "} if"
      << endl << endl;

  // Write the actual pairing as a curve.
  out << tab << "x1 y1 moveto x1 y1 controlX controlY x2 y2 curveto stroke"
      << endl;

}

///////////////////////////////////////////////////////////////////////////////
// Write the changeable variables for the circular plot.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Circle_Structure::writeVariables( ofstream &out ) {

  // Write the separator for the variables section.
  out << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
      << "% Set variables that change with each circular plot." << endl
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl;

  // Write variables for font size and type.
  int fontSize = 24;
  out << "% Set font size and type." << endl
      << "/fontSize " << fontSize << " def" << endl
      << "/quarterFont fontSize 4 div def" << endl
      << "/halfFont fontSize 2 div def" << endl
      << "/Courier findfont fontSize scalefont setfont" << endl << endl;

  // Write variables concerning the number and placement of nucleotides.
  out << "% Set variables handling number, placement of nucleotides." << endl
      << "/numBases " << strandLength << " def" << endl
      << "/currentBase 0 def" << endl
      << "/basePoints numBases array def" << endl << endl;

  // Write variables concerning position and scaling of the circular backbone.
  int radius = ( 3 * strandLength ) + fontSize;
  double scale = (double)( 270 - fontSize ) / (double)radius;

  int extraBases =
    ( strandLength >= 10000 ) ? 4 :
    ( strandLength >= 1000 ) ? 3 :
    ( strandLength >= 100 ) ? 2 :
    1;

  // Write variables handling scaling and translation of circular backbone.
  out << "% Set variables handling scaling, translation of circular backbone."
      << endl
      << "/scaleFactor " << scale << " def" << endl
      << "/scaleFactorX scaleFactor def" << endl
      << "/scaleFactorY scaleFactor def" << endl
      << "/translateFactorX 0 def" << endl
      << "/translateFactorY 0 def"
      << endl << endl;

  // Write variables handling properties of circular backbone.
  out << "% Set variables handling properties of circular backbone." << endl
      << "/radius " << radius << " def" << endl
      << "/center 306 scaleFactor div def" << endl
      << "/angle 360 numBases " << extraBases << " add div -1 mul def" << endl
      << "/labelSpace 20 def"
      << endl << endl;

  // Write variables concerning number and placement of pairings and bases.
  writeBasesArray( out );
  writePairingsArray( out );
  writePairingsVariables( out );

}
