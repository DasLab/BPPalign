/*
 * An implementation file for a class that writes a radial Postscript
 * structure. The "radial" in the class's name refers to the underlying
 * algorithm; a structure is not displayed as a circle.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#include "Postscript_Radial_Structure.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Radial_Structure::Postscript_Radial_Structure( string file,
							  string out ) :
  Postscript_Structure( file, out ) {

  useThickPairs = true;

}

///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Radial_Structure::~Postscript_Radial_Structure() {}

///////////////////////////////////////////////////////////////////////////////
// Place an annotation legend for this structure.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Radial_Structure::placeLegend( Postscript_Legend* legend ) {

  legend->setLocation( 40, 710 );

}

///////////////////////////////////////////////////////////////////////////////
// Write out the backbone of the structure.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Radial_Structure::writeBackbone( ofstream &out ) {

  // Write opening of loop that writes the backbone.
  out << "% Use repeat loop to write the backbone." << endl
      << "0 1 numBases {" << endl
      << tab << "clear" << endl << endl;

  // Save the nucleotide location and show the nucleotide.
  // Also, if necessary, draw lines on which a future label will sit.
  out << tab << "% Save nucleotide location and write next backbone section."
      << endl
      << tab << "/current coordinates currentBase get def" << endl
      << tab << "basePoints currentBase current put" << endl
      << tab << "/base currentBase 1 add def" << endl
      << tab << "/x current 0 get def" << endl
      << tab << "/y current 1 get def" << endl << endl
      << tab << "currentBase 0 ne {" << endl
      << tab << tab << "/previous coordinates currentBase 1 sub get def"
      << endl
      << tab << tab << "/prevX previous 0 get def" << endl
      << tab << tab << "/prevY previous 1 get def" << endl
      << tab << tab << "newpath prevX prevY moveto x y lineto closepath stroke"
      << endl << endl
      << tab << tab
      << "% If necessary, write a line to the future nucleotide label." << endl
      << tab << tab << "base 10 mod 0 eq {" << endl
      << tab << tab << tab << "% Determine where the number label is." << endl
      << tab << tab << tab << "/index base 10 idiv 1 sub def" << endl
      << tab << tab << tab << "/labelCoords coordinateNumbers index get def"
      << endl
      << tab << tab << tab << "/labelX labelCoords 0 get def" << endl
      << tab << tab << tab << "/labelY labelCoords 1 get def"
      << endl << endl
      << tab << tab << tab << "% Write the line to the number label, if valid."
      << endl
      << tab << tab << tab << "labelX 0 ne labelY 0 ne and {" << endl
      << tab << tab << tab << tab
      << "newpath x y moveto labelX labelY lineto stroke" << endl
      << tab << tab << tab << "} if" << endl
      << tab << tab << "} if" << endl
      << tab << "} if" << endl << endl;

  // Write incrementing of current base.
  out << tab << "% Increment base." << endl
      << tab << "/currentBase currentBase 1 add def" << endl;

  // Write closing of loop that writes the backbone.
  out << "} repeat" << endl << endl;

}

///////////////////////////////////////////////////////////////////////////////
// Write the descriptor (file name) of the radial structure.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Radial_Structure::writeDescriptor( ofstream &out ) {

  // Write the separator for the plot writing section.
  out << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
      << "% Write out the radial structure." << endl
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << endl << endl;

  // Write the structure description string.
  out << "% Write the structure description string." << endl
      << "horizontalKeyLocation verticalDescriptionLocation moveto" << endl
      << "descriptionString show" << endl << endl;

  // Write annotation data, if necessary.
  Postscript_Structure::writeDescriptor( out );

}

///////////////////////////////////////////////////////////////////////////////
// Write nucleotides and number labels on the structure.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Radial_Structure::writeNucleotidesAndLabels( ofstream &out ) {

  // Write opening of loop that places nucleotides.
  out << "% Use repeat loop to write nucleotides on the backbone." << endl
      << "0.5 setlinewidth" << endl
      << "0 1 numBases {" << endl
      << tab << "clear" << endl << endl;

  // Write saving of the nucleotide location and number.
  out << tab << "% Save the nucleotide location and number." << endl
      << tab << "/current coordinates currentBase2 get def" << endl
      << tab << "/x current 0 get def" << endl
      << tab << "/y current 1 get def" << endl
      << tab << "/base currentBase2 1 add def"
      << endl << endl;

  // Write the white spacer box behind the nucleotide.
  out << tab << "% Write the white spacer box behind the nucleotide." << endl
      << tab << WHITE << " setrgbcolor" << endl
      << tab << "newpath" << endl
      << tab << "x quarterFont sub 1 sub y quarterFont sub 1 sub moveto"
      << endl
      << tab << "0 halfFont quarterFont add rlineto" << endl
      << tab << "halfFont quarterFont add 0 rlineto" << endl
      << tab << "0 halfFont quarterFont add -1 mul rlineto" << endl
      << tab << "closepath fill" << endl << endl;

  // Set the color of the nucleotide.
  if( isProbabilityAnnotated || isSHAPEAnnotated ) {
    out << tab << "% Set the color of the nucleotide based on the annotation."
	<< endl
	<< tab << "/color annotations currentBase2 get def" << endl
	<< tab << "color 0 get color 1 get color 2 get setrgbcolor"
	<< endl << endl;
  } else {
    out << tab << BLACK << " setrgbcolor" << endl << endl;
  }

  // Write placing of the nucleotide in the box.
  out << tab << "% Place the nucleotide in the box." << endl
      << tab << "x quarterFont sub y quarterFont sub moveto" << endl
      << tab << "bases currentBase2 get show"
      << endl << endl;

  // Write the beginning of the loop that writes a number label, as necessary.
  out << tab << "% If a condition is met for a number label, write a label."
      << endl
      << tab << "base 10 mod 0 eq {" << endl;

  // Determine the proper length of the number label.
  out << tab << tab << "% Determine proper length of number label." << endl
      << tab << tab << "/stringChars 2 def" << endl
      << tab << tab << "base 10 ge { /stringChars 3 def } if" << endl
      << tab << tab << "base 100 ge { /stringChars 4 def } if" << endl
      << tab << tab << "base 1000 ge { /stringChars 5 def } if"
      << endl << endl;

  // Write creation of the number label.
  out << tab << tab << "% Create the number label." << endl
      << tab << tab << "/numString stringChars string def" << endl
      << tab << tab << "base numString cvs"
      << endl << endl;

  // Determine where the number label should be.
  out << tab << tab << "% Determine where the number label should be." << endl
      << tab << tab << "/index base 10 idiv 1 sub def" << endl
      << tab << tab << "/labelCoords coordinateNumbers index get def" << endl
      << tab << tab << "/labelX labelCoords 0 get def" << endl
      << tab << tab << "/labelY labelCoords 1 get def"
      << endl << endl;

  // Write the beginning of the loop that places a valid label.
  out << tab << tab << "% Place the label if it's valid." << endl
      << tab << tab << "labelX 0 ne labelY 0 ne and {" << endl;

  // Write the white spacer box behind a valid label.
  out << tab << tab << tab << "% Write the white spacer box behind the label."
      << endl
      << tab << tab << tab << WHITE << " setrgbcolor" << endl
      << tab << tab << tab << "newpath" << endl
      << tab << tab << tab
      << "labelX quarterFont sub 1 sub labelY quarterFont sub 1 sub moveto"
      << endl
      << tab << tab << tab << "0 halfFont quarterFont add rlineto" << endl
      << tab << tab << tab
      << "halfFont quarterFont add stringChars 1 sub mul 0 rlineto" << endl
      << tab << tab << tab << "0 halfFont quarterFont add -1 mul rlineto"
      << endl
      << tab << tab << tab << "closepath fill"
      << endl << endl;

  // Place the label in the box.
  out << tab << tab << tab << "% Place the label in the box." << endl
      << tab << tab << tab << BLACK << " setrgbcolor" << endl
      << tab << tab << tab
      << "labelX quarterFont sub labelY quarterFont sub moveto" << endl
      << tab << tab << tab << "numString show"
      << endl;

  // Write the end of the loop that places a valid label.
  out << tab << tab << "} if" << endl;

  // Write the end of the loop that writes a number label, as necessary.
  out << tab << "} if" << endl << endl;

  // Write incrementing of current base.
  out << tab << "% Increment base." << endl
      << tab << "/currentBase2 currentBase2 1 add def" << endl;

  // Write closing of loop that writes nucleotides on the backbone.
  out << "} repeat" << endl << endl;

}

///////////////////////////////////////////////////////////////////////////////
// Write radial structure Postscript that handles how a pair is written.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Radial_Structure::writePair( ofstream &out ) {

  out << tab << "newpath x1 y1 moveto x2 y2 lineto closepath stroke" << endl;

}

///////////////////////////////////////////////////////////////////////////////
// Write the changeable variables for the radial structure.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Radial_Structure::writeVariables( ofstream &out ) {

  // Write the separator for the variables section.
  out << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
      << "% Set variables for structure " << structureNumber << "." << endl
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << endl << endl;

  // Write variables for font size and type.
  int fontSize = 12;
  if( structureNumber == 1 ) {
    out << "% Set font size and type." << endl
	<< "/fontSize " << fontSize << " def" << endl
	<< "/quarterFont fontSize 4 div def" << endl
	<< "/halfFont fontSize 2 div def" << endl
	<< "/Courier-Bold findfont fontSize scalefont setfont"
	<< endl << endl;
  }

  // Write variables for the main description string.
  out << "% Set variables for the main description string." << endl;

  if( structureNumber == 1 ) {
    out << "/horizontalKeyLocation 40 def" << endl
	<< "/verticalDescriptionLocation 725 def" << endl;
  }

  out << "/descriptionString ("
      << escapeAndTrim( strand->GetCommentString( structureNumber ).c_str(),
			TRUNC_END )
      << ") def"
      << endl << endl;

  // Write variables concerning the number and placement of nucleotides.
  out << "% Set variables handling number, placement of nucleotides." << endl;

  if( structureNumber == 1 ) {
      out << "/numBases " << strandLength << " def" << endl
	  << "/labelSpace fontSize def" << endl;
  }

  out << "/currentBase 0 def" << endl
      << "/currentBase2 0 def" << endl
      << "/basePoints numBases array def"
      << endl << endl;

  // Write variables concerning number and placement of pairings and bases.
  if( structureNumber == 1 ) { writeBasesArray( out ); }
  writePairingsArray( out );
  writePairingsVariables( out );

  // Check to make sure valid drawing coordinates can be found. If not, set the
  // error flag. If so, write the array of coordinates.
  int errorCode =
    strand->DetermineDrawingCoordinates( fontSize, fontSize, structureNumber );

  if( !( error = checker->isErrorStatus( errorCode ) ) ) {

    // Initialize the drawing coordinates vector.
    vector< vector<double> > coords;
    coords.resize( strandLength );

    for( int i = 1; i <= strandLength; i++ ) {
      vector<double> row;
      row.resize( 2 );
      if( i == 1 || i == strandLength || i % 10 == 0 ) { row.resize( 4 ); }
      coords[i-1] = row; 
    }

    // Fill the drawing coordinates array with the nucleotide coordinates.
    for( int i = 1; i <= strandLength; i++ ) {
      for( int j = 0; j <= 1; j++ ) {
	int number =
	  ( j == 0 ) ? strand->GetNucleotideXCoordinate( i ) :
	  strand->GetNucleotideYCoordinate( i );

	errorCode = strand->GetErrorCode();
	if( !( error = checker->isErrorStatus( errorCode ) ) ) {
	  coords[i-1][j] = number;
	} else {
	  i = strandLength + 1;
	  j = 2;
	} 
      }

      if( i % 10 == 0 ) {
	int labelCoordX = strand->GetLabelXCoordinate( i );
	if( checker->isErrorStatus() ) {
	  error = true;
	  break;
	}

	int labelCoordY = strand->GetLabelYCoordinate( i );
	if( checker->isErrorStatus() ) {
	  error = true;
	  break;
	}

	if( !isError() ) {
	  coords[i-1][2] = labelCoordX;
	  coords[i-1][3] = labelCoordY;
	}
      }
    }

    // If coordinates array was sucessfully filled, translate coordinates
    // as necessary, determine minimum, maximum bounds, and write the
    // coordinates array to a file.
    double maxX = 0, maxY = 0;
    if( !error ) {
      double minX = numeric_limits<double>::infinity();
      maxX = numeric_limits<double>::infinity() * -1;
      double minY = minX;
      maxY = maxX;

      // Translate coordinates so the least point is at (3,3) to make all   
      // coordinates positive                                               
      for( int i = 0; i < strandLength; i++ ) {
	double x = coords[i][0], y = coords[i][1], labelX = x, labelY = y;
	if( ( ( ( i + 1 ) % 10 ) == 0 ) &&
	    ( coords[i][2] != numeric_limits<int>::infinity() ) ) {
	  labelX = coords[i][2];
	  labelY = coords[i][3];
	}

	minX = min( minX, min( x, labelX ) );
	minY = min( minY, min( y, labelY ) );
	maxX = max( maxX, max( x, labelX ) );
	maxY = max( maxY, max( y, labelY ) );
      }

      minX -= 3;
      minY -= 3;

      for( int i = 1; i <= strandLength; i++ ) {
	if( minX < 0 ) { coords[i-1][0] -= minX; }
	if( minY < 0 ) { coords[i-1][1] -= minY; }

	if( i % 10 == 0 ) {
	  if( minX < 0 && coords[i-1][2] != 0 ) { coords[i-1][2] -= minX; }
	  if( minY < 0 && coords[i-1][3] != 0 ) { coords[i-1][3] -= minY; }
	}
      }

      if( minX < 0 ) { maxX -= minX; }
      if( minY < 0 ) { maxY -= minY; }

      // Write the coordinates array.
      out << "% Create the array of nucleotide coordinates." << endl
	  << "/coordinates [" << endl;

      for( int i = 1; i <= strandLength; i++ ) {
	out << tab << "[" << coords[i-1][0] << " " << coords[i-1][1] << "]"
	    << endl;
      }

      out << "] def" << endl << endl;

      // Write the labels array.
      out << "% Create the array of nucleotide number labels." << endl
	  << "/coordinateNumbers [" << endl;

      for( int i = 1; i <= strandLength; i++ ) {
	if( i % 10 == 0 ) {
	  out << tab << "[" << coords[i-1][2] << " " << coords[i-1][3] << "]"
	      << endl;
	}
      }

      out << "] def" << endl << endl;
    }

    // If an error did not occur, write variables handling scaling and
    // translation of the structure. In the case of the radial structure, they
    // directly depend on the bounds calculated in the coordinates array.
    if( !error ) {
      double scaleX = 512 / maxX;
      double scaleY = 590 / maxY;
      double scaleFactor = ( scaleX < scaleY ) ? scaleX : scaleY;

      out << "% Set variables handling scaling, translation of the structure."
	  << endl
	  << "/scaleFactor " << scaleFactor << " def" << endl
	  << "/scaleFactorX scaleFactor def" << endl
	  << "/scaleFactorY scaleFactor def" << endl
	  << "/translateFactorX 45 def" << endl
	  << "/translateFactorY 90 def"
	  << endl << endl;
    }
  }
  
}
