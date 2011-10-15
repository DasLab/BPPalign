/*
 * An implementation file for a class that writes a triangular Postscript dot
 * plot.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#include "Postscript_Dot_Plot.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Dot_Plot::Postscript_Dot_Plot( DotPlotHandler* handler )
  : Postscript_Image( handler->getInputFile(), handler->getOutputFile() ) {

  // Initialize the dot plot handler.
  plotHandler = handler;

}

///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Dot_Plot::~Postscript_Dot_Plot() {}

///////////////////////////////////////////////////////////////////////////////
// Get the dot plot handler that backs this image.
///////////////////////////////////////////////////////////////////////////////
DotPlotHandler* Postscript_Dot_Plot::getHandler() {

  return plotHandler;

}

///////////////////////////////////////////////////////////////////////////////
// Write the descriptors (file name and key) of the plot.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Dot_Plot::writeDescriptor( ofstream &out ) {

  // Write the separator for the plot writing section.
  out << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
      << "% Write out the dot plot." << endl
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl;

  // Write the plot key.
  writeSpecificKey( out );

  // Write the plot description string.
  out << "% Write the plot description string." << endl
      << "/Courier-Bold findfont 12 scalefont setfont" << endl
      << "36 72 moveto" << endl
      << "(" << escapeAndTrim( plotHandler->getInfoString(), TRUNC_BEGIN )
      << ") show"
      << endl << endl;

}

///////////////////////////////////////////////////////////////////////////////
// Write the dots on the dot plot.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Dot_Plot::writeDots( ofstream &out ) {

  // Get the plot minimum and maximum, and number of entries in its legend.
  double minimum = plotHandler->getMinimum();
  double maximum = plotHandler->getMaximum();
  int entries = plotHandler->getEntries();

  // Write the dot placement routine.
  out << "% Write the dot placement routine." << endl
      << "/dot {" << endl
      << tab << "/x exch 4 mul 1 sub def" << endl
      << tab << "/y exch 4 mul 1 sub def" << endl
      << tab << "newpath" << endl
      << tab << "/y seqLength y sub def" << endl
      << tab << "x y moveto" << endl
      << tab << "x y dotSize sub lineto" << endl
      << tab << "x dotSize add y dotSize sub lineto" << endl
      << tab << "x dotSize add y lineto" << endl
      << tab << "closepath" << endl
      << tab << "setrgbcolor" << endl
      << tab << "fill" << endl
      << "} def" << endl << endl;

  // Open the values array.
  out << "% Write the values array." << endl
      << "/values [" << endl;

  // If the values fall over a range, do all the values in the range.
  if( minimum != maximum ) {

    // Determine the dot range.
    double range = maximum - minimum;
    if( ( minimum < 0 ) && ( maximum < 0 ) ) { range = minimum - maximum; }
    double increment = abs( range / (double)plotHandler->getEntries() );

    // Write the dots.
    for( int i = 1; i <= strandLength; i++ ) {
      for( int j = 1; j <= strandLength; j++ ) {
	double value = plotHandler->getValue( i, j );

	if( ( value >= minimum ) && ( value <= maximum ) ) {
	  out << "  ["
	      << plotHandler->getColorAsPostscript( i, j ) << " "
	      << i << " " << j
	      << "]" << endl;
	}
      }
    }
  }

  // Otherwise, if only a specific value should be shown, pick it out.
  else {

    // Write dots that match the value asked for.
    for( int i = 1; i <= strandLength; i++ ) {
      for( int j = 1; j <= strandLength; j++ ) {
	double value = plotHandler->getValue( i, j );

	// Both the minimum and maximum are the specific value, search for min.
	if( value == minimum ) {
	  out << "  [1 0 0 " << i << " " << j << "]" << endl;
	}
      }
    }
  }

  // Close the dot array.
  out << "] def" << endl << endl;

  // Write the dot placement routine.
  out << "% Write the for loop that shows the dots." << endl
      << "values { aload pop dot } forall" << endl << endl;

}

///////////////////////////////////////////////////////////////////////////////
// Write the grid that is the backbone of the dot plot.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Dot_Plot::writeGrid( ofstream &out ) {

  // Write the triangle plot border.
  out << "% Write the triangle plot border." << endl
      << "newpath" << endl
      << "0 seqLength moveto" << endl
      << "seqLength seqLength lineto" << endl
      << "seqLength 0 lineto" << endl
      << "closepath" << endl
      << "1 setlinewidth" << endl
      << "stroke" << endl << endl;

  // Write the gridline routine.
  out << "% Write the gridline routine." << endl
      << "/gridline { newpath moveto lineto stroke } def" << endl << endl;

  // Write the scale number placement routine.
  out << "% Write the scale number placement routine." << endl
      << "/scaleNumber {" << endl
      << tab << "gsave" << endl
      << tab << "moveto" << endl
      << tab << "rotate" << endl
      << tab << "0 -5 rmoveto" << endl
      << tab << "/increment 12 string def" << endl
      << tab << "increment cvs" << endl
      << tab << "increment show" << endl
      << tab << "grestore" << endl
      << "} def" << endl << endl;

  // Write gridlines.
  out << "% Write evenly spaced gridlines." << endl
      << "/Courier-Bold findfont 5 scalefont setfont" << endl
      << "/counter 0 def" << endl
      << "/number 1 def" << endl
      << "divisions {" << endl
      << tab << "/x 40 counter mul def" << endl
      << tab << "/x2 seqLength x sub def" << endl
      << tab << "/y seqLength 20 add def" << endl
      << tab << "x y x y 20 sub x sub gridline" << endl
      << tab << "y seqLength x sub x x2 gridline" << endl
      << tab << "counter 0 ne { /number counter 10 mul def } if" << endl
      << tab << "number 90 x seqLength 1 add scaleNumber" <<  endl
      << tab << "number 0 seqLength 1 add x2 scaleNumber" << endl
      << tab << "/counter counter 1 add def" << endl
      << "} repeat"
      << endl << endl;

  // Write extra end gridline if sequence length isn't a multiple of 10.
  out << "% Write end gridline if seq length isn't divisible by 10." << endl
      << "seqLength 10 mod 0 ne {" << endl
      << tab << "seqLength y seqLength y 20 sub gridline" << endl
      << tab << "seqLength 4 idiv 90 seqLength seqLength 1 add scaleNumber"
      << endl
      << tab << "y 0 seqLength 0 gridline" << endl
      << tab << "seqLength 4 idiv 0 seqLength 1 add 0 scaleNumber" << endl
      << "} if"
      << endl << endl;

}

///////////////////////////////////////////////////////////////////////////////
// Write the key for the dot plot.
// The base class version of this is a stub and does nothing.
// This method MUST be overridden in subclasses.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Dot_Plot::writeSpecificKey( ofstream &out ) {}

///////////////////////////////////////////////////////////////////////////////
// Write the key for the dot plot. This is the unchanging workhorse method.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Dot_Plot::writeKey( Postscript_RGB_Legend* legend,
				    ofstream &out ) {

    // Set specifics of legend and write it.
  legend->setBounds( plotHandler->getMinimum(), plotHandler->getMaximum() );
  legend->setLocation( 36, 300 );
  legend->writeLegend( out );

  // Delete the legend.
  delete legend;

}

///////////////////////////////////////////////////////////////////////////////
// Write a specific sequence of commands that writes the plot.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Dot_Plot::writeSpecificImageType( ofstream &out ) {

  // Write the plot border and grid lines.
  if( !error ) { writeGrid( out ); }

  // Write the dots on the plot.
  if( !error ) { writeDots( out ); }

}

///////////////////////////////////////////////////////////////////////////////
// Write any necessary variables for the dot plot.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Dot_Plot::writeVariables( ofstream &out ) {

  // Write the header for the variables section.
  out << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << "%%%%%%%%%%%" << endl
      << "% Write the variables necessary to create the dot plot." << endl
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << "%%%%%%%%%%%"
      << endl << endl;

  // Define the adjusted sequence length, and by extension the plot scale
  // and translation factors.
  int adjustedLength = strandLength * 4;
  double scale = 540.0 / (double)( adjustedLength + 20 );

  double reciprocal = 1.0 / scale;
  double translateX = 36.0 * reciprocal;
  double translateY = 216.0 * reciprocal;

  // Write the sequence length, scale, and translation factors.
  out << "% Define sequence length, scale, and translation factors." << endl
      << "/seqLength " << adjustedLength << " def" << endl
      << "/scaleFactorX " << scale << " def" << endl
      << "/scaleFactorY " << scale << " def" << endl
      << "/translateFactorX " << translateX << " def" << endl
      << "/translateFactorY " << translateY << " def" << endl
      << endl;

  // Define the size of the dots.
  out << "% Define the size of the dots." << endl
      << "/dotSize " << 2 << " def" << endl << endl;

  // Determine the number of gridlines in each direction.
  int gridLines = ( adjustedLength / 40 ) + 1;

  out << "% Determine the number of grid lines in each direction." << endl
      << "/divisions " << gridLines << " def"
      << endl << endl;

}
