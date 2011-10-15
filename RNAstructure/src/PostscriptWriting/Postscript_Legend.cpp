/*
 * An implementation file for a class that writes a Postscript color coded
 * legend.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#include "Postscript_Legend.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Legend::Postscript_Legend( int number ) {

  // Initialize the vector of colors.
  colors.resize( number );

  // Initialize the vector of texts.
  texts.resize( number );

}

///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Legend::~Postscript_Legend() {}

///////////////////////////////////////////////////////////////////////////////
// Get the colors used in the legend.
///////////////////////////////////////////////////////////////////////////////
vector<string> Postscript_Legend::getColors() {

  return colors;

}

///////////////////////////////////////////////////////////////////////////////
// Set the colors used in the legend.
// The base class version of this is a stub.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Legend::setColors() {}

///////////////////////////////////////////////////////////////////////////////
// Set the location of the legend on the Postscript page.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Legend::setLocation( int x, int y ) {

  locationX = x;
  locationY = y;

}

///////////////////////////////////////////////////////////////////////////////
// Set the texts used in the legend.
// The base class version of this is a stub.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Legend::setTexts() {}

///////////////////////////////////////////////////////////////////////////////
// Write the legend.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Legend::writeLegend( ostream &out ) {

  // Set the colors and texts appropriately.
  setColors();
  setTexts();

  // If the number of colors and texts aren't the same, set error and return.
  if( texts.size() != colors.size() ) {
    error = true;
    return;
  }

  // Write the legend array.
  out << "% Create the legend array." << endl
      << "/legend [" << endl;

  int size = texts.size();
  for( int i = 1; i <= size; i++ ) {
    out << tab << "[(" << texts.at( i - 1 ) << ") "
	<< colors.at( i - 1 ) << "]" << endl;
  }

  out << "] def" << endl << endl;

  // Write the loop that loads and shows the legend array.
  out << "% Show the legend." << endl
      << "/Courier-Bold findfont 8 scalefont setfont" << endl
      << locationX << " " << locationY << " moveto"
      << endl << endl
      << "legend { gsave aload pop setrgbcolor show grestore "
      << "0 -10 rmoveto } forall" << endl << endl;

}
