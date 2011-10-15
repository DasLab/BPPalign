/*
 * An implementation file for a class that writes a Postscript color coded
 * legend whose color increments change evenly from red (lowest values) to
 * green (median values) to blue (highest values), displayed with red first.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#include "Postscript_RGB_Legend.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Postscript_RGB_Legend::Postscript_RGB_Legend( int entries, string description )
  : Postscript_Legend( entries ) {

  divider = description;

}

///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
Postscript_RGB_Legend::~Postscript_RGB_Legend() {}

///////////////////////////////////////////////////////////////////////////////
// Set the lower and upper bounds of the legend.
///////////////////////////////////////////////////////////////////////////////
void Postscript_RGB_Legend::setBounds( double min, double max ) {

  minimum = min;
  maximum = max;

}

///////////////////////////////////////////////////////////////////////////////
// Set the colors for the legend.
///////////////////////////////////////////////////////////////////////////////
void Postscript_RGB_Legend::setColors() {

  // If the minimum and maximum bounds are the same, highlighting an exact
  // value rather than a range, set the color as just red, then return.
  if( minimum == maximum ) {
    colors.resize( 1 );
    colors[0] = RED;

    return;
  }

  // Get the size of the legend and set the lowest and highest colors, which
  // are always the same, red and blue.
  int size = colors.size();
  colors[0] = RED;
  colors[size - 1] = BLUE;

  // If the number of colors is odd, set the central color to green.
  if( size % 2 == 1 ) { colors[size / 2] = BRIGHT_GREEN; }

  // Determine the number of dynamic colors (red to green on the top half of
  // the legend, green to blue on the bottom half).
  int dynamics = ( size - 2 ) / 2;

  // Determine the legend color increment.
  double increment = 1.00 / ( (double)dynamics + 1.0 );

  // Set colors for the red-green section of the legend.
  for( int i = 1; i <= dynamics; i++ ) {
    double green = increment * (double)i;
    double red = 1.0 - green;

    stringstream stream( "" );
    stream << setprecision( 2 ) << red << " " << green << " " << 0.0;
    colors[i] = stream.str();
  }

  // Set colors for the green-blue section of the legend.
  int start = ( size % 2 == 0 ) ? dynamics + 1 : dynamics + 2;
  int end = size - 2;
  for( int i = start; i <= end; i++ ) {
    double blue = increment * (double)( i - ( start - 1 ) );
    double green = 1.0 - blue;

    stringstream stream( "" );
    stream << setprecision( 2 ) << 0.0 << " " << green << " " << blue;
    colors[i] = stream.str();
  }

}

///////////////////////////////////////////////////////////////////////////////
// Set the segment texts for this legend.
///////////////////////////////////////////////////////////////////////////////
void Postscript_RGB_Legend::setTexts() {

  // If the minimum and maximum bounds are the same, highlighting an exact
  // value rather than a range, set a simple text descriptor, then return.
  if( minimum == maximum ) {
    texts.resize( 1 );

    stringstream descriptor;
    descriptor << divider << " = " << minimum; 
    texts[0] = descriptor.str();

    return;
  }

  // Determine the range of the legend.
  double range = maximum - minimum;
  if( ( minimum < 0.0 ) && ( maximum < 0.0 ) ) {
    range = ( minimum - maximum ) * -1;
  }

  // Determine the increment between legend segments.
  int size = texts.size();
  double increment = range / (double)size;

  // Determine the proper indenting for each segment.
  string::size_type maxFirst = string::npos;
  for( int i = 1; i <= size; i++ ) {
    double value = minimum + ( increment * (double)( i - 1 ) );
    if( i == 1 ) { value = minimum; }

    stringstream stream( "" );
    stream << value;
    string val = stream.str();
    if( ( maxFirst == string::npos ) || ( val.length() > maxFirst ) ) {
      maxFirst = val.length();
    }
  }

  // Generate the texts.
  string lessThanEqual = "<=";
  string secondDivider = " <  ";
  for( int i = 1; i <= size; i++ ) {
    double first = minimum + ( increment * (double)( i - 1 ) );
    double second = minimum + ( increment * (double)i );

    if( i == 1 ) { first = minimum; }
    if( i == size ) {
      second = maximum;
      secondDivider = " " + lessThanEqual + " ";
    }

    stringstream firstStream( "" );
    bool firstScientific =
      ( ( first != 0.0 ) && ( ( first > -0.01 ) && ( first < 0.01 ) ) );
    if( firstScientific ) { firstStream << scientific; }
    firstStream << setw( maxFirst ) << setfill( ' ' ) << first;
    if( firstScientific ) { maxFirst = firstStream.str().length(); }

    stringstream secondStream( "" );
    bool secondScientific =
      ( ( second != 0.0 ) && ( ( second > -0.01 ) && ( second < 0.01 ) ) );
    if( secondScientific ) { secondStream << scientific; }
    secondStream << second;

    texts[i-1] = firstStream.str() + " " + lessThanEqual + " " +
      divider + secondDivider + secondStream.str();
  }

}
