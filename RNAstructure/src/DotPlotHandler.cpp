/*
 * An implementation file for a class that holds a dot plot data structure and
 * handles all dot plot manipulation. This class can be used throughout
 * RNAstructure, for any type of dot plot.
 *
 * Copyright 2010, Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "DotPlotHandler.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
DotPlotHandler::DotPlotHandler( string inFile, string outFile ) {

  // Set the input and output file names.
  input = inFile;
  output = outFile;

  // Set boolean flags that track handler progress.
  error = false;
  dataRead = false;

  // Initialize the plot type to undefined, since it hasn't been specified yet.
  plotType = TYPE_UNDEFINED;

  // Initialize the default number of thresholds (entries) in the plot, as will
  // be shown on the plot's legend, if drawn.
  entries = DEFAULT_ENTRIES;

  // Initialize the default minimum and maximum.
  defaultMin = numeric_limits<double>::infinity();
  defaultMax = defaultMin * -1;

  // Initialize the current minimum and maximum.
  minimum = defaultMin;
  maximum = defaultMax;

}

///////////////////////////////////////////////////////////////////////////////
// Destructor.
///////////////////////////////////////////////////////////////////////////////
DotPlotHandler::~DotPlotHandler() {}

///////////////////////////////////////////////////////////////////////////////
// Build and/or rebuild the array of cutoff values for dots.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::buildThresholds() {

  // Clear the numeric values threshold array and size it properly.
  thresholds.clear();
  thresholds.resize( entries + 1 );

  // Calculate the separation between each value in the numeric values
  // threshold array, and place the minimum and maximum values in the array.
  double valueIncrement = ( maximum - minimum ) / entries;
  thresholds[0] = minimum;
  thresholds[entries] = maximum;

  // Place all the other increments that aren't minimum or maximum in the
  // numeric values threshold array.
  for( int i = 1; i < entries; i++ ) {
    thresholds[i] = minimum + ( valueIncrement * i );
  }

  // Clear and resize the color thresholds array as appropriate.
  colorInfo.clear();
  colorInfo.resize( entries );

  // Determine how many elements in the color thresholds array go from red to
  // green, and how many go from green to blue.
  bool odd = ( ( entries % 2 ) == 1 );
  int spaces = ( odd ) ? ( entries - 3 ) / 2 : ( entries - 2 ) / 2;
  double colorIncrement = 1.0 / (double)( spaces + 1 );
  int changePoint = ( ( odd ) ? spaces + 1 : spaces );

  // Place the color increments in the color thresholds array, working toward
  // green from both directions.
  for( int i = 0; i <= changePoint; i++ ) {

    // Create a vector to hold the next red-to-green row, and another to hold
    // the next green-to-blue row.
    vector<double> redGreenRow( 3 );
    vector<double> greenBlueRow( 3 );

    // Determine the percentage of either red or blue relative to the
    // percentage of green in this threshold.
    double next = colorIncrement * (double)i;
    double reciprocal = 1.0 - next;

    // Build the next red-to-green row.
    redGreenRow[0] = reciprocal;
    redGreenRow[1] = next;
    redGreenRow[2] = 0.0;
    colorInfo[i] = redGreenRow;

    // Build the text green-to-blue row.
    greenBlueRow[0] = 0.0;
    greenBlueRow[1] = next;
    greenBlueRow[2] = reciprocal;
    colorInfo[entries-1-i] = greenBlueRow;
  }

}

///////////////////////////////////////////////////////////////////////////////
// Make sure the bounds read from a data file are valid.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::checkBounds() {

  // If the minimum value has not been previously set, or is outside the
  // acceptable range, set it to the default minimum.
  if( !( minimum >= defaultMin && minimum <= defaultMax ) ) {
    minimum = defaultMin;
  }

  // If the maximum value has not been previously set, or is outside the
  // acceptable range, set it to the default maximum.
  if( !( maximum >= defaultMin && maximum <= defaultMax ) ) {
    maximum = defaultMax;
  }

}

///////////////////////////////////////////////////////////////////////////////
// Get a color string from the color thresholds array coded as Postscript.
///////////////////////////////////////////////////////////////////////////////
string DotPlotHandler::getColorAsPostscript( const int i, const int j ) {

  // Get the index of the proper color for the dot.
  int index = getThresholdIndex( i, j );

  // Create the color string based on the index in the color thresholds array.
  stringstream color;
  color << colorInfo[index][0] << " "
	<< colorInfo[index][1] << " "
	<< colorInfo[index][2];
  return color.str();

}

///////////////////////////////////////////////////////////////////////////////
// Get a color string from the color thresholds array coded as RGB.
///////////////////////////////////////////////////////////////////////////////
string DotPlotHandler::getColorAsRGB( const int i, const int j ) {

  // Get the index of the proper color for the dot.
  int index = getThresholdIndex( i, j );

  // Create the color string based on the index in the color thresholds array.
  stringstream color;
  color << (int)( colorInfo[index][0] * 255.0 ) << " "
	<< (int)( colorInfo[index][1] * 255.0 ) << " "
	<< (int)( colorInfo[index][2] * 255.0 );
  return color.str();

}

///////////////////////////////////////////////////////////////////////////////
// Get the number of entries/thresholds used with the plot.
///////////////////////////////////////////////////////////////////////////////
int DotPlotHandler::getEntries() {

  return entries;

}

///////////////////////////////////////////////////////////////////////////////
// Get the information string describing the plot.
///////////////////////////////////////////////////////////////////////////////
string DotPlotHandler::getInfoString() {

  return infoString;

}

///////////////////////////////////////////////////////////////////////////////
// Get the name of the input file that plot data comes from.
///////////////////////////////////////////////////////////////////////////////
string DotPlotHandler::getInputFile() {

  return input;

}

///////////////////////////////////////////////////////////////////////////////
// Get the maximum value in the plot.
///////////////////////////////////////////////////////////////////////////////
double DotPlotHandler::getMaximum() {

  return maximum;

}

///////////////////////////////////////////////////////////////////////////////
// Get the minimum value in the plot.
///////////////////////////////////////////////////////////////////////////////
double DotPlotHandler::getMinimum() {

  return minimum;

}

///////////////////////////////////////////////////////////////////////////////
// Get the name of the output file the plot is written to.
///////////////////////////////////////////////////////////////////////////////
string DotPlotHandler::getOutputFile() {

  return output;

}

///////////////////////////////////////////////////////////////////////////////
// Get the character code for the plot type.
///////////////////////////////////////////////////////////////////////////////
char DotPlotHandler::getPlotType() {

  return plotType;

}

///////////////////////////////////////////////////////////////////////////////
// Get the index in the numeric values threshold array of a particular dot.
///////////////////////////////////////////////////////////////////////////////
int DotPlotHandler::getThresholdIndex( const int i, const int j ) {

  // Get the value of the dot.
  double value = getValue( i, j );

  // Loop through the thresholds array until the right index is found.
  for( int i = 1; i <= entries; i++ ) {
    if( thresholds[i] > value ) { return i - 1; }
  }
  return entries - 1;

}

///////////////////////////////////////////////////////////////////////////////
// Get the value of a particular dot.
///////////////////////////////////////////////////////////////////////////////
double DotPlotHandler::getValue( const int i, const int j ) {

  return values[i-1][j-1];

}

///////////////////////////////////////////////////////////////////////////////
// Initialize the 2D vector of dots with default values.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::initializeValuesArray( int length ) {

  // Based on the length of the sequence, create a length x length matrix whose
  // cell values are all infinity.
  double infinity = numeric_limits<double>::infinity();
  for( int i = 1; i <= length; i++ ) {
    vector<double> row;
    for( int j = 1; j <= length; j++ ) { row.push_back( infinity ); }
    values.push_back( row );
  }

}

///////////////////////////////////////////////////////////////////////////////
// Get the error state of the dot plot handler.
///////////////////////////////////////////////////////////////////////////////
bool DotPlotHandler::isError() {

  return error;

}

///////////////////////////////////////////////////////////////////////////////
// Read data from a Dynalign save file, based on one specific sequence in it.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::readDynalignData( int sequence ) {

  // Create a Dynalign object and an error checker for it.
  Dynalign_object* dynalign = new Dynalign_object( input.c_str() );
  ErrorChecker<Dynalign_object>* checker =
    new ErrorChecker<Dynalign_object>( dynalign );

  // If the Dynalign object and error checker were created correctly, continue.
  if( !checker->isErrorStatus() ) {

    // Determine the length of the sequence under consideration, either the
    // first or second.
    int length = 0;

    if( sequence == 1 ) {
      length = dynalign->GetRNA1()->GetSequenceLength();
    } else {
      length = dynalign->GetRNA2()->GetSequenceLength();
    }

    // If an error occurred getting the length, set the error flag and return.
    // If the length was found sucessfully, initialize the 2D dot value vector
    // with default values.
    if( checker->isErrorStatus() ) {
      error = true;
      return;
    } else {
      initializeValuesArray( length );
    }

    // Loop through all the possible dots and set their values.
    for( int i = 1; i <= length; i++ ) {
      for( int j = i; j <= length; j++ ) {

	// Get the raw energy value of the dot.
	double value = dynalign->GetBestPairEnergy( i, j, sequence );

	// Disregard any energy values over 0.0; change it to infinity so it
	// isn't considered.
        if( value > 0.0 ) {
          value = numeric_limits<double>::infinity();
        }

	// If the energy value is valid, check to see if it's a new minimum or
	// maximum to set bounds with.
        if( value != numeric_limits<double>::infinity() ) {
          if( value > defaultMax ) { defaultMax = value; }
          else if( value < defaultMin ) { defaultMin = value; }
        }

	// Put the energy value in its proper place in the vector.
        values[i-1][j-1] = value;

	// If an error occurs, stop reading in data and return.
        if( checker->isErrorStatus() ) {
          error = true;
          return;
        }
      }
    }

    // Check the bounds of the read file for validity.
    checkBounds();
  }

  // If the Dynalign object and the error checker weren't contructed correctly,
  // set the error flag.
  else { error = true; }

  // If no error ocurred, set the info string by getting the comment from
  // either first or second sequence in the Dynalign save file, as necessary.
  if( !error ) {
    stringstream infoStream;

    infoStream << input << ", Sequence " << sequence << ": ";
    if( sequence == 1 ) {
      infoStream << dynalign->GetRNA1()->GetCommentString();
    } else {
      infoStream << dynalign->GetRNA2()->GetCommentString();
    }

    infoString = infoStream.str();
  }

  // Delete the error checker and the Dynalign object.
  delete checker;
  delete dynalign;

  // If no error occurred, build the color and value thresholds arrays.
  if( !isError() ) { buildThresholds(); }

  // If thresholds were set correctly, set the data read flag to true.
  if( !isError() ) { dataRead = true; }

}

///////////////////////////////////////////////////////////////////////////////
// Read in the dot data for Dynalign sequence 1.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::readDynalignSeq1Data() {

  plotType = TYPE_DYNALIGN1;
  readDynalignData( 1 );

}

///////////////////////////////////////////////////////////////////////////////
// Read in the dot data for Dynalign sequence 2.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::readDynalignSeq2Data() {

  plotType = TYPE_DYNALIGN2;
  readDynalignData( 2 );

}

///////////////////////////////////////////////////////////////////////////////
// Read in the dot data for a folding free energy plot.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::readFoldingData() {

  plotType = TYPE_ENERGY;
  label = "kcal/mol";
  readSingleStrandData( 4 );

}

///////////////////////////////////////////////////////////////////////////////
// Read in the dot data for a base pairing probability plot.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::readPartitionData() {

  plotType = TYPE_PROBABILITY;
  label = "-log10(BP Probability)";
  readSingleStrandData( 3 );

}

///////////////////////////////////////////////////////////////////////////////
// Read in the data for a plot derived from a single strand.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::readSingleStrandData( int type ) {

  // Check to see if the strand type is either partition function (3) or
  // folding (4) save file; if not, set the error flag and return.
  if( !( ( type == 3 ) || ( type == 4 ) ) ) {
    error = true;
    return;
  }

  // Set the info string as the input file.
  infoString = input;

  // Create the strand that reads the data, and the error checker for it.
  RNA* strand = new RNA( input.c_str(), type );
  ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );

  // If the strand and error checker were created successfully, get the
  // sequence length and initialize the dots vector.
  if( !checker->isErrorStatus() ) {
    int length = strand->GetSequenceLength();
    if( checker->isErrorStatus() ) {
      error = true;
      return;
    } else {
      initializeValuesArray( length );
    }

    // Loop through the possible pairings to fill the dots vector.
    for( int i = 1; i <= length; i++ ) {
      for( int j = i; j <= length; j++ ) {

	// Get the value for a pair between nucleotides i and j, depending on
	// the type of plot that's being made.
	double value = ( type == 3 ) ?
	  -log10( strand->GetPairProbability( i, j ) ):
	  strand->GetPairEnergy( i, j );

	// If the plot is a folding free energy plot, transform values above
	// 0.0 to infinity, so they aren't considered as dots.
	if( ( type == 4 ) && ( value > 0.0 ) ) {
	  value = numeric_limits<double>::infinity();
	}

	// If the value of the dot is not infinity, check to see if it is a new
	// minimum or maximum.
	if( value != numeric_limits<double>::infinity() ) {
	  if( value > defaultMax ) { defaultMax = value; }
	  else if( value < defaultMin ) { defaultMin = value; }
	}

	// Set the dot value in the array.
	values[i-1][j-1] = value;

	// If an error occurred, set the error flag and return.
	if( checker->isErrorStatus() ) {
	  error = true;
	  return;
	}
      }
    }

    // Check the bounds of the read file for validity.
    checkBounds();
  }

  // If the strand and error checker weren't created successfully, set the
  // error flag.
  else { error = true; }

  // Delete the strand and error checker.
  delete checker;
  delete strand;

  // If no error has occurred, set the color and value thresholds for the plot.
  if( !isError() ) { buildThresholds(); }

  // If thresholds were set successfully, set the data flag to true.
  if( !isError() ) { dataRead = true; }

}

///////////////////////////////////////////////////////////////////////////////
// Read dot plot data from a text file.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::readTextData() {

  // Open the input stream for the text file.
  string line;
  ifstream inFile( input.c_str() );

  // If the file was opened correctly, continue.
  if( inFile.is_open() ) {

    // Get the sequence length from the first line.
    // If the sequence length is valid, initialize the dots vector. If not,
    // set the error flag.
    getline( inFile, line );
    stringstream lengthStream( line );
    int length = 0;
    if( lengthStream >> length ) { initializeValuesArray( length ); }
    else { error = true; }

    // Determine the type of dot plot using the label on the second line.
    // Note that Dynalign dot plots are both energy plots and have the
    // same plot label, so these cannot be explicitly distinguished as
    // of yet and are both called energy plots.
    getline( inFile, line );
    stringstream labelStream( line );
    string label;
    for( int i = 1; i <= 3; i++ ) { labelStream >> label; }

    if( label == "kcal/mol" ) { plotType = TYPE_ENERGY; }
    else { plotType = TYPE_PROBABILITY; }

    // While the file has lines and no error has occurred, read data.
    stringstream dotStream( stringstream::in | stringstream::out );
    while( !inFile.eof() && !isError() ) {

      // Get the next line in the file.
      getline( inFile, line );
      dotStream << line;

      // Read the ith nucleotide. If i is not a valid number, set the
      // error flag and stop reading data.
      int ith = 0;
      if( dotStream >> ith ) { ith--; }
      else {
	error = true;
	break;
      }

      // Read the jth nucleotide. If j is not a valid number, set the
      // error flag and stop reading data.
      int jth = 0;
      if( dotStream >> jth ) { jth--; }
      else {
        error = true;
        break;
      }

      // Read the dot value for the pair between i and j.
      // If the value is valid, place it in the dots vector and check to see
      // if it's a new maximum or minimum.
      // If the value isn't valid, set the error flag and stop reading data.
      double dotVal = 0.0;
      if( dotStream >> dotVal ) {
	values[ith][jth] = dotVal;

	if( dotVal > defaultMax ) { defaultMax = dotVal; }
	else if( dotVal < defaultMin ) { defaultMin = dotVal; }
      } else {
	error = true;
	break;
      }
    }

    // Check the bounds of the read file for validity.
    checkBounds();

    // If no error occurred in reading, set the data read flag to true.
    dataRead = true;

    // Close the input stream.
    inFile.close();
  }

  // If the file was not opened correctly, set the error flag.
  else { error = true; }

}

///////////////////////////////////////////////////////////////////////////////
// Reset the range of the plot to the default min and max.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::resetRange() {

  minimum = defaultMin;
  maximum = defaultMax;

}

///////////////////////////////////////////////////////////////////////////////
// Set the number of entries in the dot plot.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::setEntries( int newEntries ) {

  if( ( newEntries >= MIN_COLORS ) && ( newEntries <= MAX_COLORS ) ) {
    entries = newEntries;
    buildThresholds();
  }

}

///////////////////////////////////////////////////////////////////////////////
// Set the maximum value in the dot plot.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::setMaximum( double newMaximum ) {

  if( ( ( newMaximum <= defaultMax ) && ( newMaximum >= minimum ) ) ||
      ( dataRead == false ) ) {
    maximum = newMaximum;
    buildThresholds();
  }

}

///////////////////////////////////////////////////////////////////////////////
// Set the minimum value in the dot plot.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::setMinimum( double newMinimum ) {

  if( ( ( newMinimum >= defaultMin ) && ( newMinimum <= maximum ) ) ||
      ( dataRead == false ) ) {
    minimum = newMinimum;
    buildThresholds();
  }

}

///////////////////////////////////////////////////////////////////////////////
// Set the dot plot type.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::setPlotType( char type ) {

  // Determine if the given type is valid.
  bool okType =
    ( type == TYPE_DYNALIGN1 ) ||
    ( type == TYPE_DYNALIGN2 ) ||
    ( type == TYPE_ENERGY ) ||
    ( type == TYPE_PROBABILITY );

  // If the given type is valid, set it as the plot type. If not, set the plot
  // type as undefined.
  if( okType ) {
    plotType = type;
  } else {
    type == TYPE_UNDEFINED;
  }

}

///////////////////////////////////////////////////////////////////////////////
// Write a dot plot text file.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::writePlotFile() {

  // Open the output stream to the text file and write the header.
  ofstream out( output.c_str() );
  int length = values.size();
  out << length << endl
      << "i\tj\t" << label << endl;

  // Loop through each possible pair/dot.
  for( int i = 1; i <= length; i++ ) {
    for( int j = 1; j <= length; j++ ) {

      // Get the value for a dot at i and j.
      double value = getValue( i, j );

      // If the dot value is within the acceptable range, write it into the
      // text file.
      if( ( value >= minimum ) && ( value <= maximum ) ) {
	out << i << "\t" << j << "\t" << value << endl;
      }
    }
  }

  // Close the written text file.
  out.close();

}
