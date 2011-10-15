/*
 * An implementation file for a class that writes annotations on a Postscript
 * RNA secondary structure.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#include "Postscript_Annotation_Handler.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Annotation_Handler::Postscript_Annotation_Handler( int len,
							      int structs ) {

  // Initialize the length of the sequence, and the number of structures
  // under consideration.
  length = len;
  structures = structs;

}

///////////////////////////////////////////////////////////////////////////////
// Get the color of a base pair probability annotated nucleotide.
///////////////////////////////////////////////////////////////////////////////
string Postscript_Annotation_Handler::getProbabilityColor( int structure,
							   int nucleotide ) {

  // Get the annotation code.
  char value = probabilityAnnotations[structure - 1][nucleotide - 1];

  // Return the proper color, based on the annotation code.
  // If the annotation code is not recognized, the color black is returned,
  // resulting in an unannotated nucleotide.
  return
    ( value == 'a' ) ? RED :
    ( value == 'b' ) ? ORANGE :
    ( value == 'c' ) ? DARK_YELLOW :
    ( value == 'd' ) ? GREEN :
    ( value == 'e' ) ? BRIGHT_GREEN :
    ( value == 'f' ) ? LIGHT_BLUE :
    ( value == 'g' ) ? BLUE :
    ( value == 'h' ) ? DARK_PINK :
    ( value == 'i' ) ? BLACK :
    BLACK;

}

///////////////////////////////////////////////////////////////////////////////
// Get the color of a SHAPE-annotated nucleotide.
///////////////////////////////////////////////////////////////////////////////
string Postscript_Annotation_Handler::getSHAPEColor( int index ) {

  // Get the annotation code.
  char value = shapeAnnotations[index - 1];

  // Return the proper color, based on the annotation code.
  // If the annotation code is not recognized, the color black is returned,
  // resulting in an unannotated nucleotide.
  return
    ( value == 'a' ) ? RED :
    ( value == 'b' ) ? ORANGE :
    ( value == 'c' ) ? BLACK :
    ( value == 'd' ) ? GRAY :
    BLACK;

}

///////////////////////////////////////////////////////////////////////////////
// Read base pair probabilities from a partition function save file.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Annotation_Handler::readPartition( string file,
						   RNA* structureStrand ) {

  // Initialize the RNA strand and error checker that reads partition data.
  RNA* partStrand = new RNA( file.c_str(), PFS_TYPE );
  ErrorChecker<RNA>* partChecker = new ErrorChecker<RNA>( partStrand );

  // If the RNA strand and error checker were created successfully, read in the
  // annotation data.
  if( !( error = partChecker->isErrorStatus() ) ) {

    // If there are no structures in the strand, print out an error message.
    // Otherwise, initialize the annotation array to handle the appropriate
    // amount of structures.
    if( structures == 0 ) {
      cerr << "No structures or pairs are present to annotate." << endl;
      error = true;
    } else {
      probabilityAnnotations.resize( structures );
      for( int i = 1; i <= structures; i++ ) {
	vector<char> row;
	row.resize( length );
	for( int j = 1; j <= length; j++ ) { row[j-1] = 'i'; }
	probabilityAnnotations[i - 1] = row;
      }
    }

    // For each structure possible, read in its base pair probability data.
    for( int i = 1; i <= structures; i++ ) {

      // If an error has occurred, stop reading data.
      if( error ) { break; }

      // Loop through the structure to find pairs.
      for( int j = 1; j <= length; j++ ) {

	// If an error has occurred, stop reading data.
	if( error ) { break; }

	// Get the next pair. If an error occurred, stop reading data.
	int pair = structureStrand->GetPair( j, i );
	int code = structureStrand->GetErrorCode();
	if( code != 0 ) {
	  cerr << endl << structureStrand->GetErrorMessage( code ) << endl;
	  error = true;
	  break;
	}

	// If the next nucleotide is in fact paired, determine the proper
	// color code for it.
	if( ( pair != 0 ) && ( pair > j ) ) {

	  // Get the probability for this pair.
	  // If an error occurred, stop reading data.
	  double bp = partStrand->GetPairProbability( j, pair );
	  if( ( error = partChecker->isErrorStatus() ) ) { break; }

	  // Set the proper values for the color code.
	  probabilityAnnotations[i-1][j-1] =
	    ( bp >= 0.99 ) ? 'a' :
	    ( bp > 0.95 ) ? 'b' :
	    ( bp > 0.90 ) ? 'c' :
	    ( bp > 0.80 ) ? 'd' :
	    ( bp > 0.70 ) ? 'e' :
	    ( bp > 0.60 ) ? 'f' :
	    ( bp > 0.50 ) ? 'g' :
	    'h';

	  probabilityAnnotations[i-1][pair-1] =
	    probabilityAnnotations[i-1][j-1];
	}
      }
    }
  }

  // If an error occurred, print out an extra error message to make sure the
  // user knows the error came from reading the partition function annotation
  // file in.
  if( error ) {
    cerr << "Partition function save file not read successfully." << endl;
  }

  // Delete the RNA strand and error checker when they're no longer needed.
  delete partStrand;
  delete partChecker;

}

///////////////////////////////////////////////////////////////////////////////
// Read a SHAPE data file.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Annotation_Handler::readSHAPE( string file ) {

  // Clear the SHAPE annotations vector and resize it to the proper length.
  shapeAnnotations.clear();
  shapeAnnotations.resize( length );

  // Initialize all codes in the annotations vector to gray.
  for( int i = 1; i <= length; i++ ) { shapeAnnotations[i - 1] = 'd'; }

  // Open the file stream for the SHAPE data file.
  ifstream in( file.c_str() );
  string line;

  // If the file stream was opened successfully, begin reading data.
  bool error = false;
  if( in ) {

    // Read the data from the SHAPE file.
    while( !in.eof() ) {

      // Read the next line.
      getline( in, line );

      // Create a stringstream to extract integer values from the line.
      stringstream lineStream( stringstream::in | stringstream::out );
      lineStream << line;

      // Get the index.
      int index = 0;
      lineStream >> index;

      // Check to make sure the given index isn't out of range. If it is, set
      // the error flag and stop reading the file.
      if( index <= 0 || index > length ) {
	if( !in.eof() ) { error = true; }
	break;
      }

      // Get the SHAPE value.
      double shape = 0.0;
      lineStream >> shape;

      // Determine the color in this index in the annotations vector.
      char colorCode =
	( shape > 0.7 ) ? 'a' :
	( shape > 0.3 ) ? 'b' :
	( shape > -500 ) ? 'c' :
	'd';

      // Put the color code in the annotations vector.
      shapeAnnotations[index - 1] = colorCode;
    }

    // Close the file stream.
    in.close();
  }

  // Otherwise, set the error flag.
  else { error = true; }

  // If an error occurred, show an error message.
  if( error ) {
    cerr << "Error reading SHAPE data annotation file." << endl;
  }

}
