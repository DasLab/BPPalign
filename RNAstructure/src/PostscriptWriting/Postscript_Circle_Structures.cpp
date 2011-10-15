/*
 * An implementation file for a class that writes a single circular Postscript
 * structure whose pairs show a comparison between two distinct structures.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#include "Postscript_Circle_Structures.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Circle_Structures::Postscript_Circle_Structures( string in,
							    string out )
  : Postscript_Circle_Structure( in, out ) {

  // Initialize the second pairings vector.
  pairings2.resize( 0 );

  // Initialize the number of combined unique pairs to 0.
  uniquePairs = 0;

  // Set the pairings annotated.
  setPairingsAnnotated();

  // Set the default structure IDs.
  structureID1 = "Structure 1";
  structureID2 = "Structure 2";

}

///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Circle_Structures::~Postscript_Circle_Structures() {}

///////////////////////////////////////////////////////////////////////////////
// Create the data structures for the second structure.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Circle_Structures::addSecondStructure( string second,
						       int number ) {

  // Save any relevant data from the first strand.
  // However, only do this if the data hasn't been saved already.
  if( pairings2.size() == 0 ) {
    filename2 = filename;
    pairings2.resize( strandLength );

    for( int i = 1; i <= strandLength; i++ ) {
      pairings2[i-1] = pairings[i-1];
    }

    comment2 = strand->GetCommentString();
  }

  // Delete the old RNA strand and error checker.
  delete strand;
  delete checker;

  // Initialize the second strand.
  strand = new RNA( second.c_str(), CT_TYPE );

  // Re-initialize the strand helpers.
  initializeStrandHelpers();

  // Set the second structure filename.
  filename = second;

  // Parse the appropriate structure number to fill data of the new strand.
  parsePairs( number );

}

///////////////////////////////////////////////////////////////////////////////
// Get the number of combined pairs between the structures.
///////////////////////////////////////////////////////////////////////////////
int Postscript_Circle_Structures::getNumPairs() {

  return uniquePairs;

}

///////////////////////////////////////////////////////////////////////////////
// Set the colors which identify pairing types.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Circle_Structures::setColors( const string& both,
					      const string& bothN,
					      const string& one,
					      const string& oneN,
					      const string& two,
					      const string& twoN ) {

  // Set color and name for pairs in both structures.
  bothColor = both;
  bothName = bothN + ":";

  // Set color and name for pairs in structure 1.
  oneColor = one;
  oneName = oneN + ":";

  // Set color and name for pairs in structure 2.
  twoColor = two;
  twoName = twoN + ":";

  // Add padding so color names are all the same length.
  size_t maxLength = max( bothName.length(),
			  max( oneName.length(), twoName.length() ) );

  while( bothName.length() != maxLength ) { bothName.append( " " ); }
  while( oneName.length() != maxLength ) { oneName.append( " " ); }
  while( twoName.length() != maxLength ) { twoName.append( " " ); }

}

///////////////////////////////////////////////////////////////////////////////
// Set new identifiers for structures 1 and 2.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Circle_Structures::setIdentifiers( const string id1,
						   const string id2 ) {

  structureID1 = id1;
  structureID2 = id2;


}

///////////////////////////////////////////////////////////////////////////////
// Write the combined structure file names and color key.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Circle_Structures::writeDescriptor( ofstream &out ) {

  // Create and pad the file name descriptor strings for use later.
  string fileDesc1 = structureID2 + " file name: ";
  string fileDesc2 = structureID1 + " file name: ";

  size_t maxFileDescLength = max( fileDesc1.length(), fileDesc2.length() );
  while( fileDesc1.length() != maxFileDescLength ) { fileDesc1.append( " " ); }
  while( fileDesc2.length() != maxFileDescLength ) { fileDesc2.append( " " ); }

  // Create and pad new structure ID strings for use later.
  string padID1 = structureID2 + ": ";
  string padID2 = structureID1 + ": ";

  size_t maxIDLength = max( padID1.length(), padID2.length() );
  while( padID1.length() != maxIDLength ) { padID1.append( " " ); }
  while( padID2.length() != maxIDLength ) { padID2.append( " " ); }

  // Write the header for the combined structure section.
  out << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << "%%%%%%%%%%%" << endl
      << "%% Write out the combined circular structure." << endl
      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      << "%%%%%%%%%%%" << endl << endl;

  // Write the font size and type for the descriptor.
  out << "% Set the font size and type for the descriptor." << endl
      << "/Courier-Bold findfont fontSize 2 div scalefont setfont"
      << endl << endl;

  // Create the beginning of the descriptor array.
  out << "/numDescriptors 7 def" << endl
      << "/descriptors [" << endl;

  // Write the input file names in the array.
  out << tab << "[(" << fileDesc1 << escapeAndTrim( filename, TRUNC_BEGIN )
      << ") " << BLACK << " 40 740]"
      << endl
      << tab << "[(" << fileDesc2 << escapeAndTrim( filename2, TRUNC_BEGIN )
      << ") " << BLACK << " 40 725]"
      << endl;

  // Write the structure comments in the array.
  out << tab << "[(" << padID1
      << escapeAndTrim( strand->GetCommentString( structureNumber ),
			TRUNC_END )
      << ") "
      << BLACK << " 40 695]" << endl
      << tab << "[(" << padID2
      << escapeAndTrim( comment2, TRUNC_END ) << ") " << BLACK << " 40 680]"
      << endl;

  // Write the color key values in the array.
  out << tab << "[(" << bothName << " Pair in both structures) " << bothColor
      << " 40 650]" << endl
      << tab << "[(" << oneName << " Pair in " << structureID1 << " only) "
      << oneColor << " 40 635]" << endl
      << tab << "[(" << twoName << " Pair in " << structureID2 << " only) "
      << twoColor << " 40 620]" << endl;

  // Write the end of the descriptor array.
  out << "] def" << endl << endl;

  // Write the contents of the descriptor array as Postscript.
  out << "% Write the contents of the descriptor array." << endl
      << "descriptors { aload pop moveto setrgbcolor show } forall"
      << endl << endl;

  // Write annotation data, if necessary.
  Postscript_Structure::writeBaseAnnotations( out );

  // Reset the font back to its original size and type.
  out << "% Reset the font to its original size, type, and color." << endl
      << "/Courier findfont fontSize scalefont setfont" << endl
      << "0 setgray" << endl << endl;

}

///////////////////////////////////////////////////////////////////////////////
// Write the array of pairings as a Postscript variable.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Circle_Structures::writePairingsArray( ofstream &out ) {

  // Reset the number of unique pairs to 0 beforevwriting this array, because
  // the number of pairs depends on how this array is written.
  uniquePairs = 0;

  // Write the array of pairings.
  out << "% Write the array of pairings." << endl
      << "/pairings [" << endl;

  for( int i = 1; i <= strandLength; i++ ) {
    int pair1 = pairings2[i-1];
    int pair2 = pairings[i-1];

    bool onePair = ( pair1 != 0 );
    bool twoPair = ( pair2 != 0 );
    bool pairEqual = ( pair1 == pair2 );

    if( onePair || twoPair ) {
      if( onePair ) {
	uniquePairs++;

	if( pairEqual ) {
	  out << tab << "[" << i << " " << pair1 << " " << bothColor << "]"
	      << endl;
	} else {
	  out << tab << "[" << i << " " << pair1 << " " << oneColor << "]"
	      << endl;
	}
      }

      if( twoPair && !pairEqual ) {
	uniquePairs++;

	out << tab << "[" << i << " " << pair2 << " " << twoColor << "]"
	    << endl;
      }
    }
  }

  out << "] def" << endl << endl;

}
