/*
 * An implementation file for a class that writes a Postscript RNA secondary
 * structure.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#include "Postscript_Structure.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Structure::Postscript_Structure( string in, string out ) :
  Postscript_Image( in, out ) {

  // Create the proper CT strand.
  strand = new RNA( in.c_str(), 1 );

  // Initialize the error checker and sequence length.
  initializeStrandHelpers();

  // Initialize the pair types vector.
  pairTypes.resize( 2 );

  // Inititalize the pairings vector.
  pairings.resize( strandLength );

  // Initialize the type of pairing lines to not use thick pairs.
  useThickPairs = false;

  // Initialize the structure number to be drawn.
  structureNumber = 1;

  // Initialize the colorization allowed to none at all.
  typeColorCode = 0;

  // Initialize annotator, and annotation legends, to a default null pointer.
  annotator = 0;
  probabilityLegend = 0;
  shapeLegend = 0;

  // Initialize both boolean annotation flags to false.
  isProbabilityAnnotated = false;
  isSHAPEAnnotated = false;

}

///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Structure::~Postscript_Structure() {

  // If an annotator was used, delete it.
  if( annotator != 0 ) { delete annotator; }

  // If a probability legend was created, delete it.
  if( probabilityLegend != 0 ) { delete probabilityLegend; }

  // If a SHAPE legend was created, delete it.
  if( shapeLegend != 0 ) { delete shapeLegend; }

}

///////////////////////////////////////////////////////////////////////////////
// Add base pair probability annotation to the structure.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::addProbabilityAnnotation( string file ) {

  // Set nucleotides able to be annotated.
  setNucleotidesAnnotated();

  // Set the flag that says the structure is probability annotated.
  isProbabilityAnnotated = true;

  // Get the number of structures in the strand as a reference for annotation.
  int structures = strand->GetStructureNumber();

  // Create the annotator.
  annotator = new Postscript_Annotation_Handler( strandLength, structures );

  // Read the base pair probability data into the annotator.
  annotator->readPartition( file, strand );

  // Create the probability annotation legend.
  probabilityLegend = new Postscript_Probability_Legend();
  placeLegend( probabilityLegend );

}

///////////////////////////////////////////////////////////////////////////////
// Add SHAPE annotation to the structure.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::addSHAPEAnnotation( string file ) {

  // Set nucleotides able to be annotated.
  setNucleotidesAnnotated();

  // Set the flag that says the structure is SHAPE annotated.
  isSHAPEAnnotated = true;

  // Create the annotator.
  annotator = new Postscript_Annotation_Handler( strandLength );

  // Read the SHAPE data into the annotator.
  annotator->readSHAPE( file );

  // Create the SHAPE annotation legend.
  shapeLegend = new Postscript_SHAPE_Legend();
  placeLegend( shapeLegend );

}

///////////////////////////////////////////////////////////////////////////////
// Calculate the number of pseudoknotted pairs in the structure.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::calculatePseudoknottedPairs() {

  // Break pseudoknots in the strand.
  int errorCode = strand->BreakPseudoknot( false );
  if( ( error = checker->isErrorStatus( errorCode ) ) ) { return; }

  // If no error occurred, count the number of pairs in the strand without
  // pseudoknots.
  int pairs = 0;
  for( int i = 1; i <= strandLength; i++ ) {
    if( strand->GetPair( i, structureNumber ) > i ) { pairs++; }
  }

  // Subtract the number of pairs from the strand without pseudoknots from the
  // number of pairs in the strand with pseudoknots to get the number of
  // pseudoknotted pairs.
  pairTypes[1] = pairTypes[0] - pairs;

}

///////////////////////////////////////////////////////////////////////////////
// Get the number of pairs in the structure.
///////////////////////////////////////////////////////////////////////////////
int Postscript_Structure::getNumPairs() {

  return pairTypes[0];

}

///////////////////////////////////////////////////////////////////////////////
// Get the number of pseudoknotted pairs in the structure.
///////////////////////////////////////////////////////////////////////////////
int Postscript_Structure::getNumPseudoknottedPairs() {

  return pairTypes[1];

}

///////////////////////////////////////////////////////////////////////////////
// Get a particular pair from the structure.
///////////////////////////////////////////////////////////////////////////////
int Postscript_Structure::getPair( int pair ) {

  return pairings[pair-1];

}

///////////////////////////////////////////////////////////////////////////////
// Get the RNA strand backing this structure.
///////////////////////////////////////////////////////////////////////////////
RNA* Postscript_Structure::getStrand() {

  return strand;

}

///////////////////////////////////////////////////////////////////////////////
// Save pairs from the structure.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::parsePairs( int structure ) {

  // If an error has not occurred, begin parsing pairs in the given structure.
  if( !error ) {

    // Initialize the pairing array to its proper length.
    pairings.clear();
    pairings.resize( strandLength );

    // Initialize the pairing type array to its proper length.
    pairTypes.clear();
    pairTypes.resize( 2 );

    // For each nucleotide in the sequence, check if it's paired.
    for( int i = 1; i <= strandLength; i++ ) {

      // Check if the nucleotide is paired, then check for an error.
      int pair = strand->GetPair( i, structure );
      if( ( error = checker->isErrorStatus() ) ) { break; }

      // If no error occurred, check to make sure the pair hasn't already been
      // registered; if it has, don't re-register it in the pairings array.
      pairings[i-1] = 0;
      if( pair != 0 && pair > i ) {
	pairings[i-1] = pair;
	pairTypes[0]++;
      }
    }

    // Set the structure number to be drawn as the structure that was parsed.
    structureNumber = structure;
  }

}

///////////////////////////////////////////////////////////////////////////////
// Set the location of an annotation legend for this structure, if applicable.
// The base class version of this method is a stub.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::placeLegend( Postscript_Legend* legend ) {}

///////////////////////////////////////////////////////////////////////////////
// Set the type color code so no color annotation is allowed.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::setNoAnnotation() {

  // Set the color code to no annotation.
  typeColorCode = 0;

  // Set both boolean annotation flags to false.
  isProbabilityAnnotated = false;
  isSHAPEAnnotated = false;

}

///////////////////////////////////////////////////////////////////////////////
// Set the type color code so nucleotides can be color annotated if desired.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::setNucleotidesAnnotated() {

  typeColorCode = 1;

}

///////////////////////////////////////////////////////////////////////////////
// Set the type color code so pairs can be color annotated if desired.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::setPairingsAnnotated() {

  typeColorCode = 2;

}

///////////////////////////////////////////////////////////////////////////////
// Set the type color code so both pairs and nucleotides can be color annotated
// if desired.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::setPairingsAndNucleotidesAnnotated() {

  typeColorCode = 3;

}

///////////////////////////////////////////////////////////////////////////////
// Write all structures in a file, not just one.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::writeAllStructures() {

  // Get the number of structures in the file.
  int structures = strand->GetStructureNumber();

  // For each structure in the file, parse the pairs and write the file.
  for( int i = 1; i <= structures; i++ ) {
    bool append = ( i != 1 );

    parsePairs( i );
    writeOutput( append );
  }
}

///////////////////////////////////////////////////////////////////////////////
// Write the backbone of the structure.
// The base class version of this is a stub.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::writeBackbone( ofstream &out ) {}

///////////////////////////////////////////////////////////////////////////////
// Write the array of base annotation colors.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::writeBaseAnnotations( ofstream &out ) {

  // Only write a base annotations array if one of the annotation flags is set.
  if( ( isProbabilityAnnotated == true ) || ( isSHAPEAnnotated == true ) ) {

    // Write the opening of the annotation array.
    out << "% Write the array of annotation colors." << endl
	<< "/annotations [" << endl;

    // If the structure is base pair probability annotated, build its array.
    if( isProbabilityAnnotated ) {
      for( int i = 1; i <= strandLength; i++ ) {
	out << tab << "["
	    << annotator->getProbabilityColor( structureNumber, i )
	    << "]" << endl;
      }
    }

    // Otherwise, if the structure is SHAPE annotated, build its array.
    else {
      for( int i = 1; i <= strandLength; i++ ) {
	out << tab << "[" << annotator->getSHAPEColor( i ) << "]" << endl;
      }
    }

    // Write the closing of the annotation array.
    out << "] def" << endl << endl;

    // Write the annotation legend.
    if( isProbabilityAnnotated ) { probabilityLegend->writeLegend( out ); }
    else if( isSHAPEAnnotated ) { shapeLegend->writeLegend( out ); }
  }

}

///////////////////////////////////////////////////////////////////////////////
// Write the array of bases as a Postscript variable.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::writeBasesArray( ofstream &out ) {

  // Write the array of nucleotides.
  out << "% Create the array of nucleotides." << endl
      << "/bases [" << endl << tab << flush;

  for( int i = 1; i <= strandLength; i++ ) {
    out << "(" << strand->GetNucleotide( i ) << ") " << flush;
  }

  out << endl << "] def" << endl << endl;

}

///////////////////////////////////////////////////////////////////////////////
// Write the descriptor for the structure. In the base class, the descriptor is
// the annotation data, if necessary.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::writeDescriptor( ofstream &out ) {

  // If necessary, write the base annotation data.
  if( !error ) { writeBaseAnnotations( out ); }

} 

///////////////////////////////////////////////////////////////////////////////
// Write the nucleotides on the backbone, and nucleotide labels as necessary.
// The base class version of this is a stub.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::writeNucleotidesAndLabels( ofstream &out ) {}

///////////////////////////////////////////////////////////////////////////////
// Write structure-dependent Postscript that handles how a pair is drawn.
// The base class version of this is a stub.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::writePair( ofstream &out ) {}

///////////////////////////////////////////////////////////////////////////////
// Write the pairings in the structure.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::writePairings( ofstream &out ) {

  // Write the opening of the loop that writes pairings.
  out << "% Write the pairing lines inside the backbone." << endl;

  if( typeColorCode != 0 ) { out << "0 setgray" << endl; }
  if( useThickPairs ) { out << "3 setlinewidth" << endl; }

  out << "0 1 numPairings {" << endl
      << tab << "/pair pairings currentPairing get def" << endl << endl;

  // Write the x,y coordinates of the first nucleotide in the pair.
  out << tab << "% Determine coordinates of first nucleotide in pair." << endl
      << tab << "/base1 pair 0 get 1 sub def" << endl
      << tab << "/point1 basePoints base1 get def" << endl
      << tab << "/x1 point1 0 get def" << endl
      << tab << "/y1 point1 1 get def" << endl << endl;

  // Write the x,y coordinates of the second nucleotide in the pair.
  out << tab << "% Determine coordinates of second nucleotide in pair." << endl
      << tab << "/base2 pair 1 get 1 sub def" << endl
      << tab << "/point2 basePoints base2 get def" << endl
      << tab << "/x2 point2 0 get def" << endl
      << tab << "/y2 point2 1 get def" << endl << endl;

  // If the pair is supposed to be colored, set its appropriate color.
  if( typeColorCode >= 2 ) {
    out << tab << "% If the pair should be colored, set its appropriate color."
	<< endl
	<< tab << "pair 2 get pair 3 get pair 4 get setrgbcolor"
	<< endl << endl;
  }

  // Write drawing of the new pair.
  out << tab << "% Draw current pair, then increment current pairing." << endl;
  writePair( out );

  // If a pair color was set, reset color to black.
  if( typeColorCode >= 2 ) { out << tab << "0 setgray" << endl; }

  // Increment the current pair.
  out << tab << "/currentPairing currentPairing 1 add def" << endl;

  // Write the closing of the loop that writes pairings.
  out << "} repeat" << endl << endl;

}

///////////////////////////////////////////////////////////////////////////////
// Write the array of pairings as a Postscript variable.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::writePairingsArray( ofstream &out ) {

  // Write the array of pairings.
  out << "% Write the array of pairings." << endl
      << "/pairings [" << endl;

  for( int i = 0; i < strandLength; i++ ) {
    if( pairings[i] != 0 ) {
      int nucleotide = i + 1;
      int pair = pairings[i];
      out << tab << "[" << nucleotide << " " << pair << flush;
      if( typeColorCode >= 2 ) { out << " -1" << flush; }
      out << "]" << endl;
    }
  }

  out << "] def" << endl << endl;

}

///////////////////////////////////////////////////////////////////////////////
// Write specific variables pertaining to pairings.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::writePairingsVariables( ofstream &out ) {

  // Write variables concerning number and placement of pairings.
  out << "% Set variables handling number, placement of pairings." << endl
      << "/numPairings " << getNumPairs() << " def" << endl
      << "/numPseudoknotted " << getNumPseudoknottedPairs() << " def" << endl
      << "/currentPairing 0 def" << endl << endl;

}

///////////////////////////////////////////////////////////////////////////////
// Write specific Postscript commands to write a structure.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Structure::writeSpecificImageType( ofstream &out ) {

  // Write the backbone of the structure.
  if( !error ) { writeBackbone( out ); }

  // Write the pairings in the structure.
  if( !error ) { writePairings( out ); }

  // Write the nucleotides and their labels.
  if( !error ) { writeNucleotidesAndLabels( out ); }

}

