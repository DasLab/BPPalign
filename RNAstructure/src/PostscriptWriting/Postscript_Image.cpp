/*
 * An implementation file for a class that writes a Postscript image file.
 *
 * (c) 2010  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "Postscript_Image.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Image::Postscript_Image( string input, string output ) {

  filename = input;
  psFile = output;

  strand = 0;
  checker = 0;

  strandLength = 0;

}

///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Image::~Postscript_Image() {

  if( strand != 0 ) { delete strand; }
  if( checker != 0 ) { delete checker; }

}

///////////////////////////////////////////////////////////////////////////////
// Initialize important strand helper data pieces.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Image::initializeStrandHelpers() {

  // Initialize the error checker.
  checker = new ErrorChecker<RNA>( strand );

  // If there is something wrong with the strand, set the error flag.
  if( checker->isErrorStatus() ) { error = true; }

  // Otherwise, set the sequence length with the valid strand.
  else { strandLength = strand->GetSequenceLength(); }

}

///////////////////////////////////////////////////////////////////////////////
// Set the vector of additional data.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Image::setAdditionalData( vector<string> data ) {

  additional = data;

}

///////////////////////////////////////////////////////////////////////////////
// Write any additional data after the object.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Image::writeAdditionalData( ofstream &out ) {

  int size = additional.size();
  for( int i = 1; i <= size; i++ ) { out << additional.at( i - 1 ) << endl; }

}

///////////////////////////////////////////////////////////////////////////////
// Write the description string for the image.
// The base class version of this is a stub.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Image::writeDescriptor( ofstream &out ) {}

///////////////////////////////////////////////////////////////////////////////
// Write the Postscript output file.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Image::writeOutput( bool append ) {

  // If an error has not occurred, begin writing the Postscript file.
  if( !error ) {

    // Set the mode of the output stream.
    ios_base::openmode mode = ( append ) ? ios_base::app : ios_base::out;

    // Open the output file for writing and write the Postscript shebang.
    ofstream out( psFile.c_str(), mode );
    out << "%!" << endl << endl;

    // Write the main changeable Postscript variables into the file.
    writeVariables( out );

    // Write the image file name/description string.
    writeDescriptor( out );

    // Write the translation and scaling of the main image.
    out << "% Write the translation and scaling of the main image." << endl
	<< "gsave" << endl
	<< "scaleFactorX scaleFactorY scale" << endl
	<< "translateFactorX translateFactorY translate" << endl << endl;

    // Write the unique commands that create a particular type of image.
    writeSpecificImageType( out );

    // Remove scaling after the main image is written.
    out << "% Restore the original scaling to 100%." << endl
	<< "grestore" << endl << endl;

    // Write any additional data (descriptors, captions, values) in the normal
    // coordinate space after the object itself has been written.
    writeAdditionalData( out );

    // Write the show page command and close the output file.
    out << "showpage" << endl;
    out.close();
  }

}

///////////////////////////////////////////////////////////////////////////////
// Write specific commands to create a specific image type.
// The base class version of this is a stub.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Image::writeSpecificImageType( ofstream &out ) {}

///////////////////////////////////////////////////////////////////////////////
// Write variables that are specific to an image.
// The base class version of this is a stub.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Image::writeVariables( ofstream &out ) {}
