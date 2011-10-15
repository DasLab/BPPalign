/*
 * A program that draws a structure and writes image output to a Postscript
 * image file.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "DrawStructure.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
DrawStructure::DrawStructure() {

  // Initialize the calculation type description.
  calcType = "Structure drawing";

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool DrawStructure::parse( int argc, char** argv ) {

  // Initialize arrays of flags without and with parameters.
  const char* flagsNoParams[] = {
    "-c", "-C", "--circle",
  };

  const char* flagsParams[] = {
    "-p", "-P", "--probability",
    "-s", "-S", "--SHAPE"
  };

  // Initialize usage string.
  string usageString =
    "USAGE: draw <ct file> <ps file> [options]";

  // Initialize required parameters string.
  string ctHelp = "\
<ct file>\n\
     The name of a file containing CT data for the structure to be drawn.\
";

  string psHelp = "\
<ps file>\n\
     The name of a Postscript image file to which output will be written.\
";

  string paramsString =
    ctHelp + "\n\n" + psHelp;

  // Initialize flags without parameters string.
  string circleHelp = "\
-c, -C, --circle\n\
     Specify that the structure should be drawn with its backbone stretched\n\
     around a circle.\n\
     Default is to show a collapsed structure.\
";

  string helpLabel = "\
-h, -H, --help\n\
     Display the usage details message.\
";

  string flagsNoParamsString =
    circleHelp + "\n\n" + helpLabel;

  // Initialize the flags with parameters string.
  string probHelp = "\
-p, -P, --probability\n\
     Specify the name of the file from which base pairing probability data\n\
     will be read for annotation.\n\
     This file must be a partition function save file.\n\
     Default is no probability annotation file used.\
";

  string shapeHelp = "\
-s, -S, --SHAPE\n\
     Specify the name of the file from which SHAPE data will be read for\n\
     annotation.\n\
     Default is no SHAPE annotation file used.\
";

  string flagsParamsString =
    probHelp + "\n\n" + shapeHelp;

  // Initialize parser and check the command line.
  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 2,
						   3, flagsNoParams,
						   6, flagsParams );
  parser->setUsageStrings( usageString.c_str(), paramsString.c_str(),
			   flagsNoParamsString.c_str(),
			   flagsParamsString.c_str() );
  bool noError = parser->checkLine();

  // If command line structure is OK, proceed with argument retrieval, parsing.
  if( noError ) {

    // Get required parameters.
    inputFile = parser->getParameter( 1 );
    outputFile = parser->getParameter( 2 );

    // Set non-parameterized flag(s).
    circular = parser->containsInGroup( "-c -C --circle" );

    // Set parameterized options
    // If one of these options is not found on the command line, the action to
    // set that option is skipped.
    // If one throws an error, an error message is shown and parsing continues,
    // so all errors are shown if there are more than one. 
    bool flagSet = true;

    // Probability annotation file option
    const char* probabilityString = probabilityFile.c_str();
    flagSet = parser->setOption( probabilityString, "-p -P --probability",
				 "Probability annotation file not found.",
				 true );
    probabilityFile = probabilityString;
    if( flagSet == false ) { noError = false; }

    // SHAPE annotation file option
    const char* SHAPEstring = SHAPEFile.c_str();
    flagSet = parser->setOption( SHAPEstring, "-s -S --SHAPE",
				 "SHAPE annotation file not found.", true );
    SHAPEFile = SHAPEstring;
    if( flagSet == false ) { noError = false; }
  }

  // Delete the parser and return the error state.
  delete parser;
  return noError;

}

///////////////////////////////////////////////////////////////////////////////
// Run structure drawing.
///////////////////////////////////////////////////////////////////////////////
void DrawStructure::run() {

  // Show a message that structure drawing has begun.
  cout << "Drawing structure(s)..." << flush;

  // Initialize the Postscript wrapper for the structure drawing, which is
  // the same whatever the number of type of structure drawings done.
  Postscript_Wrapper* wrap = new Postscript_Wrapper();
  bool noError = true;

  // If the user didn't want to have any unannotated output, write an
  // unannotated output file. If the user requested a circular structure
  // drawing, the drawing will be a circular structure. If not, the drawing
  // will be a radial structure.
  if( probabilityFile == "" && SHAPEFile == "" ) {

    if( circular ) {
      noError = wrap->structureCircular( inputFile, outputFile );
    } else {
      noError = wrap->structureRadial( inputFile, outputFile );
    }
  }


  // If the user wants probability annotated structures, write a probability
  // annotated output file. If the user requested a circular structure drawing,
  // the drawing will be a circular structure. If not, the drawing will be a
  // radial structure.
  if( probabilityFile != "" ) {

    if( circular ) {
      noError = wrap->structureCircular_Probability( inputFile, outputFile,
						     probabilityFile );
    } else {
      noError = wrap->structureRadial_Probability( inputFile, outputFile,
						   probabilityFile );
    }
  }

  // If the user wants SHAPE annotated structures, write a SHAPE annotated
  // output file. If the user requested a circular structure drawing, the
  // drawing will be a circular structure. If not, the drawing will be a radial
  // structure.
  if( SHAPEFile != "" ) {

    if( circular ) {
      noError = wrap->structureCircular_SHAPE( inputFile, outputFile,
					       SHAPEFile );
    } else {
      noError = wrap->structureRadial_SHAPE( inputFile, outputFile,
					     SHAPEFile );
    }
  }

  // Delete the Postscript wrapper for the structure drawing, which is the same
  // whatever the number of type of structure drawings done.           
  delete wrap;

  // Print confirmation of run finishing
  if( noError ) {
    cout << "done." << endl
	 << calcType << " complete." << endl;
  } else {
    cerr << calcType << " complete with errors." << endl;
  }

}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

  DrawStructure* runner = new DrawStructure();
  bool parseable = runner->parse( argc, argv );
  if( parseable == true ) { runner->run(); }
  delete runner;
  return 0;

}
