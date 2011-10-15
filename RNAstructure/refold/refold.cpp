/*
 * A program that folds a strand of nucleic acids from a previously folded save
 * file.
 * These nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "refold.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
//////////////////////////////////////////////////////////////////////////////
refold::refold() {

  // Initialize the calculation type description.
  calcType = "Refolding";

  // Initialize the maximum percent energy difference.
  percent = 10;

  // Initialize the maximum number of structures.
  maxStructures = 20;

  // Initialize the folding window size.
  windowSize = -1;

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool refold::parse( int argc, char** argv ) {

  // Initialize array of flags with parameters.
  const char* flagsParams[] = {
    "-m", "-M", "--maximum",
    "-p", "-P", "--percent",
    "-w", "-W", "--window",
  };

  // Initialize usage string.
  string usageString =
    "USAGE: refold <save file> <ct file> [options]";

  // Initialize required parameters string.
  string saveHelp = "\
<save file>\n\
     The name of a folding save file that contains structure data.\
";

  string ctHelp = "\
<ct file>\n\
     The name of a CT file to which output will be written.\
";

  string paramsString =
    saveHelp + "\n\n" + ctHelp;

  // Initialize flags without parameters string.
  string dnaHelp = "\
-d, -D, --DNA\n\
     Specify that the sequence is DNA, and DNA parameters are to be used.\n\
     Default is to use RNA parameters.\
";

  string helpLabel = "\
-h, -H, --help\n\
     Display the usage details message.\
";

  string flagsNoParamsString =
    dnaHelp + "\n\n" + helpLabel;

  // Initialize flags with parameters string.
  string maxStructuresHelp = "\
-m, -M, --maximum\n\
     Specify a maximum number of structures.\n\
     Default is 20 structures.\
";

  string percentHelp = "\
-p, -P, --percent\n\
     Specify a maximum percent energy difference.\n\
     Default is 10 percent (specified as 10, not 0.1).\
";

  string windowHelp = "\
-w, -W, --window\n\
     Specify a window size.\n\
     Default is determined by the length of the sequence.\
";

  string flagsParamsString =
    maxStructuresHelp + "\n\n" + percentHelp + "\n\n" + windowHelp;

  // Initialize parser and check the command line.
  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 2,
						   0, NULL,
						   9, flagsParams );
  parser->setUsageStrings( usageString.c_str(), paramsString.c_str(),
                           flagsNoParamsString.c_str(),
			   flagsParamsString.c_str() );
  bool noError = parser->checkLine();

  // If command line structure is OK, proceed with argument retrieval, parsing.
  if( noError ) {

    // Get required parameters.
    saveFile = parser->getParameter( 1 );
    ctFile = parser->getParameter( 2 );

    // Set parameterized options
    // If one of these options is not found on the command line, the action to
    // set that option is skipped.
    // If one throws an error, an error message is shown and parsing continues,
    // so all errors are shown if there are more than one.
    bool flagSet = true;

    // Maximum number of structures option
    flagSet = parser->setOption( maxStructures, "-m -M --maximum",
                                 "Invalid maximum structures given." "> 0" );
    if( flagSet == false ) { noError = false; }

    // Percent difference option
    flagSet = parser->setOption( percent, "-p -P --percent",
                                 "Invalid percent difference given.", ">= 0" );
    if( flagSet == false ) { noError = false; }

    // Window size option
    flagSet = parser->setOption( windowSize, "-w -W --window",
                                 "Invalid window size given.", ">= 0" );
    if( flagSet == false ) { noError = false; }
  }

  // Delete the parser and return the error state.
  delete parser;
  return noError;

}

///////////////////////////////////////////////////////////////////////////////
// Run refolding calculations.
///////////////////////////////////////////////////////////////////////////////
void refold::run() {

  // Create a variable to handle errors.
  int error = 0;

  /*
   * Use the constructor for RNA that specifies a filename.
   * Specify type = 4 (folding save file).
   * The save file handles the type of nucleic acid being folded, so it can be
   * set as the default (RNA) in the constructor.
   *
   * After construction of the strand data structure, create the error checker
   * which monitors for errors.  
   * Throughout, the error status of the calculation is checked with a variant
   * of the isErrorStatus method, which returns 0 if no error occurred. The
   * calculation proceeds as long as error = 0.
   */
  cout << "Initializing nucleic acids..." << flush;
  RNA* strand = new RNA( saveFile.c_str(), 4 );
  ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
  error = checker->isErrorStatus();
  if( error == 0 ) { cout << "done." << endl; }

  /*
   * Set the window size, based on the length of the sequence given as input.
   * Only do this, though, if window size hasn't been set on the command line.
   * Use method GetSequenceLength to get the length of the sequence.
   * The window sizes in relation to the length are hardcoded values.
   */
  if( windowSize == -1 && error == 0 ) {
    int length = strand->GetSequenceLength();
    windowSize = ( length > 1200 ) ? 20 :
      ( length > 800 ) ? 15 :
      ( length > 500 ) ? 11 :
      ( length > 300 ) ? 7 :
      ( length > 120 ) ? 5 :
      ( length > 50 ) ? 3 : 2;
  }

  /*
   * Refold the regenerated strand using the ReFoldSingleStrand method.
   * After the main calculation is complete, use the error checker's
   * isErrorStatus method to check for errors.
   */
  if( error == 0 ) {

    // Show a message that the main calculation has started.
    cout << "Refolding nucleic acids..." << flush;

    // Refold strand and check for errors.
    int mainCalcError =
        strand->ReFoldSingleStrand( percent, maxStructures, windowSize );
    error = checker->isErrorStatus( mainCalcError );

    // If no error occurred, print message that main calculation is done.
    if( error == 0 ) { cout << "done." << endl; }
  }

  /*
   * Write a CT output file using the WriteCt method.
   * After writing is complete, use the error checker's isErrorStatus method to
   * check for errors.
   */
  if( error == 0 ) {

    // Show a message saying that the CT file is being written.
    cout << "Writing output ct file..." << flush;

    // Write the CT file and check for errors.
    int writeError = strand->WriteCt( ctFile.c_str() );
    error = checker->isErrorStatus( writeError );

    // If no errors occurred, show a CT file writing completion message.
    if( error == 0 ) { cout << "done." << endl; }
  }

  // Delete the error checker and data structure.
  delete checker;
  delete strand;

  // Print confirmation of run finishing.
  if( error == 0 ) { cout << calcType << " complete." << endl; }
  else { cerr << calcType << " complete with errors." << endl; }

}

////////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

  refold* runner = new refold();
  bool parseable = runner->parse( argc, argv );
  if( parseable == true ) { runner->run(); }
  delete runner;
  return 0;

}
