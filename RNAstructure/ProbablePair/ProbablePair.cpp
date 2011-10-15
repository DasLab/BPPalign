/*
 * A program that generates structures based on levels of probable pairs.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "ProbablePair.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
ProbablePair::ProbablePair() {

  // Initialize the calculation type description.
  calcType = "Calculation of probable structures";

  // Initialize the probable pairing threshold.
  threshold = 0.0;

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool ProbablePair::parse( int argc, char** argv ) {

  // Initialize array of flags with parameters.
  const char* flagsParams[] = {
    "-t", "-T", "--temperature"
  };

  // Initialize usage string.
  string usageString =
    "USAGE: ProbablePair <pfs file> <ct file> [options]";

  // Initialize required parameters string.
  string pfsHelp = "\
<pfs file>\n\
     The name of the input partition function save file.\
";

  string ctHelp = "\
<ct file>\n\
     The name of a CT file to which output will be written.\
";

  string paramsString =
    pfsHelp + "\n\n" + ctHelp;

  // Initialize flags without parameters string.
  string helpLabel = "\
-h, -H, --help\n\
     Display the usage details message.\
";

  string flagsNoParamsString =
    helpLabel;

  // Initialize flags with parameters string.
  string thresholdHelp = "\
-t, -T, --threshold\n\
     Specify pairing threshold for pairs to be included.\n\
     This should be expressed as a number >= 0.5 and <= 1.0.\n\
     Default is 0, which signifies structures should be generated at\n\
     multiple thresholds:\n\
          >= 0.99, >= 0.97, >= 0.95, >= 0.90, >= 0.80, >= 0.70, >= 0.60,\n\
          and >= 0.50.\
";

  string flagsParamsString =
    thresholdHelp;

  // Initialize parser and check the command line.
  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 2,
						   0, NULL,
						   3, flagsParams );
  parser->setUsageStrings( usageString.c_str(), paramsString.c_str(),
                           flagsNoParamsString.c_str(),
                           flagsParamsString.c_str() );
  bool noError = parser->checkLine();

  // If command line structure is OK, proceed with argument retrieval, parsing.
  if( noError ) {

    // Get required parameters.
    input = parser->getParameter( 1 );
    output = parser->getParameter( 2 );

    // Set parameterized options
    // If one of these options is not found on the command line, the action to
    // set that option is skipped.
    // If one throws an error, an error message is shown and parsing continues,
    // so all errors are shown if there are more than one.
    bool flagSet = true;

    // Threshold option
    flagSet = parser->setOption( threshold, "-t -T --threshold",
				 "Invalid threshold given.", ">= 0" );
    if( flagSet == false ) { noError = false; }
  }

  // Delete the parser and return the error status.
  delete parser;
  return noError;

}

///////////////////////////////////////////////////////////////////////////////
// Run probable structure generation.
///////////////////////////////////////////////////////////////////////////////
void ProbablePair::run() {

  // Create a variable to handle errors.
  int error = 0;

  /*
   * Use the constructor for RNA that specifies a filename.
   * Specify type = 3 (pfs file).
   *
   * After construction of the strand data structure, create the error checker
   * which monitors for errors.  
   * Throughout, the error status of the calculation is checked with a variant
   * of the isErrorStatus method, which returns 0 if no error occurred. The
   * calculation proceeds as long as error = 0.
   */
  cout << "Initializing nucleic acids..." << flush;
  RNA* strand = new RNA( input.c_str(), 3 );
  ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
  error = checker->isErrorStatus();
  if( error == 0 ) { cout << "done." << endl; }

  /*
   * Predict probable pairs using the PredictProbablePairs method.
   * After the main calculation is complete, use the error checker's
   * isErrorStatus method to check for errors.
   */
  if( error == 0 ) {

    // Show a message that the main calculation has started.
    cout << "Determining probable structures..." << flush;

    // Calculate maximum expected accuracy structures and check for errors.
    int mainCalcError = strand->PredictProbablePairs( threshold );
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
    int writeError = strand->WriteCt( output.c_str() );
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

  ProbablePair* runner = new ProbablePair();
  bool parseable = runner->parse( argc, argv );
  if( parseable == true ) { runner->run(); }
  delete runner;
  return 0;

}
