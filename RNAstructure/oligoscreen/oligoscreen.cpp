/*
 * A program that calculates thermodynamic parameters from a set of
 * oligonucleotides.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "oligoscreen.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
oligoscreen::oligoscreen() {

  // Initialize the calculation type description.
  calcType = "OligoScreen";

  // Initialize the nucleic acid type.
  isRNA = true;

  // Initialize the calculation temperature.
  temperature = 310.15;

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool oligoscreen::parse( int argc, char** argv ) {

  // Initialize arrays of flags without and with parameters.
  const char* flagsNoParams[] = {
    "-d", "-D", "--DNA"
  };

  const char* flagsParams[] = {
    "-t", "-T", "--temperature"
  };

  // Initialize usage string.
  string usageString =
    "USAGE: oligoscreen <list file> <report file> [options]";

  // Initialize required parameters string.
  string listHelp = "\
<list file>\n\
     The name of a file that holds the oligo list.\
";

  string reportHelp = "\
<report file>\n\
     Specify a report file to which output will be written.\
";

  string paramsString =
    listHelp + "\n\n" + reportHelp;

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

  // Initialize the flags with parameters string.
  string tempHelp = "\
-t, -T, --temperature\n\
     Specify the temperature at which calculation takes place in Kelvin.\n\
     Default is 310.15 K, which is 37 degrees C.\
";

  string flagsParamsString =
    tempHelp;

  // Initialize parser and check the command line.
  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 2,
						   3, flagsNoParams,
						   3, flagsParams );
  parser->setUsageStrings( usageString.c_str(), paramsString.c_str(),
                           flagsNoParamsString.c_str(),
			   flagsParamsString.c_str() );
  bool noError = parser->checkLine();

  // If command line structure is OK, proceed with argument retrieval, parsing.
  if( noError ) {

    // Get required parameters.
    listFile = parser->getParameter( 1 );
    reportFile = parser->getParameter( 2 );

    // Set non-parameterized flag(s).
    isRNA = !parser->containsInGroup( "-d -D --DNA" );

    // Set parameterized options
    // If one of these options is not found on the command line, the action to
    // set that option is skipped.
    // If one throws an error, an error message is shown and parsing continues,
    // so all errors are shown if there are more than one.
    bool flagSet = true;

    // Temperature option
    flagSet = parser->setOption( temperature, "-t -T --temperature",
				 "Invalid temperature given.", ">= 0" );
    if( flagSet == false ) { noError = false; }
  }

  // Delete the parser and return the error state.
  delete parser;
  return noError;

}

///////////////////////////////////////////////////////////////////////////////
// Run OligoScreen.
///////////////////////////////////////////////////////////////////////////////
void oligoscreen::run() {

  // Create a variable that handles errors.
  int error = 0;

  /*
   * Use the constructor for Oligowalk_object that specifies a nucleic acid
   * type.
   * isRNA identifies whether the strand is RNA (true) or DNA (false).
   *
   * After construction of the strand data structure, create the error checker
   * which monitors for errors.  
   * Throughout, the error status of the calculation is checked with a variant
   * of the isErrorStatus method, which returns 0 if no error occurred. The
   * calculation proceeds as long as error = 0.
   */
  cout << "Initializing nucleic acids..." << flush;
  Oligowalk_object* strand = new Oligowalk_object( isRNA );
  ErrorChecker<Oligowalk_object>* checker =
    new ErrorChecker<Oligowalk_object>( strand );
  error = checker->isErrorStatus();
  if( error == 0 ) { cout << "done." << endl; }

  /*
   * Set the temperature using the SetTemperature method.
   * Only set the temperature if a given temperature doesn't equal the default.
   * If the temperature does need to be set, use the error checker's
   * isErrorStatus method to check for errors.
   */
  if( ( error == 0 ) && ( temperature != 310.15 ) ) {

    // Show a message saying that the temperature is being set.
    cout << "Setting temperature..." << flush;

    // Set the temperature and check for errors.
    int tempError = strand->SetTemperature( temperature );
    error = checker->isErrorStatus( tempError );

    // If no error occurred, print a message saying that temperature is set.
    if( error == 0 ) { cout << "done." << endl; }
  }

  /*
   * Run OligoScreen using the OligoScreen method.
   * After the main calculation is complete, use the error checker's
   * isErrorStatus method to check for errors.
   */
  if( error == 0 ) {

    // Show a message that the main calculation has started.
    cout << "Running OligoScreen..." << flush;

    // Calculate maximum expected accuracy structures and check for errors.
    int mainCalcError =
      strand->OligoScreen( listFile.c_str(), reportFile.c_str() );
    error = checker->isErrorStatus( mainCalcError );

    // If no error occurred, print message that main calculation is done.
    if( error == 0 ) { cout << "done." << endl; }
  }

  // Delete the error checker and data structure.
  delete checker;
  delete strand;

  // Print confirmation of run finishing.
  if( error == 0 ) { cout << calcType << " complete." << endl; }
  else { cerr << calcType << " complete with errors." << endl; }

}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

  oligoscreen* runner = new oligoscreen();
  bool parseable = runner->parse( argc, argv );
  if( parseable == true ) { runner->run(); }
  delete runner;
  return 0;

}
