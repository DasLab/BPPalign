/*
 * A program that calculates the partition function for two strands of nucleic
 * acids.
 * These nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "bipartition.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
bipartition::bipartition() {

  // Initialize the calculation type description.
  calcType = "Bimolecular partition function";

  // Initialize the nucleic acid type.
  isRNA = true;

  // Initialize the calculation temperature.
  temperature = 310.15;

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool bipartition::parse( int argc, char** argv ) {

  // Initialize arrays of flags without and with parameters.
  const char* flagsNoParams[] = {
    "-d", "-D", "--DNA"
  };

  const char* flagsParams[] = {
    "-t", "-T", "--temperature"
  };

  // Initialize usage string.
  string usageString =
    "USAGE: bipartition <seq file 1> <seq file 2> <pfs file> [options]";

  // Initialize required parameters string.
  string seq1Help = "\
<seq file 1>\n\
     The name of a file containing a first input sequence.\
";

  string seq2Help = "\
<seq file 2>\n\
     The name of a file containing a second input sequence.\
";

  string pfsHelp = "\
<pfs file>\n\
     The name of a partition function save file to which output will be\n\
     written.\
";

  string paramsString =
    seq1Help + "\n\n" + seq2Help + "\n\n" + pfsHelp;

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
  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 3,
						   3, flagsNoParams,
						   3, flagsParams );
  parser->setUsageStrings( usageString.c_str(), paramsString.c_str(),
                           flagsNoParamsString.c_str(),
			   flagsParamsString.c_str() );
  bool noError = parser->checkLine();

  // If command line structure is OK, proceed with argument retrieval, parsing.
  if( noError ) {

    // Get required parameters.
    seqFile1 = parser->getParameter( 1 );
    seqFile2 = parser->getParameter( 2 );
    pfsFile = parser->getParameter( 3 );

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
// Run the bimolecular partition function calculation.
///////////////////////////////////////////////////////////////////////////////
void bipartition::run() {

  // Create a variable that handles errors.
  int error = 0;

  /*
   * Use the constructor for HybridRNA that specifies two filenames and types.
   * For both sequences, specify type = 2 (sequence file).
   * isRNA identifies whether the strand is RNA (true) or DNA (false).
   *
   * After construction of the strand data structure, create the error checker
   * which monitors for errors.
   * Throughout, the error status of the calculation is checked with a variant
   * of the isErrorStatus method, which returns 0 if no error occurred. The
   * calculation proceeds as long as error = 0.
   */
  cout << "Initializing nucleic acids..." << flush;
  HybridRNA* strand =
    new HybridRNA( seqFile1.c_str(), 2, seqFile2.c_str(), 2, isRNA );
  ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
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
   * Calculate the bimolecular partition function using the
   * PartitionFunctionBimolecular method.
   * During calculation, monitor progress using the TProgressDialog class and
   * the Start/StopProgress methods of the RNA class. Neither of these methods
   * require any error checking.
   * After the main calculation is complete, use the error checker's
   * isErrorStatus method to check for errors in the main calculation.
   */
  if( error == 0 ) {

    // Show a message saying that the main calculation has started.
    cout << "Calculating bimolecular partition function..." << endl;

    // Create the progress monitor.
    TProgressDialog* progress = new TProgressDialog();
    strand->SetProgress( *progress );

    // Calculate the bimolecular partition function and check for errors.
    int mainCalcError =
      strand->PartitionFunctionBimolecular( pfsFile.c_str() );
    error = checker->isErrorStatus( mainCalcError );

    // Delete the progress monitor.
    strand->StopProgress();
    delete progress;

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

  bipartition* runner = new bipartition();
  bool parseable = runner->parse( argc, argv );
  if( parseable == true ) { runner->run(); }
  delete runner;
  return 0;

}
