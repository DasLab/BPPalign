/*
 * A program that finds all suboptimal structures within a predefined small
 * increment for a strand of nucleic acids.
 * These nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "AllSub.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
AllSub::AllSub() {

  // Initialize the calculation type description.
  calcType = "Generation of suboptimal structures";

  // Initialize the absolute energy difference.
  absolute = -1.0;

  // Initialize the nucleic acid type.
  isRNA = true;

  // Initialize the maximum percent energy difference.
  percent = -1.0;

  // Initialize the calculation temperature.
  temperature = 310.15;

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool AllSub::parse( int argc, char** argv ) {

  // Initialize arrays of flags without and with parameters.
  const char* flagsNoParams[] = {
    "-d", "-D", "--DNA"
  };

  const char* flagsParams[] = {
    "-a", "-A", "--absolute",
    "-c", "-C", "--constraint",
    "-p", "-P", "--percent",
    "-t", "-T", "--temperature"
  };

  // Initialize usage string.
  string usageString =
    "USAGE: AllSub <seq file> <ct file> [options]";

  // Initialize required parameters string.
  string seqHelp = "\
<seq file>\n\
     The name of a file containing an input sequence.\
";

  string ctHelp = "\
<ct file>\n\
     The name of a CT file to which output will be written.\
";

  string paramsString =
    seqHelp + "\n\n" + ctHelp;

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
  string absoluteHelp = "\
-a, -A, --absolute\n\
     Specify a maximum absolute energy difference.\n\
     Default is determined by the length of the sequence.\
";

  string constraintHelp = "\
-c, -C, --constraint\n\
     Specify a constraints file to be applied.\n\
     Default is to have no constraints applied.\
";

  string percentHelp = "\
-p, -P, --percent\n\
     Specify a maximum percent energy difference.\n\
     Default is determined by the length of the sequence.\
";

  string tempHelp = "\
-t, -T, --temperature\n\
     Specify the temperature at which calculation takes place in Kelvin.\n\
     Default is 310.15 K, which is 37 degrees C.\
";

  string flagsParamsString = 
    absoluteHelp + "\n\n" + constraintHelp + "\n\n" + percentHelp + "\n\n" +
    tempHelp; 

  // Initialize parser and check the command line.
  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 2,
						   3, flagsNoParams,
						   12, flagsParams );
  parser->setUsageStrings( usageString.c_str(), paramsString.c_str(),
			   flagsNoParamsString.c_str(),
			   flagsParamsString.c_str() );
  bool noError = parser->checkLine();

  // If command line structure is OK, proceed with argument retrieval, parsing.
  if( noError ) {

    // Get required parameters.
    seqFile = parser->getParameter( 1 );
    ctFile = parser->getParameter( 2 );

    // Set non-parameterized flag(s).
    isRNA = !parser->containsInGroup( "-d -D --DNA" );

    // Set parameterized options
    // If one of these options is not found on the command line, the action to
    // set that option is skipped.
    // If one throws an error, an error message is shown and parsing continues,
    // so all errors are shown if there are more than one.
    bool flagSet = true;

    // Absolute energy difference option
    flagSet = 
      parser->setOption( absolute, "-a -A --absolute",
			 "Invalid absolute energy difference given.", ">= 0" );
    if( flagSet == false ) { noError = false; }

    // Constraint file option
    const char* constraintString = "";
    flagSet = parser->setOption( constraintString, "-c -C --constraint",
				 "Constraint file not found.", true );
    constraintFile = constraintString;
    if( flagSet == false ) { noError = false; }

    // Percent difference option
    flagSet = parser->setOption( percent, "-p -P --percent",
				 "Invalid percent difference given.", ">= 0" );
    if( flagSet == false ) { noError = false; }

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
// Run the suboptimal structure calculations.
///////////////////////////////////////////////////////////////////////////////
void AllSub::run() {

  // Create a variable that handles errors.
  int error = 0;

  /*
   * Use the constructor for RNA that specifies a filename.
   * Specify type = 2 (sequence file).
   * isRNA identifies whether the strand is RNA (true) or DNA (false).
   *
   * After construction of the strand data structure, create the error checker
   * which monitors for errors.  
   * Throughout, the error status of the calculation is checked with a variant
   * of the isErrorStatus method, which returns 0 if no error occurred. The
   * calculation proceeds as long as error = 0.
   */
  cout << "Initializing nucleic acids..." << flush;
  RNA* strand = new RNA( seqFile.c_str(), 2, isRNA );
  ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
  error = checker->isErrorStatus();
  if( error == 0 ) { cout << "done." << endl; }

  /*
   * Set the percent difference and the maximum absolute energy difference.
   * Both differences are based on the length of the sequence, and their values
   * in relation to the sequence length are hardcoded.
   * Only set these values if they were not set on the command line. As one or
   * the other of them may or may not be set, each one is checked individually.
   *
   * Get the length of the sequence using the GetSequenceLength method.
   * Since this method only returns a length, error checking is not necessary.
   */
  if( error == 0 ) {

    // Get the sequence length to identify the hardcoded thresholds.
    int length = strand->GetSequenceLength();

    // Set the maximum percent difference, if applicable.
    if( percent == -1 ) {
      percent =
	( length > 1200 ) ? 5 :
	( length > 800 ) ? 8 :
	( length > 500 ) ? 10 :
	( length > 300 ) ? 15 :
	( length > 120 ) ? 20 :
	( length > 50 ) ? 25 :
	50;
    }

    // Set the absolute energy difference, if applicable.
    if( absolute == -1 ) {
      absolute =
	( length > 1200 ) ? 0.25 :
	( length > 800 ) ? 0.5 :
	( length > 500 ) ? 0.75 :
	( length > 300 ) ? 1 :
	( length > 120 ) ? 1.5 :
	( length > 50 ) ? 3 :
	10;
    }
  }

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
   * Add constraints if a file has been given for their inclusion.
   * Read in this constraints with method ReadConstraints.
   * After constraints are read, use the error checker's isErrorStatus method
   * to check for errors.
   */
  if( error == 0 && constraintFile != "" ) {

    // Show a message saying that constraints are being applied.
    cout << "Applying constraints..." << flush;

    // Apply constraints and check for errors.
    int constraintError = strand->ReadConstraints( constraintFile.c_str() );
    error = checker->isErrorStatus( constraintError );

    // If no error occurred, print a message saying constraints were included.
    if( error == 0 ) { cout << "done." << endl; }
  }

  /*
   * Generate structures using the GenerateAllSuboptimalStructures method.
   * During calculation, monitor progress using the TProgressDialog class and
   * the Start/StopProgress methods of the RNA class. Neither of these methods
   * require any error checking.
   * After the main calculation is complete, use the error checker's
   * isErrorStatus method to check for errors.
   */
  if( error == 0 ) {

    // Show a message saying that the main calculation has started.
    cout << "Generating suboptimal structures..." << endl;

    // Create the progress monitor.
    TProgressDialog* progress = new TProgressDialog();
    strand->SetProgress( *progress );

    // Generate the suboptimal structures and check for errors.
    int mainCalcError =
      strand->GenerateAllSuboptimalStructures( (float)percent, absolute );
    error = checker->isErrorStatus( mainCalcError );

    // Delete the progress monitor.
    strand->StopProgress();
    delete progress;

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

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

  AllSub* runner = new AllSub();
  bool parseable = runner->parse( argc, argv );
  if( parseable == true ) { runner->run(); }
  delete runner;
  return 0;

}
