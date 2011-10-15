/*
 * A program that calculates the partition function for a strand of nucleic
 * acids.
 * These nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "partition.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
//////////////////////////////////////////////////////////////////////////////
partitionInterface::partitionInterface() {

  // Initialize the calculation type description.
  calcType = "Single strand partitioning";

  // Initialize the SHAPE intercept.
  intercept = -0.8;

  // Initialize the nucleic acid type.
  isRNA = true;

  // Initialize the maximum pairing distance between nucleotides.
  maxDistance = -1;

  // Initialize the SHAPE slope.
  slope = 2.6;

  // Initialize the calculation temperature.
  temperature = 310.15;

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool partitionInterface::parse( int argc, char** argv ) {

  // Initialize arrays of flags without and with parameters.
  const char* flagsNoParams[] = {
    "-d", "-D", "--DNA"
  };

  const char* flagsParams[] = {
    "-c", "-C", "--constraint",
    "-md", "-MD", "maxdistance",
    "-sh", "-SH", "--SHAPE",
    "-si", "-SI", "--SHAPEintercept",
    "-sm", "-SM", "--SHAPEslope",
    "-t", "-T", "--temperature"
  };

  // Initialize usage string.
  string usageString =
    "USAGE: partition <seq file> <pfs file> [options]";

  // Initialize required parameters string.
  string seqHelp = "\
<seq file>\n\
     The name of a file containing an input sequence.\
";

  string pfsHelp = "\
<pfs file>\n\
     The name of a partition function save file to which output will be\n\
     written.\
";

  string paramsString =
    seqHelp + "\n\n" + pfsHelp;

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
  string constraintHelp = "\
-c, -C, --constraint\n\
     Specify a constraints file to be applied.\n\
     Default is to have no constraints applied.\
";

  string maxDistanceHelp = "\
-md, -MD, --maxdistance\n\
     Specify a maximum pairing distance between nucleotides.\n\
     Default is no restriction on distance between pairs.\
";

   string SHAPEHelp = "\
-sh, -SH, --SHAPE\n\
     Specify a SHAPE constraints file to be applied.\n\
     Default is no SHAPE constraints applied.\n\
     This constraint utilizes the SHAPE pseudoenergy constraints.\
";

  string SHAPEinterceptHelp = "\
-si, -SI, --SHAPEintercept\n\
     Specify an intercept used with SHAPE constraints.\n\
     Default is -0.8 kcal/mol.\
";

  string SHAPEslopeHelp = "\
-sm, -SM, --SHAPEslope\n\
     Specify a slope used with SHAPE constraints.\n\
     Default is 2.6 kcal/mol.\
";

  string tempHelp = "\
-t, -T, --temperature\n\
     Specify the temperature at which calculation takes place in Kelvin.\n\
     Default is 310.15 K, which is 37 degrees C.\
";

  string flagsParamsString =
    constraintHelp + "\n\n" + maxDistanceHelp + "\n\n" + SHAPEHelp + "\n\n" +
    SHAPEinterceptHelp + "\n\n" + SHAPEslopeHelp + "\n\n" + tempHelp;

  // Initialize parser and check command line
  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 2,
						   3, flagsNoParams,
						   18, flagsParams );
  parser->setUsageStrings( usageString.c_str(), paramsString.c_str(),
                           flagsNoParamsString.c_str(),
			   flagsParamsString.c_str() );
  bool noError = parser->checkLine();


  // If command line structure is OK, proceed with argument retrieval, parsing.
  if( noError ) {

    // Get required parameters.
    seqFile = parser->getParameter( 1 );
    pfsFile = parser->getParameter( 2 );

    // Set non-parameterized flag(s).
    isRNA = !parser->containsInGroup( "-d -D --DNA" );

    // Set parameterized options
    // If one of these optionss is not found on the command line, the action to
    // set that option is skipped.
    // If one throws an error, an error message is shown and parsing continues,
    // so all error are shown if there are more than one.
    bool flagSet = true;

    // Constraint (folding) file option
    const char* constraintString = constraintFile.c_str();
    flagSet = parser->setOption( constraintString, "-c -C --constraint",
                                 "Folding constraint file not found.", true );
    constraintFile = constraintString;
    if( flagSet == false ) { noError = false; }

    // Maximum pairing distance option
    flagSet =
      parser->setOption( maxDistance, "-m -M --maxdistance",
                         "Invalid maximum pairing distance given.", "> 0" );
    if( flagSet == false ) { noError = false; }

    // SHAPE constraint file option
    const char* shapeString = SHAPEFile.c_str();
    flagSet = parser->setOption( shapeString, "-sh -SH --SHAPE",
                                 "SHAPE constraint file not found.", true );
    SHAPEFile = shapeString;
    if( flagSet == false ) { noError = false; }

    // SHAPE intercept option
    flagSet = parser->setOption( intercept, "-si -SI --SHAPEintercept",
                                 "Invalid SHAPE intercept given." );
    if( flagSet == false ) { noError = false; }

    // SHAPE slope option
    flagSet = parser->setOption( slope, "-sm -SM -SHAPEslope",
                                 "Invalid SHAPE slope given." );
    if( flagSet == false ) { noError = false; }

    // Temperature option
    flagSet = parser->setOption( temperature, "-t -T --temperature",
                                 "Invalid temperature given.", ">= 0" );
    if( flagSet == false ) { noError = false; }
  }

  // Delete the parser and return the error state
  delete parser;
  return noError;

}

///////////////////////////////////////////////////////////////////////////////
// Run the partition function calculation.
///////////////////////////////////////////////////////////////////////////////
void partitionInterface::run() {

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
   * Set maximum pairing distance using the ForceMaximumPairingDistance method.
   */
  if( error == 0 && maxDistance != -1 ) {

    // Show a message saying that the maximum pairing distance is being set.
    cout << "Setting maximum distance between paired nucleotides..." << flush;

    // Set the maximum pairing distance and check for errors.
    int distError = strand->ForceMaximumPairingDistance( maxDistance );
    error = checker->isErrorStatus( distError );

    // If no error occurred, print a message saying that maximum pairing
    // distance was set.
    if( error == 0 ) { cout << "done." << endl; }
  }

  /*
   * Add constraints if files have been given for their inclusion.
   * For folding constraints, use the ReadConstraints method.
   * For SHAPE constraints, use the ReadSHAPE method.
   * After each constraint type, check the error status of the strand as above.
   */

  bool applyConstraints =
    ( constraintFile != "" ) ||
    ( SHAPEFile != "" );

  if( error == 0 && applyConstraints ) {

    // Show a message saying that constraints are being applied.
    cout << "Applying constraints..." << flush;
    int constraintError = 0;

    // Read folding constraints, if applicable.
    if( constraintFile != "" ) {
      constraintError = strand->ReadConstraints( constraintFile.c_str() );
      error = checker->isErrorStatus( constraintError );
    }

    // Read SHAPE constraints, if applicable.
    if( error == 0 && SHAPEFile != "" ) {
      constraintError = strand->ReadSHAPE( SHAPEFile.c_str(), slope, intercept );
      error = checker->isErrorStatus( constraintError );
    }

    // If no error occurred, print a message saying that constraints were
    // applied.
    if( error == 0 ) { cout << "done." << endl; }
  }

  /*
   * Calculate the partition function using the  PartitionFunction method.
   * During calculation, monitor progress using the TProgressDialog class and
   * the Start/StopProgress methods of the RNA class. Neither of these methods
   * require any error checking.
   * After the main calculation is complete, use the error checker's
   * isErrorStatus method to check for errors in the main calculation.
   */
  if( error == 0 ) {

    // Show a message saying that the main calculation has started.
    cout << "Calculating partition function..." << endl;

    // Create the progress monitor.
    TProgressDialog* progress = new TProgressDialog();
    strand->SetProgress( *progress );

    // Calculate the partition function and check for errors.
    int mainCalcError = strand->PartitionFunction( pfsFile.c_str() );
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

  partitionInterface* runner = new partitionInterface();
  bool parseable = runner->parse( argc, argv );
  if( parseable == true ) { runner->run(); }
  delete runner;
  return 0;

}