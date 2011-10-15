/*
 * A program that analyzes a strand of nucleic acids using stochastic
 * probability sampling.
 * These nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "stochastic.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
stochastic::stochastic() {

  // Initialize the type description.
  calcType = "Stochastic sampling";

  // Initialize the ensemble size.
  ensemble = 1000;

  // Initialize the random seed.
  seed = 1234;

  // Set the sequence flag to false.
  isSequence = false;

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool stochastic::parse( int argc, char** argv ) {

  // Initialize arrays of flags without and with parameters.
  const char* flagsNoParams[] = {
    "-d", "-D", "--DNA",
    "--sequence"
  };

  const char* flagsParams[] = {
    "-e", "-E", "--ensemble",
    "-s", "-S", "--seed"
  };

  // Initialize usage string.
  string usageString =
    "USAGE: stochastic <input file> <ct file> [options]";

  // Initialize required parameters string.
  string inputHelp = "\
<input file>\n\
     Specifies the input file. This file can be in two formats:\n\
     1. Partition function save file (holds probability data)\n\
     2. Sequence file (holds raw sequence: .seq or .fasta).\n\
        Note that in order to use a sequence file, the \"--sequence\" flag\n\
        must be specified (see \"--sequence\" below).\
";

  string ctHelp = "\
<ct file>\n\
     Specifies the ct file to which structures are written.\
";

  string paramString =
    inputHelp + "\n\n" + ctHelp;

  // Initialize required parameters string.
  string dnaHelp = "\
-d, -D, --DNA\n\
     This flag only matters if the input file is a sequence file and has\n\
     been specified as such (see \"--sequence\" below).\n\
     Specify that the sequence is DNA, and DNA parameters are to be used.\n\
     Default is to use RNA parameters.\
";

  string helpLabel = "\
-h, -H, --help\n\
     Display the usage details message.\
";

  string seqLabel = "\
--sequence\n\
     Specify that the input file is a sequence file.\
";

  string flagsNoParamsString =
    dnaHelp + "\n\n" + helpLabel + "\n\n" + seqLabel;

  // Initialize flags with parameters string.
  string ensembleHelp = "\
-e, -E, --ensemble\n\
     Specify the ensemble size.\n\
     Default is 1000 structures.\
";

  string seedHelp = "\
-s, -S, --seed\n\
     Specify random seed.\n\
     Default is 1234.\
";

  string flagsParamsString =
    ensembleHelp + "\n\n" + seedHelp;

  // Initialize parser and check the command line.
  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 2,
						   4, flagsNoParams,
						   6, flagsParams );
  parser->setUsageStrings( usageString.c_str(), paramString.c_str(),
                           flagsNoParamsString.c_str(),
			   flagsParamsString.c_str() );
  bool noError = parser->checkLine();

  // If command line structure is OK, proceed with argument retrieval, parsing.
  if( noError ) {

    // Get required parameters.
    inFile = parser->getParameter( 1 );
    ctFile = parser->getParameter( 2 );

    // Set non-parameterized flag(s).
    isSequence = parser->contains( "--sequence" );
    if( isSequence ) { isRNA = !parser->containsInGroup( "-d -D --DNA" ); }

    // Set parameterized options
    // If one of these options is not found on the command line, the action to
    // set that option is skipped.
    // If one throws an error, an error message is shown and parsing continues,
    // so all errors are shown if there are more than one.
    bool flagSet = true;

    // Ensemble size option
    flagSet = parser->setOption( ensemble, "-e -E --ensemble",
				 "Invalid ensemble size given.", "> 0" );
    if( flagSet == false ) { noError = false; }

    // Random seed option
    flagSet = parser->setOption( seed, "-s -S --seed",
                                 "Invalid random seed given.", "> 0" );
    if( flagSet == false ) { noError = false; }
  }

  // Delete the parser and return the error status.
  delete parser;
  return noError;

}

///////////////////////////////////////////////////////////////////////////////
// Run stochastic sampling.
///////////////////////////////////////////////////////////////////////////////
void stochastic::run() {

  // Create a variable to handle errors.
  int error = 0;

  /*
   * Use the constructor for RNA that specifies a filename.
   * If the input file is a pfs file, specify type = 3 (pfs file).
   * If the input file is a sequence file, specify type = 2 (sequence file).
   *
   * After construction of the strand, create the error checker which monitors
   * the strand for errors. Then, check for errors with the isErrorStatus
   * function, which returns 0 if no error occurs.
   * Throughout, the calculation proceeds as long as error = 0.
   */
  cout << "Initializing nucleic acids..." << flush;
  int type = ( !isSequence ) ? 3 : 2;
  RNA* strand = new RNA( inFile.c_str(), type, isRNA );
  ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
  error = checker->isErrorStatus();
  if( error == 0 ) { cout << "done." << endl; }

  /*
   * If the input file is a sequence file, calculate the partition function
   * first before calculating maximum expected accuracy.
   */
  if( error == 0 && isSequence ) {
    cout << "Calculating partition function..." << flush;

    // Run the partition function.
    error = strand->PartitionFunction();

    // Check the error status.
    if( !checker->isErrorStatus( error ) ) { cout << "done." << endl; }
  }

  /*
   * Do stochastic sampling using the SampleStochastic method.
   * After the main calculation is complete, use the error checker's
   * isErrorStatus method to check for errors.
   */
  if( error == 0 ) {

    // Show a message that the main calculation has started.
    cout << "Analyzing stochastic samples..." << flush;

    // Calculate stochastic sampling and check for errors.
    int mainCalcError =
        strand->Stochastic( ensemble, seed );
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

  stochastic* runner = new stochastic();
  bool parseable = runner->parse( argc, argv );
  if( parseable == true ) { runner->run(); }
  delete runner;
  return 0;

}
