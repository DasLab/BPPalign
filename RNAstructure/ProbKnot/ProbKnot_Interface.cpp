/*
 * A program that finds pseudoknots in a strand of nucleic acids.
 * These nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2010  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "ProbKnot_Interface.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
ProbKnot::ProbKnot() {

  // Initialize the calculation type description.
  calcType = "Prediction of Pseudoknots";

  // Initialize the nucleic acid type.
  isRNA = true;

  // Set the sequence flag to false.
  isSequence = false;

  // Initialize the number of iterations the calculation should undergo.
  iterations = 1;

  // Initialize the minimum helix length.
  minHelixLength = 3;

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool ProbKnot::parse( int argc, char** argv ) {

  // Initialize arrays of flags without and with parameters.
  const char* flagsNoParams[] = {
    "-d", "-D", "--DNA",
    "--sequence"
  };

  const char* flagsParams[] = {
    "-i", "-I", "--iterations",
    "-m", "-M", "--minimum"
  };

  // Initialize usage string.
  string usageString =
    "USAGE: ProbKnot <input file> <ct file> [options]";

  // Initialize required parameters string.
  string inputHelp = "\
<input file>\n\
     The name of the input file. This file can be in two formats:\n\
     1. Partition function save file (holds probability data)\n\
     2. Sequence file (holds raw sequence: .seq or .fasta).\n\
        Note that in order to use a sequence file, the \"--sequence\" flag\n\
        must be specified (see \"--sequence\" below).\
";

  string ctHelp = "\
<ct file>\n\
     The name of a CT file to which output will be written.\
";

  string paramsString =
     inputHelp + "\n\n" + ctHelp;

  // Initialize flags without parameters string.
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
  string iterationHelp = "\
-i, -I, --iterations\n\
     Specify the number of iterations the calculation will undergo.\n\
     Default is 1 iteration.\
";

  string minHelixHelp = "\
-m, -M, --minimum\n\
     Specify the minimum length accepted for a helix.\n\
     Default is 3 base pairs.\
";

  string flagsParamsString =
    iterationHelp + "\n\n" + minHelixHelp;

  // Initialize parser and check the command line.
  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 2,
                                                   4, flagsNoParams,
                                                   6, flagsParams );
  parser->setUsageStrings( usageString.c_str(), paramsString.c_str(),
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

    // Iterations option
    flagSet =
      parser->setOption( iterations, "-i -I --iterations",
			 "Invalid number of iterations given.", "> 0" );
    if( flagSet == false ) { noError = false; }

    // Minimum helix length option
    flagSet =
      parser->setOption( minHelixLength, "-m -M --minimum",
			 "Invalid minimum helix length given.", "> 0" );
    if( flagSet == false ) { noError = false; }
  }

  // Delete the parser and return the error state.
  delete parser;
  return noError;

}

///////////////////////////////////////////////////////////////////////////////
// Run ProbKnot pseudoknot prediction.
///////////////////////////////////////////////////////////////////////////////
void ProbKnot::run() {

  // Create a variable to handle errors.
  int error = 0;

  /*
   * Use the constructor for RNA that specifies a filename.
   * If the input file is a pfs file, specify type = 3 (pfs file).
   * If the input file is a sequence file, specify type = 2 (sequence file).
   * isRNA identifies whether the strand is RNA (true) or DNA (false).
   *
  * After construction of the strand data structure, create the error checker
   * which monitors for errors.
   * Throughout, the error status of the calculation is checked with a variant
   * of the isErrorStatus method, which returns 0 if no error occurred. The
   * calculation proceeds as long as error = 0.
   */
  cout << "Initializing nucleic acids..." << flush;
  int type = ( !isSequence ) ? 3 : 2;
  RNA* strand = new RNA( inFile.c_str(), type, isRNA );
  ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
  error = checker->isErrorStatus();
  if( error == 0 ) { cout << "done." << endl; }

  /*
   * If the input file is a sequence file, calculate the partition function
   * first before calculating pseudoknot prediction.
   * After the partition function is complete, use the error checker's
   * isErrorStatus method to check for errors.
   */
  if( error == 0 && isSequence ) {

    // Show a message saying that the partition function has been started.
    cout << "Calculating partition function..." << flush;

    // Run the partition function and check for errors.
    int partError = strand->PartitionFunction();
    error = checker->isErrorStatus( partError );

    // If no error occurred, print a message saying that partition function
    // has been finished.
    if( error == 0 ) { cout << "done." << endl; }
  }

  /*
   * Calculate pseudoknots using the ProbKnot method.
   * After the main calculation is complete, use the error checker's
   * isErrorStatus method to check for errors.
   */
  if( error == 0 ) {

    // Show a message that the main calculation has started.
    cout << "Calculating pseudoknots..." << flush;

    // Calculate maximum expected accuracy structures and check for errors.
    int mainCalcError = strand->ProbKnot( iterations, minHelixLength );
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

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

  ProbKnot* runner = new ProbKnot();
  bool parseable = runner->parse( argc, argv );
  if( parseable == true ) { runner->run(); }
  delete runner;
  return 0;

}
