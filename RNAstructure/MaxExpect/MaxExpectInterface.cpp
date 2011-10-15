/*
 * A program that predicts structures composed of probable base pairs and
 * single stranded nucleotides.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "MaxExpect.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
MaxExpect::MaxExpect() {

  // Initialize file names.
  input = 0;
  output = 0;

  // Initialize structure parameters.
  gamma = 1.0;
  percent = 50;
  structures = 1000;
  window = 5;

  // Initialize boolean flags.
  isRNA = true;
  isSequence = false;

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool MaxExpect::parse( int argc, char** argv ) {

  // Initialize arrays of flags without and with parameters.
  const char* flagsNoParams[] = {
    "-d", "-D", "--DNA",
    "--sequence"
  };

  const char* flagsParams[] = {
    "-g", "-G", "--gamma",
    "-p", "-P", "--percent",
    "-s", "-S", "--structures",
    "-w", "-W", "--window"
  };

  // Initialize usage string.
  const char* usageString =
    "USAGE: MaxExpect <input file> <ct file> [options]";

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

  string paramString = inputHelp + "\n\n" + ctHelp;

  // Initialize required parameters string.
  string dnaHelp = "\
-d, -D, --DNA\n\
     This flag only matters if the input file is a sequence file and has\n\
     been specified as such (see \"--sequence\" below).\n\
     Specifies that the sequence is DNA, and DNA parameters are to be used.\n\
     The default is to use RNA folding parameters.\
";

  string helpLabel = "\
-h, -H, --help\n\
     Displays this usage details message.\
";

  string seqLabel = "\
--sequence\n\
     Specifies that the input file is a sequence file.\
";

  string flagsNoParamsString =
    dnaHelp + "\n\n" + helpLabel + "\n\n" + seqLabel;

  // Initialize the flags with parameters string.
  string gammaHelp = "\
-g, -G, --gamma\n\
     Specifies the weight which should be put on base pairs.\n\
     Default is 1.0.\
";

  string percentHelp = "\
-p, -P, --percent\n\
     Specifies the percent difference allowed.\n\
     Default is 50 (ie, 50, not 0.5).\
";

  string maxStructuresHelp = "\
-s, -S, --structures\n\
     Specifies the maximum number of structures.\n\
     Default is 1000 structures.\
";

  string windowHelp = "\
-w, -W, --window\n\
     Specifies the window size.\n\
     Default is 5 nucleotides.\
";

  string flagsParamsString = gammaHelp + "\n\n" + percentHelp + "\n\n" +
    maxStructuresHelp + "\n\n" + windowHelp;

  // Initialize parser and check the command line
  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 2,
						   4, flagsNoParams,
						   12, flagsParams );
  parser->setUsageStrings( usageString, paramString.c_str(),
                           flagsNoParamsString.c_str(),
			   flagsParamsString.c_str() );
  bool noError = parser->checkLine();

  // If command line structure is OK proceed with argument retrieval, parsing
  if( noError ) {
    input = parser->getParameter( 1 );
    output = parser->getParameter( 2 );

    // Set non-parameterized flags. The RNA/DNA flag is set only if a sequence
    // flag is present.
    isSequence = parser->contains( "--sequence" );
    if( isSequence ) { isRNA = !parser->containsInGroup( "-d -D --DNA" ); }

    // Set parameterized options
    // If one of these options is not found on the command line, the action to
    // set that action is skipped.
    // If one throws an error, an error message is shown and parsing continues,
    // so all errors are shown if there are more than one.
    bool flagSet = true;

    // Gamma option
    flagSet = parser->setOption( gamma, "-g -G --gamma",
				 "Invalid gamma value given." );
    if( flagSet == false ) { noError = false; }

    // Percent difference option
    flagSet = parser->setOption( percent, "-p -P --percent",
                                 "Invalid percent difference given.", ">= 0" );
    if( flagSet == false ) { noError = false; }

    // Maximum structures option
    flagSet = parser->setOption( structures, "-s -S --structures",
                                 "Invalid maximum structures given.", ">= 1" );
    if( flagSet == false ) { noError = false; }

    // Window size option
    flagSet = parser->setOption( window, "-w -W --window",
                                 "Invalid window size given.", ">= 0" );
    if( flagSet == false ) { noError = false; }
  }

  // Delete the parser and return the error status.
  delete parser;
  return noError;

}

///////////////////////////////////////////////////////////////////////////////
// Run structure prediction.
///////////////////////////////////////////////////////////////////////////////
void MaxExpect::run() {

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
  RNA* strand = new RNA( input, type, isRNA );
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
   * Calculate accurate structures using the MaximizeExpectedAccuracy method.
   * After calculation is done, check error status of strand as above.
   */
  if( error == 0 ) {
    // Calculate structures.
    cout << "Calculating maximum expected accuracy structures..." << flush;

    error =
      strand->MaximizeExpectedAccuracy( percent, structures, window, gamma );

    // Check the error status
    if( !checker->isErrorStatus( error ) ) { cout << "done." << endl; }
  }

  /*
   * Write a CT file for the strand
   * Use method WriteCt to write the file.
   *
   * After writing is complete, check the error status of the strand and print
   * a final message announcing the finish of the run.
   */
  if( error == 0 ) {
    // Write the CT file
    cout << "Writing output ct file..." << flush;
    int writeError = strand->WriteCt( output );

    // Check the error status
    error = checker->isErrorStatus( writeError );
    if( error == 0 ) { cout << "done." << endl; }
  }

  // Delete the strand and checker, then print confirmation of run finishing.
  delete checker;
  delete strand;

  const char* message = "Calculation of maximum expected accuracy structures";
  if( error == 0 ) { cout << message << " complete." << endl; }
  else { cerr << message << " complete with errors." << endl; } 
}

// Main method to run the program
int main( int argc, char* argv[] ) {
  MaxExpect* myAccurate = new MaxExpect();
  bool parseable = myAccurate->parse( argc, argv );
  if( parseable == true ) { myAccurate->run(); }
  delete myAccurate;
  return 0;
}
