/*
 * A program that folds two strands of nucleic acids.
 * These strands of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "bifold.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
bifold::bifold() {

  // Initialize the calculation type description.
  calcType = "Bimolecular folding";

  // Initialize the nucleic acid type.
  isRNA = true;

  // Initialize the allowance of intramolecular pairs.
  // False means that intramolecular pairs are allowed.
  intramolecular = false;

  // Initialize the maximum internal bulge loop size.
  maxLoop = 30;

  // Initialize the maximum number of structures.
  maxStructures = 20;

  // Initialize the maximum percent energy difference.
  percent = 50.0;

  // Initialize the calculation temperature.
  temperature = 310.15;

  // Initialize the folding window size.
  windowSize = 0;

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool bifold::parse( int argc, char** argv ) {

  // Initialize arrays of flags without and with parameters.
  const char* flagsNoParams[] = {
    "-d", "-D", "--DNA",
    "-i", "-I", "--intramolecular",
  };

  const char* flagsParams[] = {
    "-l", "-L", "--loop",
    "-m", "-M", "--maximum",
    "-p", "-P", "--percent",
    "-s", "-S", "--save",
    "-t", "-T", "--temperature",
    "-w", "-W", "--window"
  };

  // Initialize usage string.
  string usageString =
    "USAGE: bifold <seq file 1> <seq file 2> <ct file> [options]";

  // Initialize required parameters string.
  string seq1Help = "\
<seq file 1>\n\
     The name of a file containing a first input sequence.\
";

  string seq2Help = "\
<seq file 2>\n\
     The name of a file containing a second input sequence.\
";

  string ctHelp = "\
<ct file>\n\
     The name of a CT file to which output will be written.\
";

  string paramsString =
    seq1Help + "\n\n" + seq2Help + "\n\n" + ctHelp;

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

  string intramolecularHelp = "\
-i, -I, --intramolecular\n\
     Forbid intramolecular pairs (pairs within the same strand).\n\
     Default is to allow intramolecular pairs.\
";

  string flagsNoParamsString =
    dnaHelp + "\n\n" + helpLabel + "\n\n" + intramolecularHelp;

  // Initialize the flags with parameters string.
  string loopHelp = "\
-l, -L, --loop\n\
     Specify a maximum internal/bulge loop size.\n\
     Default is 30 unpaired nucleotides.\
";

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

  string saveHelp = "\
-s, -S, --save\n\
     Specify the name of a save file, needed for dotplots and refolding.\n\
     Default is not to generate a save file.\
";

  string tempHelp = "\
-t, -T, --temperature\n\
     Specify the temperature at which calculation takes place in Kelvin.\n\
     Default is 310.15 K, which is 37 degrees C.\
";

  string windowHelp = "\
-w, -W, --window\n\
     Specify a window size.\n\
     Default is 0 nucleotides.\
";

  string flagsParamsString =
    loopHelp + "\n\n" + maxStructuresHelp + "\n\n" + percentHelp + "\n\n" +
    saveHelp + "\n\n" + tempHelp + "\n\n" + windowHelp;

  // Initialize parser and check the command line.
  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 3,
						   6, flagsNoParams,
						   18, flagsParams );
  parser->setUsageStrings( usageString.c_str(), paramsString.c_str(),
                           flagsNoParamsString.c_str(),
			   flagsParamsString.c_str() );
  bool noError = parser->checkLine();

  // If command line structure is OK, proceed with argument retrieval, parsing.
  if( noError ) {

    // Get required parameters.
    seqFile1 = parser->getParameter( 1 );
    seqFile2 = parser->getParameter( 2 );
    ctFile = parser->getParameter( 3 );

    // Set non-parameterized flag(s).
    isRNA = !parser->containsInGroup( "-d -D --DNA" );
    intramolecular = !parser->containsInGroup( "-i -I --intramolecular" );

    // Set parameterized options
    // If one of these options is not found on the command line, the action to
    // set that option is skipped.
    // If one throws an error, an error message is shown and parsing continues,
    // so all errors are shown if there are more than one.
    bool flagSet = true;

    // Maximum loop size option
    flagSet = parser->setOption( maxLoop, "-l -L --loop",
				 "Invalid maximum loop size given." );
    if( flagSet == false ) { noError = false; }

    // Maximum number of structures option
    flagSet = parser->setOption( maxStructures, "-m -M --maximum",
				 "Invalid number of structures given." );
    if( flagSet == false ) { noError = false; }

    // Percent difference option
    flagSet = parser->setOption( percent, "-p -P --percent",
				 "Invalid percent difference given.", ">= 0" );
    if( flagSet == false ) { noError = false; }

    // Save file option
    const char* saveString = "";
    flagSet = parser->setOption( saveString, "-s -S --save",
				 "Error creating save file name.", false );
    if( flagSet == false ) { noError = false; }
    saveFile = saveString;

    // Temperature option
    flagSet = parser->setOption( temperature, "-t -T --temperature",
				 "Invalid temperature given.", ">= 0" );
    if( flagSet == false ) { noError = false; }

    // Window size option
    flagSet = parser->setOption( windowSize, "-w -W --window",
				 "Invalid window size.", ">= 0" );
    if( flagSet == false ) { noError = false; }
  }

  // Delete the parser and return the error state.
  delete parser;
  return noError;

}

///////////////////////////////////////////////////////////////////////////////
// Run bimolecular folding calculations.
///////////////////////////////////////////////////////////////////////////////
void bifold::run() {

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
  ErrorChecker<HybridRNA>* checker = new ErrorChecker<HybridRNA>( strand );
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
   * Set allowance or denial of intramolecular pairs using the
   * SetForbidIntramolecular method.
   * Since this method only sets a flag, error checking is not necessary.
   */
  if( error == 0 ) {
    strand->SetForbidIntramolecular( intramolecular );
  }

  /*
   * Fold the hybrid strand using the FoldBimolecular method.
   * If a save file name has been specified, the FoldBimolecular method also
   * has the ability to write a folding save file.
   * During calculation, monitor progress using the TProgressDialog class and
   * the Start/StopProgress methods of the RNA class. Neither of these methods
   * require any error checking.
   * After the main calculation is complete, use the error checker's
   * isErrorStatus method to check for errors in the main calculation.
   */
  if( error == 0 ) {

    // Show a message saying that the main calculation has started.
    cout << "Folding two strands..." << endl;

    // Create the progress monitor.
    TProgressDialog* progress = new TProgressDialog();
    strand->SetProgress( *progress );

    // Fold the hybrid strand and check for errors.
    int mainCalcError =
      strand->FoldBimolecular( percent, maxStructures, windowSize,
			       saveFile.c_str(), maxLoop );
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

  bifold* runner = new bifold();
  bool parseable = runner->parse( argc, argv );
  if( parseable == true ) { runner->run(); }
  delete runner;
  return 0;

}
