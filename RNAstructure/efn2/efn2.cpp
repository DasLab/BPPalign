/*
 * A program that calculates the free energy of a strand of nucleic acids.
 * These nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2008  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "efn2.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
efn2Interface::efn2Interface() {

  // Initialize the calculation type description.
  calcType = "efn2";

  // Initialize the nucleic acid type.
  isRNA = true;

  // Initialize flag that signifies streaming to standard output.
  stdPrint = false;

  // Initialize the calculation temperature.
  temperature = 310.15;

  // Initialize the flag that signifies writing a thermodynamic details file.
  writeTherm = false;

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool efn2Interface::parse( int argc, char** argv ) {

  // Initialize arrays of flags without and with parameters.
  const char* flagsNoParams[] = {
    "-d", "-D", "--DNA",
    "-p", "-P", "--print",
    "-w", "-w", "--writedetails"
  };

  const char* flagsParams[] = {
    "-t", "-T", "--temperature"
  };

  // Initialize usage string.
  string usageString =
    "USAGE: efn2 <ct file> <energy file> [options]";

  // Initialize required parameters string.
  string ctHelp = "\
<ct file>\n\
     The name of a file containing structure CT data.\
";

  string outHelp = "\
<output file>\n\
     The energy file to which output is written.\n\
     This output file can be written in one of two forms:\n\
     1. Simple list\n\
        Lists free energy for each structure, lowest first.\n\
     2. Thermodynamic details\n\
        Writes details of every substructure in each structure, and the\n\
        corresponding free energy of each.\
";

  string paramsString =
    ctHelp + "\n\n" + outHelp;

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

  string printHelp = "\
-p, -P, --print\n\
     Print the output file to standard output.\n\
     This won't override default behavior of writing to a file.\n\
     Thermodynamic files (if written) are not piped.\
";

  string detailHelp = "\
-w, -W, --writedetails\n\
     Write a thermodynamic details file.\n\
     The thermodynamic details file replaces a standard output (list) file.\
";

  string flagsNoParamsString =
    dnaHelp + "\n\n" + helpLabel + "\n\n" + printHelp + "\n\n" + detailHelp;

  // Initialize the flags with parameters string.
  string tempHelp = "\
-t, -T, --temperature\n\
     Specify the temperature at which calculation takes place in Kelvin.\n\
     Default is 310.15 K, which is 37 C.\
";

  string flagsParamsString =
    tempHelp;

  // Initialize parser and check the command line.
  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 2,
						   9, flagsNoParams,
						   3, flagsParams );
  parser->setUsageStrings( usageString.c_str(), paramsString.c_str(),
                           flagsNoParamsString.c_str(),
			   flagsParamsString.c_str() );
  bool noError = parser->checkLine();

  // If command line structure is OK, proceed with argument retrieval, parsing.
  if( noError ) {

    // Get required parameters.
    ctFile = parser->getParameter( 1 );
    outFile = parser->getParameter( 2 );

    // Set non-parameterized flag(s).
    isRNA = !parser->containsInGroup( "-d -D --DNA" );
    stdPrint = parser->containsInGroup( "-p -P --print" );
    writeTherm = parser->containsInGroup( "-w -W --writedetails" );

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
// Run efn2 calculations; print to standard output and/or write a thermodynamic
// details file as necessary.
///////////////////////////////////////////////////////////////////////////////
void efn2Interface::run() {

  // Create a variable to handle errors.
  int error = 0;

  /*
   * Use the constructor for RNA that specifies a filename.
   * Specify type = 1 (ct file).
   * isRNA identifies whether the strand is RNA (true) or DNA (false).
   *
   * After construction of the strand data structure, create the error checker
   * which monitors for errors.  
   * Throughout, the error status of the calculation is checked with a variant
   * of the isErrorStatus method, which returns 0 if no error occurred. The
   * calculation proceeds as long as error = 0.
   */
  cout << "Initializing nucleic acids..." << flush;
  RNA* strand = new RNA( ctFile.c_str(), 1, isRNA );
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
   * Intialize the structure number and energies vector to default values
   * for use later.
   * Regardless of what kind of file is written, use the error checker's
   * isErrorStatus method to check for errors.
   */
  int structures = 0;
  vector<double> energies;

  /*
   * If the user wants a thermodynamic details file, write that alone because
   * it is simply a more detailed version of the default output file.
   * Use the WriteThermodynamicDetails method to write this details file.
   */
  if( ( error == 0 ) && ( writeTherm ) ) {

    // Show a message saying that the details file is being written.
    cout << "Writing thermodynamic details file..." << flush;

    // Write the thermodynamic details file and check for errors.
    int thermError = strand->WriteThermodynamicDetails( outFile.c_str() );
    error = checker->isErrorStatus( thermError );
  }

  /*
   * Otherwise, if the user wants a simple output file, write it, and pipe to
   * standard output if the user requests it.
   * To get the free energy values for each structure, use the
   * CalculateFreeEnergy method.
   */
  else if( ( error == 0 ) && ( !writeTherm ) ) {

    // Show a message saying that the simple output file is being written.
    cout << "Writing free energy list file..." << flush;

    // Fetch the number of structures that are in the ct file.
    structures = strand->GetStructureNumber();

    // For each structure, calculate its energy. If no error occurred in this
    // calculation, put the energy value in the energies vector. If an error
    // did occur, set the error flag and stop calculations.
    for( int i = 1; i <= structures; i++ ) {
      double energy = strand->CalculateFreeEnergy( i );
      error = checker->isErrorStatus();

      if( error == 0 ) { energies.push_back( energy ); }
      else break;
    }

    // If all free energies were calculated correctly, write the output file.
    if( error == 0 ) {
      ofstream out( outFile.c_str() );
      for( int i = 1; i <= structures; i++ ) {
	out << "Structure: " << i << "   Energy = " << fixed
	    << setprecision( 1 ) << energies.at( i - 1 ) << endl;
      }
    }
  }

  /*
   * If the free energy calculations, whichever type, were finished correctly,
   * print a message that says so.
   */
  if( error == 0 ) { cout << "done." << endl; }

  /*
   * If the output is a simple list file and should be piped to standard output
   * then pipe it.
   */
  if( ( error == 0 ) && ( !writeTherm ) && ( stdPrint ) ) {

    // Write a header for the outputted file.
    cout << endl << "Generated output file: " << outFile << endl << endl;

    for( int i = 1; i <= structures; i++ ) {
      cout << "Structure: " << i << "   Energy = " << fixed
	   << setprecision( 1 ) << energies.at( i - 1 ) << endl;
    }

    // Add a blank line after the file is printed, for neatness.
    cout << endl;
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

  efn2Interface* runner = new efn2Interface();
  bool parseable = runner->parse( argc, argv );
  if( parseable == true ) { runner->run(); }
  delete runner;
  return 0;

}
