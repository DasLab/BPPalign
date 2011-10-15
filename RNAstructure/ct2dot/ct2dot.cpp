/*
 * A program that converts a CT file to a dot bracket file.
 *
 * (c) 2010  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "ct2dot.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
ct2dot_Interface::ct2dot_Interface() {

  // Initialize the calculation type description.
  calcType = "CT file conversion";

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool ct2dot_Interface::parse( int argc, char** argv ) {

  // Initialize usage string.
  string usageString =
    "USAGE: ct2dot <ct file> <structure number> <bracket file>";

  // Initialize required parameters string.
  string ctHelp = "\
<ct file>\n\
     The name of a file containing the CT structure to convert.\
";

  string numberHelp = "\
<structure number>\n\
     The number, one-indexed, of the structure to convert in the CT file.\
";

  string bracketHelp = "\
<bracket file>\n\
     The name of a dot bracket file to which output will be written.\
";

  string paramsString =
    ctHelp + "\n\n" + numberHelp + "\n\n" + bracketHelp;

  // Initialize flags without parameters string.
  string helpLabel = "\
-h, -H, --help\n\
     Display the usage details message.\
";

  string flagsNoParamsString =
    helpLabel;

  // Initialize the flags with parameters string.
  string flagsParamsString =
    "";

  // Initialize parser and check the command line.
  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 3,
						   0, NULL,
						   0, NULL );
  parser->setUsageStrings( usageString.c_str(), paramsString.c_str(),
                           flagsNoParamsString.c_str(),
			   flagsParamsString.c_str() );
  bool noError = parser->checkLine();

  // If command line structure is OK, proceed with argument retrieval, parsing.
  if( noError ) {

    // Get required parameters.
    // Get the file names first.
    ctFile = parser->getParameter( 1 );
    bracketFile = parser->getParameter( 3 );

    // Get the structure number, which is also a required parameter.
    int testInt = 0;
    stringstream testStream( parser->getParameter( 2 ) );
    if( testStream >> testInt ) {
      number = testInt;
    } else {
      cerr << "Invalid structure number given." << endl;
      noError = false;
    }
  }

  // Delete the parser and return the error state.
  delete parser;
  return noError;

}

///////////////////////////////////////////////////////////////////////////////
// Run the CT file conversion calculation.
///////////////////////////////////////////////////////////////////////////////
void ct2dot_Interface::run() {

  // Show a message saying that conversion has started.
  cout << "Converting CT file..." << flush;

  // Create a variable that handles errors.
  int error = 0;

  // Initialize and open the CT file to convert.
  structure ct;
  openct( &ct, ctFile.c_str() );

  // Initialize and open the output dot bracket file.
  ofstream out;
  out.open( bracketFile.c_str() );

  // Write the comment line.
  // Make sure that the line ends in a newline, if not, add a newline.
  out << ">" << ct.ctlabel[number];
  if( ct.ctlabel[number][strlen(ct.ctlabel[number])-1] != '\n' ) {
    out << "\n";
  }

  // Write the sequence.
  for( int i = 1; i <= ct.numofbases; i++ ) {
    out << ct.nucs[i];
  }
  out << "\n";

  // Write the conversion mask.
  for( int i = 1; i <= ct.numofbases; i++ ) {
    if( ct.basepr[number][i] > i ) { out << "("; }
    else if( ct.basepr[number][i] == 0 ) { out << "."; }
    else { out << ")"; }
  }
  out << "\n";

  // Close the dot bracket file.
  out.close();

  // Show a message saying conversion is done.
  cout << "done." << endl;

  // Print confirmation of run finishing.
  if( error == 0 ) { cout << calcType << " complete." << endl; }
  else { cerr << calcType << " complete with errors." << endl; }

}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

  ct2dot_Interface* runner = new ct2dot_Interface();
  bool parseable = runner->parse( argc, argv );
  if( parseable == true ) { runner->run(); }
  delete runner;
  return 0;

}
