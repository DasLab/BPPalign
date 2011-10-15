/*
 * A program that runs CircleCompare to compare two different CT file
 * structures.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "CircleCompare_Interface.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
CircleCompare_Interface::CircleCompare_Interface() {

  // Initialize the calculation type description.
  calcType = "Circular structure comparison";

  // Initialize the optional annotation file names.
  probabilityFile = "";
  SHAPEFile = "";

  // Set usage of the alternative color scheme to false.
  alternative = false;

  // Set enforcing of exact prediction to false.
  exact = false;

  // Set comparison of the first structure only to false, so all structures in
  // the predicted CT file can be compared if necessary.
  firstOnly = false;

  // Set the number of structures in the predicted strand to a default of 1. No
  // matter what, if the program is run correctly, there will be at least one
  // predicted structure compared.
  predictedCount = 1;

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool CircleCompare_Interface::parse( int argc, char** argv ) {

  // Initialize arrays of flags without and with parameters.
  const char* flagsNoParams[] = {
    "-a", "-A", "--alternative",
    "-e", "-E", "--exact",
    "-f", "-F", "--firstOnly"
  };

  const char* flagsParams[] = {
    "-p", "-P", "--probability",
    "-s", "-S", "--SHAPE"
  };

  // Initialize usage string.
  string usageString =
    "USAGE: CircleCompare <predicted ct> <accepted ct> <ps file> [options]";

  // Initialize required parameters string.
  string predictedHelp = "\
<predicted ct>\n\
     The name of a file containing CT data for the predicted structure.\
";

  string acceptedHelp = "\
<accepted ct>\n\
     The name of a file containing CT data for the accepted structure.\
";

  string outputHelp = "\
<ps file>\n\
     The name of a Postscript image file to which output will be written.\
";

  string paramsString =
    predictedHelp + "\n\n" + acceptedHelp + "\n\n" + outputHelp;

  // Initialize flags without parameters string.
  string alternativeHelp = "\
-a, -A, --alternative\n\
     Specify that an alternative color scheme should be used.\n\
     Default is not to use the alternative color scheme.\
";

  string exactHelp = "\
-e, -E, --exact\n\
     Specify exact comparison when structure comparison is scored.\n\
     Default is to allow flexible pairings.\
";

  string firstHelp = "\
-f, -F, --firstOnly\n\
     Specify that only the first structure in the predicted ct file should\n\
     be compared with the accepted structure.\n\
     Default is that all structures in the predicted ct file are compared\n\
     with the accepted structure.\
";

  string helpLabel = "\
-h, -H, --help\n\
     Display the usage details message.\
";

  string flagsNoParamsString =
    alternativeHelp + "\n\n" + exactHelp + "\n\n" + firstHelp + "\n\n" +
    helpLabel;

  // Initialize the flags with parameters string.
  string probabilityHelp = "\
-p, -P, --probability\n\
     Specify the name of the file from which base pairing probability data\n\
     will be read for annotation.\n\
     This file should describe pairing data for the predicted structure, not\n\
     the accepted structure.\n\
     Default is no probability annotation file used.\
";

  string shapeHelp = "\
-s, -S, --SHAPE\n\
     Specify the name of the file from which SHAPE data will be read for\n\
     annotation.\n\
     Default is no SHAPE annotation file used.\
";

  string flagsParamsString =
    probabilityHelp + "\n\n" + shapeHelp;

  // Initialize parser and check the command line.
  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 3,
                                                   9, flagsNoParams,
                                                   6, flagsParams );
  parser->setUsageStrings( usageString.c_str(), paramsString.c_str(),
                           flagsNoParamsString.c_str(),
			   flagsParamsString.c_str() );
  bool noError = parser->checkLine();

  // If command line structure is OK, proceed with argument retrieval, parsing.
  if( noError ) {

    // Get required parameters.
    predicted = parser->getParameter( 1 );
    accepted = parser->getParameter( 2 );
    output = parser->getParameter( 3 );

    // Set non-parameterized flag(s).
    alternative = parser->containsInGroup( "-a -A --alternative" );
    exact = parser->containsInGroup( "-e -E --exact" );
    firstOnly = parser->containsInGroup( "-f -F --firstOnly" );

    // Set parameterized options.
    // If one of these options is not found on the command line, the action to
    // set that option is skipped.
    // If one throws an error, an error message is shown and parsing continues,
    // so all errors are shown if there are more than one.
    bool flagSet = true;

    // Probability annotation file option
    const char* probabilityString = "";
    flagSet = parser->setOption( probabilityString, "-p -P --probability",
				 "Probability annotation file not found.",
				 true );
    if( flagSet == false ) { noError = false; }
    probabilityFile = probabilityString;

    // SHAPE annotation file option
    const char* SHAPEstring = "";
    flagSet = parser->setOption( SHAPEstring, "-s -S --SHAPE",
				 "SHAPE annotation file not found.", true );
    if( flagSet == false ) { noError = false; }
    SHAPEFile = SHAPEstring;
  }

  // Delete the parser and return the error state.
  delete parser;
  return noError;

}

///////////////////////////////////////////////////////////////////////////////
// Run circle comparisons.
///////////////////////////////////////////////////////////////////////////////
void CircleCompare_Interface::run() {

  // Print out message saying comparison has started.
  cout << "Assembling circular structure comparison file..." << endl;

  // If the output file already exists, delete the file so a new file with that
  // name can be generated.
  ifstream test( output.c_str() );
  bool exists = test.good();
  test.close();
  if( exists ) { remove( output.c_str() ); }

  // Create a CircleCompare object.
  CircleCompare* circle = new CircleCompare( predicted, accepted, output );

  // For each predicted structure, compare it to the accepted structure.
  bool noError = true;
  for( int i = 1; i <= predictedCount; i++ ) {
    noError = circle->writeCirclePlot( i, alternative, exact,
				       probabilityFile, SHAPEFile );
    if( noError ) {
      cout << "Comparison " << i << " of " << predictedCount << " finished."
	   << endl;
    } else break;
  }

  // Delete the CircleCompare object.
  delete circle;

  // Print confirmation of run finishing.
  if( noError == true ) { cout << calcType << " complete." << endl; }
  else { cerr << calcType << " complete with errors." << endl; }

}

///////////////////////////////////////////////////////////////////////////////
// Validate the input CT files.
///////////////////////////////////////////////////////////////////////////////
bool CircleCompare_Interface::validate() {

  // Show message that CT file validation is in progress.
  cout << "Checking input CT files..." << flush;

  // Initialize array for sequence length 1 and 2.
  int lengths[2];

  // Create a temporary strand for each CT file and get its length.
  for( int i = 1; i <= 2; i++ ) {

    // Create the temporary strand and check it for errors.
    string file = ( i == 1 ) ? predicted : accepted;
    RNA* tempStrand = new RNA( file.c_str(), Postscript_Constants::CT_TYPE );
    bool error = false;
    int code = tempStrand->GetErrorCode();

    // If an error occurred, show an error message and set the error flag.
    if( code != 0 ) {
      cerr << endl << tempStrand->GetErrorMessage( code ) << endl;
      error = true;
    }

    // If an error didn't happen, do specific checking for each strand.
    else {

      // Set the sequence length of the strand.
      lengths[i-1] = tempStrand->GetSequenceLength();

      // If checking the predicted strand, get the count of how many possible
      // structures are in it. However, only do this if all possible predicted
      // structures should be compared separately.
      if( i == 1 && !firstOnly ) {
	predictedCount = tempStrand->GetStructureNumber();
      }

      // If checking the accepted strand, make sure there's only one structure
      // in that strand. If there is not exactly one, show error, return false.
      if( ( i == 2 ) && ( tempStrand->GetStructureNumber() != 1 ) ) {
	cerr << endl
	     << "Accepted CT file must contain exactly one structure." << endl;
	return false;
      }
    }

    // Delete the temporary strand.
    // If an error occurred during the check, return false.
    delete tempStrand;
    if( error ) { return false; }
  }

  // If the two lengths calculated from the temporary strands are not the same,
  // show an error message and return false.
  if( lengths[0] != lengths[1] ) {

    cerr << endl << "Sequences in CT files are not the same length."
	 << endl;
    return false;
  }

  // Print out a message signaling that the input has passed checking stage,
  // and return true.
  cout << "done." << endl;
  return true;

}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the interface.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

  CircleCompare_Interface* compare = new CircleCompare_Interface();

  bool valid = compare->parse( argc, argv );
  if( valid ) { valid = compare->validate(); }
  if( valid ) { compare->run(); }

  delete compare;

}
