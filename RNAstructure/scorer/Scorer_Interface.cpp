/*
 * An interface file for a program that scores two structures and outputs their
 * sensitivity and PPV.
 * These structures can be composed of either DNA or RNA.
 *
 * (c) 2010  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "Scorer_Interface.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Scorer_Interface::Scorer_Interface() {

  // Set the boolean flags to their default values.
  exact = false;
  print = false;

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool Scorer_Interface::parse( int argc, char** argv ) {

  // Initialize array of flags without parameters.
  const char* flagsNoParams[] = {
    "-e", "-E", "--exact",
    "-p", "-P", "--print"
  };

  // Initialize usage string.
  string usageString =
    "USAGE: scorer <predicted ct> <accepted ct> <output file>";

  // Initialize required parameters string.
  string predictedHelp = "\
<predicted ct>\n\
     Specifies the predicted ct file: structure being tested.\
";

  string acceptedHelp = "\
<accepted ct>\n\
     Specifies the accepted ct file: structure known to be correct.\
";

  string outputHelp = "\
<output file>\n\
     Specifies the output file to which scores are written.\
";

  string paramsString = predictedHelp + "\n\n" + acceptedHelp + "\n\n" +
    outputHelp;

  // Initialize flags without parameters string.
  string exactHelp = "\
-e, -E, --exact\n\
     Specifies exact comparisons (slippage is not allowed).\
";

  string helpLabel = "\
-h, -H, --help\n\
     Displays this usage details message.\
";

  string printHelp = "\
-p, -P, --print\n\
     Prints the output file to standard output.\n\
     This won't override the default behavior of writing to a file.\
";

  string flagsNoParamsString = exactHelp + "\n\n" + helpLabel + "\n\n" +
    printHelp;

  // Initialize the flags with parameters string.
  string flagsParamsString = "";

  // Initialize parser and check the command line.
  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 3,
						   6, flagsNoParams,
						   0, NULL );
  parser->setUsageStrings( usageString.c_str(), paramsString.c_str(),
			   flagsNoParamsString.c_str(),
			   flagsParamsString.c_str() );
  bool noError = parser->checkLine();

  // If command line structure is OK, proceed with argument retrieval, parsing
  if( noError ) {

    // Get the required parameters.
    predicted = parser->getParameter( 1 );
    accepted = parser->getParameter( 2 );
    output = parser->getParameter( 3 );

    // Set non-parameterized boolean flags.
    exact = parser->containsInGroup( "-e -E --exact" );
    print = parser->containsInGroup( "-p -P --print" );
  }

  // Delete the parser and return the error state.
  delete parser;
  return noError;

}

///////////////////////////////////////////////////////////////////////////////
// Run scoring calculations.
///////////////////////////////////////////////////////////////////////////////
void Scorer_Interface::run() {

  /*
   * Create a variable to track errors.
   * Throughout, the calculation proceeds as long as error = 0.
   */
  int error = 0;

  /*
   * Create two RNA strands, one for the predicted structure and one for the
   * accepted structure.
   * For both, specify type = 1 (ct file).
   * The type of nucleic acid does not matter, so it isn't specified.
   */

  // Show message that intialization has begun.
  cout << "Initializing predicted and accepted structures..." << flush;

  // Create the first RNA strand.
  RNA* strand = new RNA( predicted.c_str(), 1 );
  ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
  error = checker->isErrorStatus();

  // Create the second RNA strand.
  RNA* strand2 = new RNA( accepted.c_str(), 1 );
  ErrorChecker<RNA>* checker2 = new ErrorChecker<RNA>( strand2 );
  if( error != 0 ) { error = checker2->isErrorStatus(); }

  // Show message that initialization is finished.
  if( error == 0 ) { cout << "done." << endl; }

  /*
   * Check both structures to make sure that only one structure is present in
   * each strand.
   */

  if( error == 0 ) {

    // Get the number of structures in each strand.
    int number1 = strand->GetStructureNumber();
    int number2 = strand2->GetStructureNumber();

    // If the number of structures in either strand is incorrect, show the
    // appropriate error message.
    const char* numError = " ct must contain only one structure.";
    if( number1 != 1 ) { cerr << "Predicted" << numError << endl; }
    if( number2 != 1 ) { cerr << "Accepted" << numError << endl; }
    if( ( number1 != 1 ) || ( number2 != 1 ) ) { error = 1; }
  }

  /*
   * Check the length of the structures to make sure they are equal.
   */
  if( error == 0 ) {

    // Get the sequence lengths.
    int length1 = strand->GetSequenceLength();
    int length2 = strand2->GetSequenceLength();

    // If the lengths aren't equal, show an error message.
    if( length1 != length2 ) {
      cerr << "Predicted and accepted structures are not the same length."
	   << endl;
      error = 1;
    }

  }

  /*
   * Calculate sensitivity and PPV.
   */

  // Create the sensitivity and PPV string streams.
  stringstream sensitivity( stringstream::in | stringstream::out );
  stringstream ppv( stringstream::in | stringstream::out );

  // If no error has previously occurred, do calculations.
  if( error == 0 ) {

    // Print message saying sensitivity and PPV are being calculated.
    cout << "Calculating sensitivity and PPV..." << flush;

    // Get the predicted and accepted structures that back the RNA strand.
    structure* predictedBack = strand->GetStructure();
    structure* acceptedBack = strand2->GetStructure();

    // Calculate sensitivity.
    int pairs1 = 0;
    int score1 = 0;
    scorer( acceptedBack, predictedBack, &score1, &pairs1, 1, exact );
    double percent1 = ( ( (double)score1 ) / ( (double)pairs1 ) ) * 100;

    // Calculate PPV.
    int pairs2 = 0;
    int score2 = 0;
    scorerppv( acceptedBack, predictedBack, &score2, &pairs2, 1, exact );
    double percent2 = ( ( (double)score2 ) / ( (double)pairs2 ) ) * 100;

    // Fill the string streams that hold sensitivity and PPV data.
    const char* ss = "Sensitivity: ";
    const char* pp = "PPV:         ";

    sensitivity << ss << score1 << " / " << pairs1 << " = " << fixed
		<< setprecision( 2 ) << percent1 << "%";
    ppv << pp << score2 << " / " << pairs2 << " = " << fixed
	<< setprecision( 2 ) << percent2 << "%";

    // Print message saying sensitivity and PPV are done.
    cout << "done." << endl;
  }

  /*
   * Output the sensitivity and PPV data to a file, and to the screen, if the
   * user requested it.
   */

  if( error == 0 ) {

    // Write the sensitivity and PPV to a file stream.
    ofstream out( output.c_str() );
    out << "Accepted Structure:  " << accepted << endl
	<< "Predicted Structure: " << predicted << endl
	<< sensitivity.str() << endl
	<< ppv.str()
	<< endl;
    out.close();

    // If the user wants data output to the screen, print to standard output.
    if( print ) {
      cout << endl
	   << "Accepted Structure:  " << accepted << endl
	   << "Predicted Structure: " << predicted << endl
	   << sensitivity.str() << endl
	   << ppv.str() << endl
	   << "Saved in output file: " << output
	   << endl << endl;
    }
  }

  /*
   * Clean up the strands and error checkers.
   */

  // Clean up strand and error checker 1.
  delete strand;
  delete checker;

  // Clean up strand and error checker 2.
  delete strand2;
  delete checker2;

  /*
   * Print out a confirmation of the run finishing.
   */

  // Compose the final message.
  string message = "Nucleic acid structure scoring done";
  if( error != 0 ) { message += " with errors"; }
  message += ".";

  // Print the message to either standard output or standard error, depending
  // on if an error occurred at some point.
  if( error == 0 ) { cout << message << endl; }
  else { cerr << message << endl; }

}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the interface.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

  Scorer_Interface* score = new Scorer_Interface();
  bool parseable = score->parse( argc, argv );
  if( parseable ) { score->run(); }
  delete score;

}
