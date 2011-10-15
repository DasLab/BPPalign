/*
 * A program that calculates a dot plot and writes image output to either a
 * Postscript file or a dot plot text file.
 * This class can also write a dot plot text file.
 *
 * This class can read three types of files:
 *      1.  Partition function save files (.pfs)
 *      2.  Folding save files (.sav)
 *      3.  Dynalign save files (.dsv)
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "DotPlots.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
DotPlots::DotPlots() {

  // Initialize the calculation type description.
  calcType = "Dot plot output";

  // Initialize initial output type to Postscript image.
  simpleOut = false;

}

///////////////////////////////////////////////////////////////////////////////
// Identify the simple name of the executable.
///////////////////////////////////////////////////////////////////////////////
string DotPlots::identify( string executableName ) {

  // Get the length of the entire path to the executable.
  int length = executableName.length();

  // Find where the executable name itself starts.
  int begin = executableName.find_last_of( "/\\" ) + 1;
  if( begin < 0 || begin > length ) { begin = 0; }

  // Slice the executable name off the path.
  string next( executableName.substr( begin, length ) );
  int length2 = next.length();

  // Chop off any extension the executable name has.
  // (This step is basically for Windows, which tacks ".exe" on executables.)
  int end = next.find( "." );
  if( end < 0 || end > length2 ) { end = length2; }

  // Return the executable name.
  return next.substr( 0, end );

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool DotPlots::parse( int argc, char* argv[] ) {

  // Identify the type name using argv[0].
  string type = identify( argv[0] );

  // Initialize array of flags without parameters. The contents are the same
  // regardless of the dot plot type, but only certain options are used with
  // different plot types.
  const char* flagsNoParams[] = {

    // Used for all plots.
    "-t", "-T", "--text",

    // Used in DynalignDotPlot only.
    "-s2", "-S2", "--sequence2"
  };

  // Initialize the array of flags with parameters.
  const char* flagsParams[] = {
    "-e", "-E", "--entries",
    "-min", "-MIN", "--minimum",
    "-max", "-MAX", "--maximum"
  };

  // Initialize usage string.
  string usageString = 
    "USAGE: " +  type +  " <save file> <ps file> [options]";

  // Initialize required parameters string.
  string inputString = "\
<save file>\n\
     The name of a file from which dot plot data will be read.\n\
     This input file is always some type of binary save file.\
";

  string outString = "\
<ps file>\n\
     The name of a Postscript image file to which output will be written.\
";

  string paramsString =
     inputString + "\n\n" + outString;

  // Initialize flags without parameters string.
  string helpLabel = "\
-h, -H, --help\n\
     Display the usage details message.\
";

  string seq2String = "";
  if( type == "DynalignDotPlot" ) {
    seq2String = "\
-s2, -S2, --sequence2\n\
     Specifies that the dot plot should be the second sequence.\n\
     If no sequence is specified, the plot is the first sequence.\
";
  }

  string textString = "\
-t, -T, --text\n\
     Specifies that output should be a dot plot (text) file.\
";

  string flagsNoParamsString = helpLabel + "\n\n";
  if( type == "DynalignDotPlot" ) {
    flagsNoParamsString = flagsNoParamsString + seq2String + "\n\n";
  }
  flagsNoParamsString += textString;

  // Initialize flags with parameters string.
  stringstream entryStream;
  entryStream << "\
-e, -E, --entries\n\
     Specifies the number of colors in the dot plot.\n\
     Default is 5 colors.\n\
     Minimum is " << MIN_COLORS << " colors.\n\
     Maximum is " << MAX_COLORS << " colors.\
";

  string maxString = "\
-max, -MAX, --maximum\n\
     Specifies the maximum value that is viewable in the plot.\n\
     Default is the largest allowable point in a given data file.\n\
     If the given value is greater than the default, it is ignored.\
";

  string minString = "\
-min, -MIN, --minimum\n\
     Specifies the minimum value that is viewable in the plot.\n\
     Default is the smallest allowable point in a given data file.\n\
     If the given value is less than the default, it is ignored.\
";

  string flagsParamsString =
    entryStream.str() + "\n\n" + maxString + "\n\n" + minString;

  // Initialize parser and check the command line.
  int numFlagsNoParams = ( type != "DynalignDotPlot" ) ? 3 : 6;

  ParseCommandLine* parser = new ParseCommandLine( argc, argv, 2,
                                                   numFlagsNoParams,
						   flagsNoParams,
                                                   9, flagsParams );
  parser->setUsageStrings( usageString.c_str(), paramsString.c_str(),
                           flagsNoParamsString.c_str(),
			   flagsParamsString.c_str() );
  bool noError = parser->checkLine();

  // Set the plot type. If it isn't valid, stop parsing right away.
  char plotType;
  if( type == "DynalignDotPlot" ) {
    plotType = ( !parser->containsInGroup( "-s2 -S2 --sequence2" ) ) ?
      TYPE_DYNALIGN1 : TYPE_DYNALIGN2;
  } else if( type == "EnergyPlot" ) {
    plotType = TYPE_ENERGY;
  } else if( type == "ProbabilityPlot" ) {
    plotType = TYPE_PROBABILITY;
  } else {
    plotType = TYPE_UNDEFINED;
    noError = false;
  }

  // Initialize the Dynalign sequence number, if it's necessary later.
  int dynalignSequence = 0;

  // If command line structure is OK, proceed with argument retrieval, parsing.
  if( noError ) {

    // Get required parameters.
    string inputFile = parser->getParameter( 1 );
    string outputFile = parser->getParameter( 2 );

    // Set text option.
    simpleOut = parser->containsInGroup( "-t -T --text" );

    // Set entries option. Since there is a range of acceptable values, the
    // option is set without bounds, and then the bounds are explicitly
    // checked later.
    int entries = 5;
    bool entriesSet = parser->setOption( entries, "-e -E --entries",
					 "Invalid plot entries value given." );
    if( entriesSet ) {
      if( entries < MIN_COLORS ) {
	cerr << "Too few plot entries given." << endl;
	noError = false;
      } else if( entries > MAX_COLORS ) {
	cerr << "Too many plot entries given." << endl;
	noError = false;
      }
    } else { noError = false; }

    // Set the minimum bound of the plot.
    double minBound = numeric_limits<double>::infinity() * -1;
    bool minSet = parser->setOption( minBound, "-min -MIN --minimum",
				     "Invalid minimum plot value given." );
    if( minSet == false ) { noError = false; }

    // Set the maximum bound of the plot.
    double maxBound = numeric_limits<double>::infinity();
    bool maxSet = parser->setOption( maxBound, "-max -MAX --maximum",
				     "Invalid maximum plot value given." );
    if( maxSet == false ) { noError = false; }

    // Check that the maximum is greater than or equal to the minimum.
    if( minBound > maxBound ) {
      cerr << "Minimum plot value cannot be greater than maximum plot value."
	   << endl;
      noError = false;
    }

    // If parsing happened successfully, create the back end dot plot handler.
    if( noError ) {

      // Initialize the dot plot handler with correct values.
      plotMaker = new DotPlotHandler( inputFile, outputFile );
      plotMaker->setEntries( entries );
      plotMaker->setMinimum( minBound );
      plotMaker->setMaximum( maxBound );
      plotMaker->setPlotType( plotType );
    }
  }

  // Delete the parser and return the error state.
  delete parser;
  return noError;

}

///////////////////////////////////////////////////////////////////////////////
// Run dot plot creation.
///////////////////////////////////////////////////////////////////////////////
void DotPlots::run() {

  // If simple text output was requested, write a text dot plot file.
  if( simpleOut ) {

    // Show a message saying that a dot plot text file is being made.
    cout << "Reading dot plot data..." << flush;

    // Depending on the type of dot plot, get its data.
    if( plotMaker->getPlotType() == TYPE_DYNALIGN1 ) {
      plotMaker->readDynalignSeq1Data();
    } else if( plotMaker->getPlotType() == TYPE_DYNALIGN2 ) {
      plotMaker->readDynalignSeq2Data();
    } else if( plotMaker->getPlotType() == TYPE_ENERGY ) {
      plotMaker->readFoldingData();
    } else if( plotMaker->getPlotType() == TYPE_PROBABILITY ) {
      plotMaker->readPartitionData();
    } else {
      plotMaker->setPlotType( TYPE_UNDEFINED );
    }

    // If an error occurred, show an error message. If an error didn't occur,
    // show that data reading is done.
    if( ( plotMaker->isError() ) ||
	( plotMaker->getPlotType() == TYPE_UNDEFINED ) ) {
      cerr << "Error reading dot plot text file data." << endl;
      return;
    } else { cout << "done." << endl; }

    // If an error didn't happen, write the text file.
    cout << "Writing dot plot text file..." << flush;
    plotMaker->writePlotFile();
    if( !plotMaker->isError() ) { cout << "done." << endl; }
    else { cerr << "Error writing dot plot text file." << endl; }

    // Write a finishing message depending on the error state.
    if( !plotMaker->isError() ) {
      cout << "Dot plot text file writing complete." << endl;
    } else {
      cerr << "Dot plot text file writing complete with errors." << endl;
    }
  }

  // Otherwise, generate a Postscript image of the plot.
  else {

    // Show a message that Postscript dot plot has been started.
    cout << "Creating Postscript dot plot image..." << flush;

    // Create the Postscript wrapper.
    Postscript_Wrapper* wrap = new Postscript_Wrapper();
    bool noError = true;

    // Depending on the plot type, draw a plot.
    if( plotMaker->getPlotType() == TYPE_DYNALIGN1 ) {
      noError = wrap->plotDynalign1( plotMaker );
    } else if( plotMaker->getPlotType() == TYPE_DYNALIGN2 ) {
      noError = wrap->plotDynalign2( plotMaker );
    } else if( plotMaker->getPlotType() == TYPE_ENERGY ) {
      noError = wrap->plotEnergy( plotMaker );
    } else if( plotMaker->getPlotType() == TYPE_PROBABILITY ) {
      noError = wrap->plotProbability( plotMaker );
    } else {
      plotMaker->setPlotType( TYPE_UNDEFINED );
      noError = false;
    }

    // If an error occurred, show an error message. If not, show it was done
    // successfully.
    if( noError ) { cout << "done." << endl; }
    else { cerr << "Error drawing Postscript dot plot." << endl; }

    // Delete the Postscript wrapper.
    delete wrap;

    // Write a finishing message depending on the error state.
    if(!plotMaker->isError() && noError ) {
      cout << calcType << " complete." << endl;
    } else {
      cerr << calcType << " complete with errors." << endl;
    }
  }

  // Delete the dot plot handler after everything is done.
  delete plotMaker;

}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

  DotPlots* runner = new DotPlots();
  bool parseable = runner->parse( argc, argv );
  if( parseable == true ) { runner->run(); }
  delete runner;

}
