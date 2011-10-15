/*
 * An implementation for a class that takes two CT files, then compares their
 * pairings and outputs the data graphically as a Postscript file.
 *
 * (c) 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#include "CircleCompare.h"

///////////////////////////////////////////////////////////////////////////////
// Arrays of Postscript colors for the default and alternative color schemes.
///////////////////////////////////////////////////////////////////////////////

// Default color scheme.
const string CircleCompare::schemeDefault[COLORS_Y][COLORS_X] = {
  { Postscript_Constants::GREEN, "Green" },
  { Postscript_Constants::RED,   "Red" },
  { Postscript_Constants::BLACK, "Black" }
};

// Alternative color scheme.
const string CircleCompare::schemeAlternative[COLORS_Y][COLORS_X] = {
  { Postscript_Constants::GREEN,  "Green" },
  { Postscript_Constants::PURPLE, "Purple" },
  { Postscript_Constants::RED,    "Red" }
};

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
CircleCompare::CircleCompare( string i1, string i2, string o ) {

  // Set the input and output file names.
  input1 = i1;
  input2 = i2;
  output = o;

}

///////////////////////////////////////////////////////////////////////////////
// Write the final output file with the circle plot.
///////////////////////////////////////////////////////////////////////////////
bool CircleCompare::writeCirclePlot( int number, bool alternative, bool exact,
				     string probability, string shape ) {

  // Choose the color scheme.
  string bothColor = ( !alternative ) ?
    schemeDefault[0][0] : schemeAlternative[0][0];
  string bothName = ( !alternative ) ?
    schemeDefault[0][1] : schemeAlternative[0][1];
  string oneColor = ( !alternative ) ?
    schemeDefault[1][0] : schemeAlternative[1][0];
  string oneName = ( !alternative ) ?
    schemeDefault[1][1] : schemeAlternative[1][1];
  string twoColor = ( !alternative ) ?
    schemeDefault[2][0] : schemeAlternative[2][0];
  string twoName = ( !alternative ) ?
    schemeDefault[2][1] : schemeAlternative[2][1];

  // Create the combined circle structure.
  Postscript_Circle_Structures* compare =
    new Postscript_Circle_Structures( input2, output );
  compare->setColors(
		     bothColor, bothName,
		     twoColor, twoName,
		     oneColor, oneName
		     );
  compare->setIdentifiers( "Accepted Structure", "Predicted Structure" );
  compare->parsePairs();
  compare->addSecondStructure( input1, number );

  // Rebuild the accepted structure, extract the structure and pairings.
  Postscript_Circle_Structure* acceptedStructure =
    new Postscript_Circle_Structure( input2, "" );
  acceptedStructure->parsePairs();
  int acceptedPairs = acceptedStructure->getNumPairs();
  structure* acceptedStructureBack =
    acceptedStructure->getStrand()->GetStructure();

  // Rebuild the predicted structure, extract the structure and pairings.
  Postscript_Circle_Structure* predictedStructure =
    new Postscript_Circle_Structure( input1, "" );
  predictedStructure->parsePairs( number );
  int predictedPairs = predictedStructure->getNumPairs();
  structure* predictedStructureBack =
    predictedStructure->getStrand()->GetStructure();

  // Add annotation to the predicted structure, if necessary.
  if( ( probability != "" ) || ( shape != "" ) ) {
    if( probability != "" ) {
      compare->addProbabilityAnnotation( probability );
      compare->setPairingsAndNucleotidesAnnotated();
    } else if( shape != "" ) {
      compare->addSHAPEAnnotation( shape );
      compare->setPairingsAndNucleotidesAnnotated();
    }

    // If an error occurred with annotation, return.
    if( compare->isError() ) {
      delete predictedStructure;
      delete acceptedStructure;
      delete compare;
      return false;
    }
  }

  // Calculate sensivity.
  int pairs1 = 0;
  int score1 = 0;
  scorer( acceptedStructureBack, predictedStructureBack, &score1, &pairs1,
	  number, exact );
  double percent1 = ( ( (double)score1 ) / ( (double)pairs1 ) ) * 100;

  stringstream sensitivity( stringstream::in | stringstream::out );
  sensitivity << score1 << " / " << pairs1 << " = " << fixed
	      << setprecision( 2 ) << percent1 << "%";
  string sens = sensitivity.str();

  // Calculate PPV.
  int pairs2 = 0;
  int score2 = 0;
  scorerppv( acceptedStructureBack, predictedStructureBack, &score2, &pairs2,
	     number, exact );
  double percent2 = ( ( (double)score2 ) / ( (double)pairs2 ) ) * 100;

  stringstream ppvStream( stringstream::in | stringstream::out );
  ppvStream << score2 << " / " << pairs2 << " = " << fixed
            << setprecision( 2 ) << percent2 << "%";
  string ppv = ppvStream.str();

  // Calculate the number of pseudoknotted pairs in each structure.
  predictedStructure->calculatePseudoknottedPairs();
  int predictedPseudoPairs = predictedStructure->getNumPseudoknottedPairs();

  acceptedStructure->calculatePseudoknottedPairs();
  int acceptedPseudoPairs = acceptedStructure->getNumPseudoknottedPairs();

  // Put numbers for accepted and predicted pairs in string streams.
  stringstream right2Stream( stringstream::in | stringstream::out );
  right2Stream << "Pairs: " << acceptedPairs << ")";
  string accPairs = right2Stream.str();

  stringstream right3Stream( stringstream::in | stringstream::out );
  right3Stream << "Pseudoknotted Pairs: " << acceptedPseudoPairs << ")";
  string accPseudo = right3Stream.str();

  stringstream predStream1( stringstream::in | stringstream::out );
  predStream1 << "(Pairs: " << predictedPairs << ")";
  string predPairs = predStream1.str();

  stringstream predStream2( stringstream::in | stringstream::out );
  predStream2 << "(Pseudoknotted Pairs: " << predictedPseudoPairs << ")";
  string predPseudo = predStream2.str();

  // Create the added data array, which is the statistics labels.
  vector<string> data;
  data.push_back( "/Courier-Bold findfont fontSize 2 div scalefont setfont" );
  data.push_back( "" );
  data.push_back( "/x1 570 (Accepted:) stringwidth pop sub def" );
  data.push_back( "/x2 570 (" + accPairs + " stringwidth pop sub def" );
  data.push_back( "/x3 570 (" + accPseudo + " stringwidth pop sub def" );
  data.push_back( "/x4 570 (Sensitivity: " + sens +
		  ") stringwidth pop sub def" );
  data.push_back( "/x5 570 (PPV: " + ppv + ") stringwidth pop sub def" );
  data.push_back( "" );
  data.push_back( "/statistics [" );
  data.push_back( "    [(Predicted:) 40 70]" );
  data.push_back( "    [" + predPairs + " 40 55]" );
  data.push_back( "    [" + predPseudo + " 40 40]" );
  data.push_back( "    [(Accepted:) x1 70]" );
  data.push_back( "    [(" + accPairs + " x2 55]" );
  data.push_back( "    [(" + accPseudo + " x3 40]" );
  data.push_back( "    [(Sensitivity: " + sens + ") x4 595]" );
  data.push_back( "    [(PPV: " + ppv + ") x5 580]" );
  data.push_back( "] def" );
  data.push_back( "" );
  data.push_back( "statistics { aload pop moveto show } forall" );

  // Set the additional data and write the structure.
  compare->setAdditionalData( data );
  compare->writeOutput( true );

  // Delete the structures.
  delete predictedStructure;
  delete acceptedStructure;
  delete compare;

  // Return true, signifying that the plot was written correctly.
  return true;

}
