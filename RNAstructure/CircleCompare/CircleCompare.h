/*
 * A header file for a class that takes two CT files, then compares their
 * pairings and outputs the data graphically as a Postscript file.
 *
 * (c) 2010 Mathews Lab, University of Rochester
 * Written by Jessica Reuter
 */

#ifndef CIRCLECOMPARE_H
#define CIRCLECOMAPRE_H

#include <iomanip>
#include <sstream>
#include <stdlib.h>

#include "../src/PostscriptWriting/Postscript_Circle_Structures.h"
#include "../src/PostscriptWriting/Postscript_Probability_Legend.h"
#include "../src/PostscriptWriting/Postscript_SHAPE_Legend.h"
#include "../src/score.h"

class CircleCompare {
 public:
  // Public constructor and methods.

  /*
   * Name:        Constructor.
   * Description: Initializes all private variables.
   * Arguments:
   *     1. i1
   *        The name of the predicted CT structure file.
   *     2. i2
   *        The name of the accepted CT structure file.
   *     3. o
   *        The name of the Postscript output file to write the image to.
   */
  CircleCompare( string i1, string i2, string o );

  /*
   * Name:        writeCirclePlot
   * Description: Using the component structure files and analyzed pairings,
   *              write the final Postscript drawing with color-coded compared
   *              base pairs.
   * Arguments:
   *     1. number
   *        The structure number to write.
   *     2. alternative
   *        True if alternative color scheme is used.
   *        Default is false.
   *     3. exact
   *        True if slippage is not allowed.
   *        Default is false (slippage allowed).
   *     4. probability
   *        The name of a file that holds probability annotation.
   *        Default is "", or no annotation.
   *     5. shape
   *        The name of a file that holds SHAPE annotation.
   *        Default is "", or no annotation.
   * Returns:
   *     True if the plot was written correctly, false if not.
   */
  bool writeCirclePlot( int number,
			bool alternative = false, bool exact = false,
			string probability = "", string shape = "" );

 private:
  // Private variables.

  // Input and output file names.
  string input1;                 // Predicted structure input file name.
  string input2;                 // Actual structure input file name.
  string output;                 // Circular plot output file name.

  // The number of elements in the Postscript color schemes.
  // For each scheme, there are three rows, each row consisting of a Postscript
  // color string, followed by the actual name of that color, with necessary
  // trailing spaces for alignment added.
  // Schemes are organized as:
  //       0: Color for pair in both predicted and accepted structure
  //       1: Color for pair in predicted structure only
  //       2: Color for pair in accepted structure only
  static const int COLORS_Y = 3;
  static const int COLORS_X = 2;

  // Postscript color schemes.
  static const string schemeDefault[COLORS_Y][COLORS_X];
  static const string schemeAlternative[COLORS_Y][COLORS_X];
};

#endif /* CIRCLECOMPARE_H */
