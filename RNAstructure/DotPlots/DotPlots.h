/*
 * A header file for a program that calculates a dot plot and writes image
 * output to either a Postscript image file or a dot plot text file.
 * This class can also write a dot plot text file.
 *
 * This class can read three types of files:
 *      1.  Partition function save files (.pfs)
 *      2.  Folding save files (.sav)
 *      3.  Dynalign save files (.ds
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#ifndef DOT_PLOTS_H
#define DOT_PLOTS_H

#include "../src/DotPlotHandler.h"
#include "../src/ParseCommandLine.h"
#include "../src/PostscriptWriting/Postscript_Wrapper.h"

class DotPlots {
 public:
  // Public constructor and methods.

  /*
   * Name:        Constructor.
   * Description: Initializes all private variables.
   */
  DotPlots();

  /*
   * Name:        parse
   * Description: Parses command line arguments to determine what options are
   *              required for a particular calculation.
   * Arguments:
   *     1.   The number of command line arguments.
   *     2.   The command line arguments themselves.
   * Returns:
   *     True if parsing completed without errors, false if not.
   */
  bool parse( int argc, char** argv );

  /*
   * Name:        run
   * Description: Run calculations.
   */
  void run();

 private:
  // Private executable identifier method.

  /*
   * Name:        identify
   * Description: Identify the name of the executable, used to distinguish a
   *              requested dot plot type.
   * Arguments:
   *     1. executableName
   *        The raw executable name, possibly with paths.
   * Returns:
   *     The name of the executable.
   */
  string identify( string executableName );

 private:
  // Private variables.

  // Description of the calculation type.
  string calcType;

  // Dot plot handler.
  DotPlotHandler* plotMaker;

  // Boolean flag signifying whether a simple text file (true) or a complex
  // Postscript image (false) is written.
  bool simpleOut;
};

#endif /* DOT_PLOTS_H */
