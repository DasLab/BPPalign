/*
 * A header file for a program that draws a structure and writes image output
 * to a Postscript image file.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#ifndef DRAWSTRUCTURE_H
#define DRAWSTRUCTURE_H

#include "../src/ParseCommandLine.h"
#include "../src/PostscriptWriting/Postscript_Wrapper.h"

class DrawStructure {
 public:
  // Public constructor and methods.

  /*
   * Name:        Constructor.
   * Description: Initializes all private variables.
   */
  DrawStructure();

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
  // Private variables.

  // Description of the calculation type.
  string calcType;

  // Input and output file names.
  string inputFile;              // The input structure ct file.
  string outputFile;             // The output Postscript image file.

  string probabilityFile;        // The optional probability annotation file.
  string SHAPEFile;              // The optional SHAPE annotation file.

  // Draw a circular structure (true) or not (false)?
  bool circular;
};

#endif /* DRAWSTRUCTURE_H */
