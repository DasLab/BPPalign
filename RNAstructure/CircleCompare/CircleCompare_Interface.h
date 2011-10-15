/*
 * A header file for a program that runs CircleCompare to compare two different
 * CT file structures.
 *
 * (c) 2010  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#ifndef CIRCLECOMPARE_INTERFACE_H
#define CIRCLECOMPARE_INTERFACE_H

#include "../src/ParseCommandLine.h"
#include "CircleCompare.h"

class CircleCompare_Interface {
 public:
  // Public constructor and methods.

  /*
   * Name:        Constructor.
   * Description: Initializes all private variables.
   */
  CircleCompare_Interface();

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

  /*
   * Name:        validate
   * Description: Check to make sure the CT files given are valid and are the
   *              same length, which is required for a successful comparison.
   * Returns:
   *     True if both CT files are valid and the same length, false if not.
   */
  bool validate();

 private:
  // Private variables.

  // Description of the calculation type.
  string calcType;

  // Input and output file names.
  string predicted;              // The input predicted ct file.
  string accepted;               // The input accepted ct file.
  string output;                 // The output Postscript image file.

  string probabilityFile;        // The optional probability annotation file.
  string SHAPEFile;              // The optional SHAPE annotation file.

  // Boolean flag signifying if the alternative color scheme is used (true)
  // or not (false).
  bool alternative;

  // Boolean flag signifying exact pairs in scoring (true) or not (false).
  bool exact;

  // Compare just the first prediction in the CT file (true) or not (false)?
  bool firstOnly;

  // The number of predicted structures to compare against.
  int predictedCount;
};

#endif /* CIRCLECOMPARE_INTERFACE_H */
