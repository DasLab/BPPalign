/*
 * A header file for a program that folds a strand of nucleic acids.
 * These nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#ifndef FOLD_H
#define FOLD_H

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class Fold {
 public:
  // Public constructor and methods.

  /*
   * Name:        Constructor.
   * Description: Initializes all private variables.
   */
  Fold();

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
  string seqFile;           // The input sequence file.
  string ctFile;            // The output ct file.
  string saveFile;          // The optional output save file.

  string constraintFile;    // The optional folding constraints file.
  string SHAPEFile;         // The optional SHAPE constraints file.

  string singleOffsetFile;  // The optional single strand offset file.
  string doubleOffsetFile;  // The optional double strand offset file.

  // The intercept for SHAPE constraints.
  double intercept;

  // The intercept for single-stranded SHAPE constraints.
  double interceptSingle;

  // Flag signifying if calculation handles RNA (true) or DNA (false).
  bool isRNA;

  // The maximum pairing distance.
  int maxDistance;

  // The maximum internal bulge loop size.
  int maxLoop;

  // The maximum number of structures to generate.
  int maxStructures;

  // The maximum percent energy difference.
  double percent;

  // The slope for SHAPE constraints.
  double slope;

  // The slope for single-stranded SHAPE constraints.
  double slopeSingle;

  // The temperature at which calculation occurs.
  double temperature;

  // The window size for calculation.
  int windowSize;
};

#endif /* FOLD_H */
