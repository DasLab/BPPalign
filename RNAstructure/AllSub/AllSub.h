/*
 * A header file for a program that finds all suboptimal structures within a
 * predefined small increment for a strand of nucleic acids.
 * These nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#ifndef ALLSUB_H
#define ALLSUB_H

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class AllSub {
 public:
  // Public constructor and methods.

  /*
   * Name:        Constructor.
   * Description: Initializes all private variables.
   */
  AllSub();

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
  string seqFile;                // The input sequence file.
  string ctFile;                 // The output ct file.
  string constraintFile;         // The optional input constraints file.

  // The absolute energy difference.
  double absolute;

  // Flag signifying if calculation handles RNA (true) or DNA (false).
  bool isRNA;

  // The maximum percent energy difference.
  double percent;

  // The temperature at which calculation occurs.
  double temperature;
};

#endif /* ALLSUB_H */
