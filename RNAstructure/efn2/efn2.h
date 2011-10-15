/*
 * A header file for a program that calculates the free energy of a strand of
 * nucleic acids.
 * These nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2008  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#ifndef EFN2_H
#define EFN2_H

#include <iomanip>
#include <vector>

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class efn2Interface {
 public:
  // Public constructor and methods.

  /*
   * Name:        Constructor.
   * Description: Initializes all private variables.
   */
  efn2Interface();

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
  string ctFile;                 // The input ct file.
  string outFile;                // The output energy text file.

  // Flag signifying if calculation handles RNA (true) or DNA (false).
  bool isRNA;

  // Flag signifying whether the output should be piped to standard output as
  // well as written (true) or not (false).
  bool stdPrint;

  // The temperature at which calculation occurs.
  double temperature;

  // Flag signifying whether a thermodynamic details file should be written
  // (true) or not (false).
  bool writeTherm;
};

#endif /* EFN2_H */
