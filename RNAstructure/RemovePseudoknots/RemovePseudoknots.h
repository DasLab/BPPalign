/*
 * A header file for a program that breaks pseudoknots to allow for easier
 * analysis of a strand of nucleic acids.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#ifndef REMOVE_PSEUDOKNOTS_H
#define REMOVE_PSEUDOKNOTS_H

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class RemovePseudoknots {
 public:
  // Public constructor and methods.

  /*
   * Name:        Constructor.
   * Description: Initializes all private variables.
   */
  RemovePseudoknots();

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
  string input;                  // The input ct file.
  string output;                 // The output ct file.

  // Flag signifying if calculation handles RNA (true) or DNA (false).
  bool isRNA;

  // The temperature at which calculation occurs.
  double temperature;
};

#endif /* REMOVE_PSEUDOKNOTS_H */
