/*
 * A header file for a program that predicts structures composed of probable
 * base pairs and single-stranded nucleotides.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#ifndef MAX_ACCURACY_H
#define MAX_ACCURACY_H

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class MaxExpect {
 public:
  // Public constructor and methods

  /*
   * Name:        Constructor
   * Description: Initalizes all private variables
   */
  MaxExpect();

  /*
   * Name:        parse
   * Description: Parses command line arguments to determine what options are
   *              required for structure prediction.
   * Arguments:
   *     1.   The number of command line arguments
   *     2.   The command line arguments themselves
   * Returns:
   *     true if parsing completed without errors, false if not
   */
  bool parse( int argc, char** argv );

  /*
   * Name:        run
   * Description: Run structure prediction
   */
  void run();

 private:
  // Input and output file names
  char* input;             // The input pfs file
  char* output;            // The output ct file

  // Boolean flags
  bool isRNA;              // Is the strand RNA (true) or DNA (false)?
  bool isSequence;         // Is input a sequence file (true) or not (false)?

  // Structure parameters
  double gamma;            // The weight given to base pairs
  double percent;          // The percent energy difference
  int structures;          // The maximum number of structures
  int window;              // The window size
};

#endif /* MAX_ACCURACY_H */
