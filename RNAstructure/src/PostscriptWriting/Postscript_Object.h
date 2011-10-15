/*
 * A header file for a class that handles an object written in Postscript.
 *
 * (c) 2010  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#ifndef POSTSCRIPT_OBJECT_H
#define POSTSCRIPT_OBJECT_H

#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "../../RNA_class/RNA.h"
#include "../ErrorChecker.h"

// A small namespace used to handle Postscript constants.
namespace Postscript_Constants {

  // Declarations of strand types
  const int CT_TYPE  = 1;
  const int SEQ_TYPE = 2;
  const int PFS_TYPE = 3;
  const int SAV_TYPE = 4;

  // Declarations of truncation modes for strings
  const int NO_TRUNC    = 0;
  const int TRUNC_BEGIN = 1;
  const int TRUNC_END   = 2;

  // Declaration of string that serves as a tab for Postscript writing
  const string tab = "    ";

  // Declarations of constants which hold Postscript colors
  const string BLACK        = "0.00 0.00 0.00";
  const string BLUE         = "0.00 0.00 1.00";
  const string BRIGHT_GREEN = "0.00 1.00 0.00";
  const string DARK_YELLOW  = "0.83 0.83 0.17";
  const string DARK_PINK    = "1.00 0.50 1.00";
  const string GRAY         = "0.67 0.67 0.67";
  const string GREEN        = "0.00 0.50 0.00";
  const string LIGHT_BLUE   = "0.00 0.67 1.00";
  const string ORANGE       = "1.00 0.50 0.00";
  const string PURPLE       = "0.40 0.00 0.60";
  const string RED          = "1.00 0.00 0.00";
  const string WHITE        = "1.00 1.00 1.00";
};

// Namespace usage statements
using namespace std;
using namespace Postscript_Constants;

class Postscript_Object {
 public:
  // Public error checking method

  /*
   * Name:        isError
   * Description: Tells if this object has encountered an error.
   * Returns:
   *     True if an error has been encountered, false if not.
   */
  bool isError();

 protected:
  // Protected constructor, destructor, and trimming method

  /*
   * Name:        Constructor
   * Descripton:  Initalizes private variables.
   */
  Postscript_Object();

  /*
   * Name:        Destructor
   * Description: Deletes all dynamically allocated variables.
   */
  ~Postscript_Object();

  /*
   * Name:        escapeAndTrim
   * Description: Escape parentheses characters "(" and ")" to ensure a string
   *              is safe to print in Postscript. In this function, "trim" has
   *              two meanings. First, it removes excess whitespace from the
   *              beginning and end of the string. Second, if the string is
   *              long enough where it will flow off the page, the string is
   *              truncated.
   *              Mode TRUNC_BEGIN should be used for file paths; it truncates
   *              a string like the following example:
   *              /dir1/dir2/dir3/dir4/dir5/file.txt -> PATH/file.txt
   *              Mode TRUNC_END should be used for any other string type; it
   *              truncates like the following example:
                  This is a string to truncate -> This is a string...
   * Arguments:
   *     1. data
   *        The string to be escaped.
   *     2. mode
   *        An int that tells what kind of truncation to do.
   *        Default is 0/NO_TRUNC (no truncation).
   * Returns:
   *     A string which holds the string with proper escapes
   */
  string escapeAndTrim( string data, int node = NO_TRUNC );

 protected:
  // Protected boolean that tells if an error occurred at any time.
  bool error;
};

#endif /* POSTSCRIPT_OBJECT_H */
