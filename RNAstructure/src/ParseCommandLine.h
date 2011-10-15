/*
 * A header file for a program that parses Unix-like command lines to determine
 * the data being given for input
 * (c) 2008 Jessica Reuter
 */

#ifndef PARSECOMMANDLINE_H
#define PARSECOMMANDLINE_H

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

using namespace std;

class ParseCommandLine {
 public:
  /*
   * Name:        Constructor
   * Description: Initalizes data in the parser
   * Arguments:
   *     1.   The number of command line arguments
   *     2.   The command line arguments themselves
   *     3.   The number of required parameters
   *     4.   The number of flags which don't take parameters
   *     5.   The flags (which don't take parameters) themselves
   *     6.   The number of flags which take parameters
   *     7.   The flags (which take parameters) themselves
   */
  ParseCommandLine( int argc, char* argv[], int parameters,
		    int numFlagsNoParams, const char* flagsNoParams[],
		    int numFlagsParams, const char* flagsParams[] );

  /*
   * Name:        Destructor
   * Description: Deletes the data array that held information about the
   *              command line during parsing.
   */
  ~ParseCommandLine();

  /*
   * Name:        checkLine
   * Description: Check the entire command line, for parameters and flags.
   * Argument:
   *     1. The help flags, as a string (see containsInGroup for description of
   *         string format); default flags are "-h", "-H", and "--help"
   * Returns:
   *     True if all components are present, false if not
   */
  bool checkLine( const char* helpFlags = "-h -H --help" );

  /*
   * Name:        contains
   * Description: Check if the command line contains a particular flag
   * Argument:
   *          The flag to check presence of
   * Returns:
   *     True if the command line contains the flag, false if not
   */
  bool contains( const char* flag );

  /*
   * Name:        containsInGroup
   * Description: Check if the command line contains one (or more) of a group
   *              of flags
   * Argument:
   *         A string holding the flags to check, separated by spaces
   * Returns:
   *     True if the command line contains one of the flags, false if not
   */
  bool containsInGroup( const char* flagGroup );

  /*
   * Name:        getOption
   * Description: Get an option that immediately follows a particular flag.
   * Argument:
   *     1.   The flag that precedes the option
   * Returns:
   *     the option
   */
  char* getOption( const char* flag );

  /*
   * Name:        getParameter
   * Description: Get a required parameter from the command line
   * Arguments:
   *     1.   The index (one-indexed) of the parameter
   * Returns:
   *     The parameter
   */
  char* getParameter( int index );

  /*
   * Name:        setOption
   * Description: Set an optional field from the command line.
   *              This method deals with char* only; a templated method exists
   *              for numbers.
   *              Note that if one flag in the flag group sets a variable, the
   *              other flags in the group will be ignored.
   * Arguments:
   *     1.   The optional field being set
   *     2.   The group of flags that mark the value to be set
   *     3.   An error message that is displayed if field setting fails
   *     4.   Boolean denoting whether the char* being set is a file name (true)
   *          or a simple char* value (false).
   * Returns:
   *     true if option was set successfully, false if not
   */
  bool setOption( const char*& variable, const char* flagGroup,
                  const char* message, bool file );

  /*
   * Name:        setUsageStrings
   * Description: Set the usage strings for a program that is using this parser.
   * Arguments:
   *     1. The usage string
   *     2. The string of required parameters
   *     2. The string of flags without parameters
   *     3. The string of flags with parameters
   */
  void setUsageStrings( const char* newUsageString,
			const char* paramString,
			const char* flagNoParamString,
			const char* flagParamString );

 public:
  // Templated method for setting an option pulled from the command line

  /*
   * Name:        setOption
   * Description: Set an optional field by its data from the command line. This
   *              method sets all types of numbers, and is a template function.
   *              Note that if one flag in the flag group sets a variable, the
   *              other flags in the group will be ignored.
   * Arguments:
   *     1.   The optional field being set
   *     2.   The group of flags that marks the value to be set
   *          Flags are separated by spaces.
   *     3.   An error message that is displayed if field setting fails
   *     4.   The cutoff type and value, as a string. The type and value must
   *          be separated by a single space, ie. ">= 0".
   *          Default value is "", signifying all values are valid
   * Returns:
   *     true if option was set successfully, false if not
   */
  template <class T>
  bool setOption( T& variable, const char* flagGroup, const char* message,
                  const char* cutoff = "" ) {
    char group[100];
    strcpy( group, flagGroup );

    char* flag = strtok( group, " " );
    while( flag != NULL ) {
      const char* safeFlag = referencePointer( flag );
      char* option = getOption( safeFlag );

      if( option != NULL ) {
        double numberTest = 0;
        istringstream stream( option );
	if( stream >> numberTest ) { variable = (T)atof( option ); }
        else {
          cerr << "Non-numeric input given for flag " << safeFlag << endl;
          return false;
	}

        if( strcmp( cutoff, "" ) != 0 ) {
          char cutoffArray[100];
          strcpy( cutoffArray, cutoff );

          char* cutoffType = strtok( cutoffArray, " " );
          double value = atof( strtok( NULL, " " ) );

          bool cutoffError =
            ( strcmp( cutoffType, "<" ) == 0 && variable >= value ) ||
            ( strcmp( cutoffType, "<=" ) == 0 && variable > value ) ||
            ( strcmp( cutoffType, "=" ) == 0 && variable != value ) ||
            ( strcmp( cutoffType, ">=" ) == 0 && variable < value ) ||
            ( strcmp( cutoffType, ">" ) == 0 && variable <= value );

          if( cutoffError == true ) {
            cerr << message << endl;
            return false;
          } else return true;
        }
      }

      flag = strtok( NULL, " " );
    }

    return true;
  }

 private:
  /*
   * Name:        referencePointer
   * Description: Ensure that a reference to a particular flag, parameter, or
   *              option connects with established data, rather than creates a
   *              completely new data point. This method is only used by the
   *              internal workings of ParseCommandLine, never by the end user.
   * Argument:
   *         The const char* whose reference must be checked.
   * Returns:
   *     The appropriate reference as a const char*
   */
  const char* referencePointer( const char* pointer );

  /*
   * Name:        usage
   * Description: A generic usage method that can be used for all programs which
   *              use this parser.
   */
  void usage();

 private:
  // Variables for parameter amounts
  int parameters;              // Number of required parameters
  int numFlagsNoParams;        // Number of flags not taking parameters
  int numFlagsParams;          // Number of flags taking parameters

  // Flag lists
  const char** flagsNoParams;  // Flags without parameters
  const char** flagsParams;    // Flags with parameters

  // Usage strings
  string mainUsage;       // The main usage string
  string paramUsage;      // The parameters usage string
  string flagUsage1;      // The flags without parameters usage string
  string flagUsage2;      // The flags with parameters usage string

  // The number of command line arguments passed in
  int argc;

  // The command line arguments themselves that were passed in
  char** argv;

  // The array of indices at which parameters can be found
  int* paramIndices;

  // The array of locations at which flags can be found
  map<const char*, int> flagIndices;
};

#endif /* PARSECOMMANDLINE_H */
