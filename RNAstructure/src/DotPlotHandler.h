/*
 * A header file for a class that holds a dot plot data structure and handles
 * all dot plot manipulation. This class can be used throughout RNAstructure,
 * for any type of dot plot.
 *
 * Copyright 2010, Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#ifndef DOT_PLOT_HANDLER_H
#define DOT_PLOT_HANDLER_H

#include <limits>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "ErrorChecker.h"
#include "../RNA_class/Dynalign_object.h"
#include "../RNA_class/RNA.h"

// Small namespace to handle dot plot constants.
namespace DotPlot_Constants {

  // Default, minimum, and maximum number of entries/colors.
  const int DEFAULT_ENTRIES = 5;
  const int MIN_COLORS = 3;
  const int MAX_COLORS = 15;

  // Constants that define the plot type.
  const char TYPE_DYNALIGN1 = '1';
  const char TYPE_DYNALIGN2 = '2';
  const char TYPE_ENERGY = 'e';
  const char TYPE_PROBABILITY = 'p';
  const char TYPE_UNDEFINED = 'X';
};

// Namespace usage declarations.
using namespace std;
using namespace DotPlot_Constants;

// Dot plot handler class.
class DotPlotHandler {
 public:

  /*
   * Name:        Constructor
   * Description: Initializes private variables.
   * Arguments:
   *     1. inFile
   *        The input file from which dot plot data is read.
   *     2. outFile
   *        The output file to which dot plot data is written.
   *        The handler by itself can only write text files, but other file
   *        types can be written when the handler is used in conjunction with
   *        other classes.
   */
  DotPlotHandler( string inFile, string outFile );

  /*
   * Name:        Destructor
   * Description: Doesn't do anything, here as a placeholder for completeness.
   */
  ~DotPlotHandler();

  /*
   * Name:        getColorAsPostscript
   * Description: Get the color of a particular dot, formatted as a Postscript
   *              value string, where each value in the string is a Postscript
   *              color value specified to two decimal places.
   *              For example: "1.00 0.50 0.00"
   * Arguments:
   *     1. i
   *        The paired nucleotide, one-indexed, in the first sequence.
   *     2. j
   *        The paired nucleotide, one-indexed, in the second sequence.
   * Returns:
   *     The color of the dot, formatted as a Postscript value string.
   */
  string getColorAsPostscript( const int i, const int j );

  /*
   * Name:        getColorAsRGB
   * Description: Get the color of a particular dot, formatted as an RGB value
   *              string, where each value in the string is an integer.
   *              For example: "255 0 0"
   * Arguments:
   *     1. i
   *        The paired nucleotide, one-indexed, in the first sequence.
   *     2. j
   *        The paired nucleotide, one-indexed, in the second sequence.
   * Returns:
   *     The color of the dot, formatted as an RGB value string.
   */
  string getColorAsRGB( const int i, const int j );

  /*
   * Name:        getEntries
   * Description: Get the number of entries/thresholds in the plot, i.e. the
   *              number of divisions in the legend for this plot.
   * Returns:
   *     The number of entries in the plot's legend.
   */
  int getEntries();

  /*
   * Name:        getInfoString
   * Description: Get the information string that describes the contents of the
   *              dot plot.
   * Returns:
   *     The dot plot information string.
   */
  string getInfoString();

  /*
   * Name:        getInputFile
   * Description: Get the name of the input file for the dot plot, as was given
   *              to the handler constructor.
   * Returns:
   *     The dot plot input file name.
   */
  string getInputFile();

  /*
   * Name:        getMaximum
   * Description: Get the maximum value displayed on the dot plot. Note that
   *              this may not be the maximum allowable value in the dot plot;
   *              if this method is called after using setMaximum (see below),
   *              it will return the maximum as it was set. Conversely, if
   *              this method is called without using setMaximum or after using
   *              resetRange (see below), it will return the default maximum.
   * Returns:
   *     The maximum viewable value in the dot plot.
   */
  double getMaximum();

  /*
   * Name:        getMinimum
   * Description: Get the minimum value displayed on the dot plot. Note that
   *              this may not be the minimum allowable value in the dot plot;
   *              if this method is called after using setMinimum (see below),
   *              it will return the minimum as it was set. Conversely, if
   *              this method is called without using setMinimum or after using
   *              resetRange (see below), it will return the default minimum.
   * Returns:
   *     The minimum viewable value in the dot plot.
   */
  double getMinimum();

  /*
   * Name:        getOutputFile
   * Description: Get the name of the output file for the dot plot, as was
   *              given to the handler constructor.
   * Returns:
   *     The dot plot output file name.
   */
  string getOutputFile();

  /*
   * Name:        getPlotType
   * Description: Get the character code that identifies a particular plot
   *              type. These character codes are dentified as type constants
   *              (see below).
   * Returns:
   *     The character code that identifies a plot type.
   */
  char getPlotType();

  /*
   * Name:        getValue
   * Description: Get the value of a particular dot.
   * Arguments:
   *     1. i
   *        The paired nucleotide, one-indexed, in the first sequence.
   *     2. j
   *        The paired nucleotide, one-indexed, in the second sequence.
   * Returns:
   *     The value of the dot.
   */
  double getValue( const int i, const int j );

  /*
   * Name:        isError
   * Description: Get the error status of a particular dot plot.
   * Returns:
   *     True if the dot plot operations encountered an error, false if not.
   */
  bool isError();

  /*
   * Name:        readDynalignSeq1Data
   * Description: Read plot data from a Dynalign save file, for sequence 1 of
   *              that particular calculation.
   */
  void readDynalignSeq1Data();

  /*
   * Name:        readDynalignSeq2Data
   * Description: Read plot data from a Dynalign save file, for sequence 2 of
   *              that particular calculation.
   */
  void readDynalignSeq2Data();

  /*
   * Name:        readFoldingData
   * Description: Read plot data from a folding save file.
   */
  void readFoldingData();

  /*
   * Name:        readPartitionData
   * Description: Read plot data from a partition function save file.
   */
  void readPartitionData();

  /*
   * Name:        readTextData
   * Description: Read plot data from a dot plot text file.
   */
  void readTextData();

  /*
   * Name:        resetRange
   * Description: Reset the visible range of the dot plot to the default.
   */
  void resetRange();

  /*
   * Name:        setEntries
   * Description: Set the number of entries in a dot plot.
   */
  void setEntries( int newEntries );

  /*
   * Name:        setMaximum
   * Description: Set the maximum value in a dot plot. If the value is not
   *              valid or is out of range, it is ignored.
   * Arguments:
   *     1. newMaximum
   *        The new maximum value to set.
   */
  void setMaximum( const double newMaximum );

  /*
   * Name:        setMinimum
   * Description: Set the minimum value in a dot plot. If the value is not
   *              valid or is out of range, it is ignored.
   * Arguments:
   *     1. newMinimum
   *        The new minimum value to set.
   */
  void setMinimum( const double newMinimum );

  /*
   * Name:        setPlotType
   * Description: Set the plot type. If the plot type isn't valid, the type is
   *              set to TYPE_UNDEFINED.
   * Arguments:
   *     1. type
   *        The plot type to set.
   */
  void setPlotType( char type );

  /*
   * Name:        writePlotFile
   * Description: Write dot plot data as a text file.
   */
  void writePlotFile();

  /*
 public:
  // Public constants that handle plot type and entry number.

  static const int DEFAULT_ENTRIES = 5;
  static const int MIN_COLORS = 3;
  static const int MAX_COLORS = 15;

  static const char TYPE_DYNALIGN1 = '1';
  static const char TYPE_DYNALIGN2 = '2';
  static const char TYPE_ENERGY = 'e';
  static const char TYPE_PROBABILITY = 'p';
  static const char TYPE_UNDEFINED = 'X';
  */

 private:
  // Private methods.

  /*
   * Name:        buildThresholds
   * Description: Build the threshold vectors for colors and raw values.
   */
  void buildThresholds();

  /*
   * Name:        checkBounds
   * Description: Check the bounds of the plot to make sure they're valid, and
   *              correct them if necessary to the default maximum or minimum.
   */
  void checkBounds();

  /*
   * Name:        getThresholdIndex
   * Description: Get the index (zero-indexed) in the threshold array that a 
   *              value falls under.
   * Arguments:
   *     1. i
   *        The paired nucleotide, one-indexed, in the first sequence.
   *     2. j
   *        The paired nucleotide, one-indexed, in the second sequence.
   * Returns:
   *     The threshold index.
   */
  int getThresholdIndex( const int i, const int j );

  /*
   * Name:        initializeValuesArray
   * Description: Initialize the 2D array (actually a vector) that holds the
   *              dot plot data.
   * Arguments:
   *     1. length
   *        The length of the sequence from which the plot is created.
   */
  void initializeValuesArray( int length );

  /*
   * Name:        readDynalignData
   * Description: Read dot plot data from a Dynalign calculation, focusing on
   *              a particular sequence to do so. Note that the sequence must
   *              be either 1 or 2.
   * Arguments:
   *     1. sequence
   *        The sequence from the Dynalign calculation to plot.
   */
  void readDynalignData( int sequence );

  /*
   * Name:        readSingleStrandData
   * Description: Read dot plot data from a save file holding information from
   *              a single strand.
   * Arguments:
   *     1. type
   *        The type of save file to read. This type echoes the types used in
   *        the RNA class constructor (see RNA.cpp).
   */
  void readSingleStrandData( int type );

 private:
  // Private variables.

  // Input and output file names.
  string input;
  string output;

  // Information labels.
  string label;               // Used as label in plot legend.
  string infoString;          // The description string for the entire plot.

  // Default maximum and minimum.
  double defaultMax;
  double defaultMin;

  // Currently visible maximum and minimum.
  double maximum;
  double minimum;

  // Boolean flags.
  bool dataRead;           // True if data has been read, false if not.
  bool error;              // True if an error occurred, false if not.

  // The number of entries in the plot legend.
  int entries;

  // Character code determining what type of plot is being made.
  char plotType;

  // Vector data structures.
  vector< vector<double> > colorInfo;   // Color thresholds vector.
  vector<double> thresholds;            // Value thresholds vector.
  vector< vector<double> > values;      // Vector grid of dot plot values.
};

#endif /* DOT_PLOT_HANDLER_H */
