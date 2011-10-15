/*
 * A header file for a class that writes a Postscript image.
 *
 * (c) 2010  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#ifndef POSTSCRIPT_IMAGE_H
#define POSTSCRIPT_IMAGE_H

#include "Postscript_Object.h"

class Postscript_Image : public Postscript_Object {
 public:
  // Public output writing methods

  /*
   * Name:        setAdditionalData
   * Description: Set additional data onto the structure in the form of
   *              captions, legends, titles, etc. The data given to this
   *              method must be strings of Postscript, laid out exactly as
   *              they should be in the final file.
   * Arguments:
   *     1. data
   *        The vector of additional Postscript data.
   */
  void setAdditionalData( vector<string> data );

  /*
   * Name:        writeOutput
   * Description: Write the Postscript output file.
   * Arguments:
   *     1. append
   *        True if the output is being appended to an existing file, false if
   *        beginning a new file.
   *        Default is false.
   */
  void writeOutput( bool append = false );

 protected:
  // Protected constructors, destructor, and methods

  /*
   * Name:        Constructor
   * Description: Initializes all protected variables.
   * Arguments:
   *     1. The input file name.
   *     2. The output Postscript file.
   */
  Postscript_Image( string input, string output );

  /*
   * Name:        Destructor
   * Description: Deletes all dynamically allocated variables.
   */
  virtual ~Postscript_Image();

  /*
   * Name:        initializeStrandHelpers
   * Description: Initialize important strand helpers, such as dynamically
   *              allocating the error checker and getting its sequence length.
   *              This method is separated from the constructor so subclasses
   *              can tailor a strand to their needs before setting it.
   */
  void initializeStrandHelpers();

  /*
   * Name:        writeAdditionalData
   * Description: Write additional information about the image object after the
   *              image has been drawn and the original coordinate system
   *              restored.
   */
  void writeAdditionalData( ofstream &out );

  /*
   * Name:        writeDescriptor
   * Description: Write the description string of the image.
   * Arguments:
   *     1. out
   *        The output stream to the output file currently being written.
   */
  virtual void writeDescriptor( ofstream &out );

  /*
   * Name:        writeSpecificImageType
   * Description: Write a series of Postscript commands that create a
   *              particular type of image.
   * Arguments:
   *     1. out
   *        The output stream to the output file currently being written.
   */
  virtual void writeSpecificImageType( ofstream &out );

  /*
   * Name:        writeVariables
   * Description: Write the variables that are unique to a particular image.
   * Arguments:
   *     1. out
   *        The output stream to the output file currently being written.
   */
  virtual void writeVariables( ofstream &out );

 protected:
  // Protected variables

  // Input and output file names.
  string filename;        // Input file from which data comes
  string psFile;          // Postscript file written from data

  // The RNA strand that handles this object.
  RNA* strand;

  // The error checker for the RNA strand that handles this object.
  ErrorChecker<RNA>* checker;

  // The length of the RNA strand.
  int strandLength;

  // The vector of additional data to write in the Postscript file after the
  // image is done. This array is, in order, line by line, the exact contents
  // of the additional data in the file.
  vector<string> additional;
};

#endif /* POSTSCRIPT_IMAGE_H */
