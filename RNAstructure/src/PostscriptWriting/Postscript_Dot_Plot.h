/*
 * A header file for a class that writes a triangular Postscript dot plot.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#ifndef POSTSCRIPT_DOT_PLOT_H
#define POSTSCRIPT_DOT_PLOT_H

#include "../DotPlotHandler.h"
#include "Postscript_Image.h"
#include "Postscript_RGB_Legend.h"

class Postscript_Dot_Plot : public Postscript_Image {
 public:
  // Public methods.

  /*
   * Name:        getHandler
   * Description: Gets the dot plot handler that handles this plot.
   */
  DotPlotHandler* getHandler();

 protected:
  // Protected constructor, destructor, and methods

  /*
   * Name:        Constructor
   * Description: Initializes private variables.
   * Arguments:
   *     1. handler
   *        The dot plot handler that holds dot plot data.
   */
  Postscript_Dot_Plot( DotPlotHandler* handler );

  /*
   * Name:        Destructor
   * Description: Deletes dot plot handler.
   */
  virtual ~Postscript_Dot_Plot();

  /*
   * Name:        writeSpecificKey
   * Description: Write a legend for a dot plot.
   *              This method is a wrapper for the main writeKey method, which
   *              ensures that the proper legend type is given for a
   *              subclass.
   * Arguments:
   *     1. out
   *        The output stream to the output file currently being written.
   */
  virtual void writeSpecificKey( ofstream &out );

  /*
   * Name:        writeSpecificKey
   * Description: Write a legend for a dot plot.
   *              This method is different than the ofstream legend because it
   *              takes a pointer to an existing legend and does common tasks
   *              on it for creation and display.
   * Arguments:
   *     1. legend
   *        The legend to write.
   *     2. out
   *        The output stream to the output file currently being written.
   */
  void writeKey( Postscript_RGB_Legend* legend, ofstream &out );

  // Methods inherited from Postscript_Object.h; see that file for descriptions
  virtual void writeDescriptor( ofstream &out );
  virtual void writeSpecificImageType( ofstream &out );
  virtual void writeVariables( ofstream &out );

 protected:
  // Protected back end dot plot handler.
  DotPlotHandler* plotHandler;

 private:
  // Private methods.

  /*
   * Name:        writeDots
   * Description: Write the dots on the plot.
   * Arguments:
   *     1. out
   *        The output stream to the output file currently being written.
   */
  void writeDots( ofstream &out );

  /*
   * Name:        writeGrid
   * Description: Write the plot border, grid lines, and numbers corresponding
   *              to them.
   * Arguments:
   *     1. out
   *        The output stream to the output file currently being written.
   */
  void writeGrid( ofstream &out );
};

#endif /* POSTSCRIPT_DOT_PLOT_H */
