/*
 * A header file for a class that writes a Postscript color coded legend.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#ifndef POSTSCRIPT_LEGEND_H
#define POSTSCRIPT_LEGEND_H

#include "Postscript_Piece.h"

class Postscript_Legend : public Postscript_Piece {
 public:
  // Public constructor, destructor, and color accessor.

  /*
   * Name:        Constructor
   * Description: Initializes private variables.
   * Arguments:
   *     1. number
   *        The number of entries in the legend.
   */
  Postscript_Legend( int number );

  /*
   * Name:        Destructor
   * Description: Actually has nothing to do, here as a placeholder
   */
  virtual ~Postscript_Legend();

  /*
   * Name:        getColors
   * Description: Get the colors from a particular legend.
   */
  vector<string> getColors();

  /*
   * Name:        setLocation
   * Description: Set the location of the legend on the Postscript page.
   * Arguments:
   *     1. x
   *        The x coordinate of the top left corner baseline of the legend.
   *     2. y
   *        The y coordinate of the top left corner baseline of the legend.
   */
  void setLocation( int x, int y );

  /*
   * Name:        writeLegend
   * Description: Write the legend into a stream, which can lead to a string
   *              or to a file.
   * Arguments:
   *     1. out
   *        A reference to a stream the legend is written to.
   */
  void writeLegend( ostream &out );

 protected:
  // Protected methods

  /*
   * Name:        setColors
   * Description: Set, in order they are shown on the legend, the colors used
   *              in the legend. The colors for a particular type of legend are
   *              detemined inside this method, not fed into it.
   */
  virtual void setColors();

  /*
   * Name:        setTexts
   * Description: Set, in order they are shown on the legend, the texts used in
   *              the legend. The texts for a particular type of legend are
   *              detemined inside this method, not fed into it.
   *              Note that this method places strings as text without changing
   *              them, so whitespace is preserved.
   */
  virtual void setTexts();

 protected:
  // Protected variables

  // Vector of colors in the legend.
  vector<string> colors;

  // The x coordinate of the legend's upper left hand corner.
  int locationX;

  // The y coordinate of the legend's upper left hand corner.
  int locationY;

  // Vector of texts in the legend.
  vector<string> texts;
};

#endif /* POSTSCRIPT_LEGEND_H */
