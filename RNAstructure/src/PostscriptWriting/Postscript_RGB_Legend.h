/*
 * A header file for a class that writes a Postscript color coded legend whose
 * color increments change evenly from red (lowest values) to green (median
 * values) to blue (highest values), displayed with red first.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#ifndef POSTSCRIPT_RGB_LEGEND_H
#define POSTSCRIPT_RGB_LEGEND_H

#include "Postscript_Legend.h"

class Postscript_RGB_Legend : public Postscript_Legend {
 public:
  // Public destructor and mutator

  /*
   * Name:        Destructor
   * Description: Actually has nothing to do, here as a placeholder
   */
  ~Postscript_RGB_Legend();

  /*
   * Name:        setBounds
   * Description: Set the minimum and maximum values contained in the legend.
   * Arguments:
   *     1. min
   *        The minimum value of the legend.
   *     2. max
   *        The maximum value of the legend.
   */
  void setBounds( double min, double max );

 protected:
  // Protected constructor and methods

  /*
   * Name:        Constructor
   * Description: Initializes private variables.
   * Arguments:
   *     1. entries
   *        The number of entries in the legend.
   *     2. description
   *        The description string for each segment of the legend.
   */
  Postscript_RGB_Legend( int entries, string description );

  // Method inherited from Postscript_Legend.h; see file for description.
  virtual void setColors();

  // Method inherited from Postscript_Legend.h; see file for description.
  virtual void setTexts();

 private:
  // The divider/identifier text used on each segment of the legend.
  string divider;

  // The minimum value of the legend.
  double minimum;

  // The maximum value of the legend.
  double maximum;
};

#endif /* POSTSCRIPT_RGB_LEGEND_H */
