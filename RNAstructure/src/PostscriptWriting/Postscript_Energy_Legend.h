/*
 * A header file for a class that writes a Postscript color coded legend whose
 * color increments change evenly from red (lowest values) to green (median
 * values) to blue (highest values), displayed with red first.
 * This particular variant displays folding free energies.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#ifndef POSTSCRIPT_ENERGY_LEGEND_H
#define POSTSCRIPT_ENERGY_LEGEND_H

#include "Postscript_RGB_Legend.h"

class Postscript_Energy_Legend : public Postscript_RGB_Legend {
 public:
  // Public constructor and destructor

  /*
   * Name:        Constructor
   * Description: Initializes unique legend description and variables in the
   *              superclass.
   * Arguments:
   *     1. entries
   *        The number of entries in the legend.
   */
  Postscript_Energy_Legend( int entries );

  /*
   * Name:        Destructor
   * Description: Actually has nothing to do, here as a placeholder
   */
  ~Postscript_Energy_Legend();
};

#endif /* POSTSCRIPT_ENERGY_LEGEND_H */
