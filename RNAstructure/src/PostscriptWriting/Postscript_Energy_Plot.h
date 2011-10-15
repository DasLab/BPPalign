/*
 * A header file for a class that writes a triangular Postscript dot plot.
 * This particular variant holds free energy values from folding of a single
 * strand of nucleic acids.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#ifndef POSTSCRIPT_ENERGY_PLOT_H
#define POSTSCRIPT_ENERGY_PLOT_H

#include "Postscript_Single_Plot.h"
#include "Postscript_Energy_Legend.h"

class Postscript_Energy_Plot : public Postscript_Single_Plot {
 public:
  // Public constructor and destructor

  /*
   * Name:        Constructor
   * Description: Wrapper to call the superclass properly.
   * Arguments:
   *     1. in
   *        The name of the file to initialize data with.
   *     2. out
   *        The name of the output file to write a dot plot to.
   */
  Postscript_Energy_Plot( DotPlotHandler* handler );

  /*
   * Name:        Destructor
   * Description: Actually does not delete anything, just here as a placeholder
   */
  ~Postscript_Energy_Plot();

  // Method inherited from Postscript_Dot_Plot.h; see that file for description
  virtual void writeSpecificKey( ofstream &out );
};

#endif /* POSTSCRIPT_ENERGY_PLOT_H */
