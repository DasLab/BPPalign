/*
 * A header file for a class that writes a triangular Postscript dot plot.
 * This particular variant holds base pair probabilities.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#ifndef POSTSCRIPT_PROBABILITY_PLOT_H
#define POSTSCRIPT_PROBABILITY_PLOT_H

#include "Postscript_Single_Plot.h"
#include "Postscript_Log10_Legend.h"

class Postscript_Probability_Plot : public Postscript_Single_Plot {
 public:
  // Public constructor and destructor

  /*
   * Name:        Constructor
   * Description: Wrapper to call the superclass properly.
   * Arguments:
   *     1. dh
   *        The dot plot handler that holds plot data.
   */
  Postscript_Probability_Plot( DotPlotHandler* dh );

  /*
   * Name:        Destructor
   * Description: Actually does not delete anything, just here as a placeholder
   */
  ~Postscript_Probability_Plot();

  // Method inherited from Postscript_Dot_Plot.h; see that file for description
  virtual void writeSpecificKey( ofstream &out );
};

#endif /* POSTSCRIPT_PROBABILITY_PLOT_H */
