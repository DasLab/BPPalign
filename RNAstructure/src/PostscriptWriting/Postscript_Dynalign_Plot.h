/*
 * A header file for a class that writes a triangular Postscript dot plot.
 * This particular variant holds free energy values from a Dynalign
 * calculation, specified as either strand 1 or strand 2.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#ifndef POSTSCRIPT_DYNALIGN_PLOT_H
#define POSTSCRIPT_DYNALIGN_PLOT_H

#include "Postscript_Dot_Plot.h"
#include "Postscript_Energy_Legend.h"

class Postscript_Dynalign_Plot : public Postscript_Dot_Plot {
 public:
  // Public constructor and destructor

  /*
   * Name:        Constructor
   * Description: Wrapper to call the superclass and create the Dynalign object
   *              necessary for this plot.
   * Arguments:
   *     1. handler
   *        The dot plot handler that holds plot data.
   */
  Postscript_Dynalign_Plot( DotPlotHandler* handler );

  /*
   * Name:        Destructor
   * Description: Deletes the Dynalign object used for this plot.
   */
  ~Postscript_Dynalign_Plot();

  // Method inherited from Postscript_Dot_Plot.h; see that file for description
  virtual void writeSpecificKey( ofstream &out );
};

#endif /* POSTSCRIPT_DYNALIGN_PLOT_H */
