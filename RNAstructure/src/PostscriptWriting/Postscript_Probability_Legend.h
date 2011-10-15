/*
 * A header file for a class that writes a Postscript color coded legend for
 * base pairing probability data.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#ifndef POSTSCRIPT_PROBABILITY_LEGEND_H
#define POSTSCRIPT_PROBABILITY_LEGEND_H

#include "Postscript_Legend.h"

class Postscript_Probability_Legend : public Postscript_Legend {
 public:
  // Public constructor and destructor

  /*
   * Name:        Constructor
   * Description: Initializes private variables.
   */
  Postscript_Probability_Legend();

  /*
   * Name:        Destructor
   * Description: Actually has nothing to do, here as a placeholder
   */
  ~Postscript_Probability_Legend();

 protected:
  // Methods inherited from Postscript_Legend.h; see file for descriptions
  virtual void setColors();
  virtual void setTexts();
};

#endif /* POSTSCRIPT_PROBABILITY_LEGEND_H */
