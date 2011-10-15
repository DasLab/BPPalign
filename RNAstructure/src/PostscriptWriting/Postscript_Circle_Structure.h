/*
 * A header file for a class that writes a circular Postscript structure from
 * varied input sources.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#ifndef POSTSCRIPT_CIRCLE_STRUCTURE_H
#define POSTSCRIPT_CIRCLE_STRUCTURE_H

#include "Postscript_Structure.h"

class Postscript_Circle_Structure : public Postscript_Structure {
 public:
  // Public constructor and destructor

  /*
   * Name:        Constructor
   * Description: Constructor which initializes the RNA data structure for this
   *              class directly, creating the RNA object inside itself.
   * Arguments:
   *     1. in
   *        The name of the file to create an internal RNA object with.
   *     2. out
   *        The name of the Postscript output file to write.
   */
  Postscript_Circle_Structure( string in, string out );

  /*
   * Name:        Destructor
   * Description: Actually has nothing to destruct, here as a placeholder
   */
  ~Postscript_Circle_Structure();

 protected:
  // Methods inherited from Postscript_Object.h; see file for descriptions
  virtual void writeDescriptor( ofstream &out );
  virtual void writeVariables( ofstream &out );

  // Methods inherited from Postscript_Structure.h; see file for descriptions
  virtual void placeLegend( Postscript_Legend* legend );
  virtual void writeBackbone( ofstream &out );
  virtual void writePair( ofstream &out );
};

#endif /* POSTSCRIPT_CIRCLE_STRUCTURE_H */
