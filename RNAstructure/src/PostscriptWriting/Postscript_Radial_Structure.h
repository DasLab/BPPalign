/*
 * A header file for a class that writes a radial Postscript structure. The
 * "radial" in the class's name refers to the underlying algorithm; a structure
 * is not displayed as a circle.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#ifndef POSTSCRIPT_RADIAL_STRUCTURE_H
#define POSTSCRIPT_RADIAL_STRUCTURE_H

#include "Postscript_Structure.h"

class Postscript_Radial_Structure : public Postscript_Structure {
 public:
  // Public constructor and destructor

  /*
   * Name:        Constructor
   * Description: Constructor which initializes the RNA data structure for this
   *              class directly, creating the RNA object inside itself.
   * Arguments:
   *     1. file
   *        The name of the file to create an internal RNA object with.
   *     2. out
   *        The name of the Postscript output file to write.
   */
  Postscript_Radial_Structure( string file, string out );

  /*
   * Name:       Destructor
   * Destructor: Actually has nothing to destruct, here as a placeholder
   */
  ~Postscript_Radial_Structure();

 protected:
  // Protected methods

  // Methods inherited from Postscript_Object.h; see file for descriptions
  virtual void writeDescriptor( ofstream &out );
  virtual void writeVariables( ofstream &out );

  // Methods inherited from Postscript_Structure.h; see file for descriptions
  virtual void placeLegend( Postscript_Legend* legend );
  virtual void writeBackbone( ofstream &out );
  virtual void writeNucleotidesAndLabels( ofstream &out );
  virtual void writePair( ofstream &out );
};

#endif /* POSTSCRIPT_RADIAL_STRUCTURE_H */
