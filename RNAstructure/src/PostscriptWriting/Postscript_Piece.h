/*
 * A header file for a class that either writes a self-contained part of a
 * Postscript image or handles preparations for such a part. This class doesn't
 * do anything actually; it is just here for organization purposes right now.
 *
 * (c) 2010  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#ifndef POSTSCRIPT_PIECE_H
#define POSTSCRIPT_PIECE_H

#include "Postscript_Object.h"

class Postscript_Piece : public Postscript_Object {
 protected:
  // Protected constructor and destructor

  /*
   * Name:         Constructor
   * Description:  Actually doesn't do anything, just here as a placeholder.
   */
  Postscript_Piece();

  /*
   * Name:        Destructor
   * Description: Actually doesn't do anything, just here as a placeholder.
   */
  ~Postscript_Piece();
};

#endif /* POSTSCRIPT_PIECE_H */
