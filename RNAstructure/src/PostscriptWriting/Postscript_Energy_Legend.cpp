/*
 * An implementation file for a class that writes a Postscript color coded
 * legend whose color increments change evenly from red (lowest values) to
 * green (median values) to blue (highest values), displayed with red first.
 * This particular variant displays folding free energies.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#include "Postscript_Energy_Legend.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
Postscript_Energy_Legend::Postscript_Energy_Legend( int entries )
  : Postscript_RGB_Legend( entries, "Free Energy (kcal/mol)" ) {

}

///////////////////////////////////////////////////////////////////////////////
// Destructor.
///////////////////////////////////////////////////////////////////////////////
Postscript_Energy_Legend::~Postscript_Energy_Legend() {}
