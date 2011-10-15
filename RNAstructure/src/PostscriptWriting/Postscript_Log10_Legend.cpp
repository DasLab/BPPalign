/*
 * An implementation file for a class that writes a Postscript color coded
 * legend whose color increments change evenly from red (lowest values) to
 * green (median values) to blue (highest values), displayed with red first.
 * This particular variant displays negative log10 base pair probabilities.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#include "Postscript_Log10_Legend.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
Postscript_Log10_Legend::Postscript_Log10_Legend( int entries )
  : Postscript_RGB_Legend( entries, "-log10(BP Probability)" ) {

}

///////////////////////////////////////////////////////////////////////////////
// Destructor.
///////////////////////////////////////////////////////////////////////////////
Postscript_Log10_Legend::~Postscript_Log10_Legend() {}
