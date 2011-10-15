/*
 * An implementation file for a class that writes a triangular Postscript dot
 * plot.
 * This particular variant holds free energy values from a Dynalign
 * calculation, specified as either strand 1 or strand 2.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#include "Postscript_Dynalign_Plot.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Dynalign_Plot::Postscript_Dynalign_Plot( DotPlotHandler* handler )
  : Postscript_Dot_Plot( handler ) {

}

///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Dynalign_Plot::~Postscript_Dynalign_Plot() {}

///////////////////////////////////////////////////////////////////////////////
// Write the key for the plot.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Dynalign_Plot::writeSpecificKey( ofstream &out ) {

  int entries = plotHandler->getEntries();
  writeKey( new Postscript_Energy_Legend( entries ), out );

}
