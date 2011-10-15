/*
 * An implementation file for a class that writes a triangular Postscript dot
 * plot.
 * This particular variant holds free energy values from folding of a single
 * strand of nucleic acids.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#include "Postscript_Energy_Plot.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Energy_Plot::Postscript_Energy_Plot( DotPlotHandler* handler ) 
  : Postscript_Single_Plot( handler, 4 ) {

}

///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Energy_Plot::~Postscript_Energy_Plot() {}

///////////////////////////////////////////////////////////////////////////////
// Write the key for the plot.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Energy_Plot::writeSpecificKey( ofstream &out ) {

  int entries = plotHandler->getEntries();
  writeKey( new Postscript_Energy_Legend( entries ), out );

}

