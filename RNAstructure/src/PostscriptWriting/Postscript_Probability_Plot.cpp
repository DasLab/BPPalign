/*
 * An implementation file for a class that writes a triangular Postscript dot
 * plot.
 * This particular variant holds base pair probabilities.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#include "Postscript_Probability_Plot.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Probability_Plot::Postscript_Probability_Plot( DotPlotHandler* dh ) 
  : Postscript_Single_Plot( dh, 3 ) {

}

///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Probability_Plot::~Postscript_Probability_Plot() {}

///////////////////////////////////////////////////////////////////////////////
// Write the key for the plot.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Probability_Plot::writeSpecificKey( ofstream &out ) {

  int entries = plotHandler->getEntries();
  writeKey( new Postscript_Log10_Legend( entries ), out );

}
