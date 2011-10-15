/*
 * An implementation file for a class that writes a Postscript color coded
 * legend for base pair probability data.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#include "Postscript_Probability_Legend.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Probability_Legend::Postscript_Probability_Legend()
  : Postscript_Legend( 8 ) {

}

///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Probability_Legend::~Postscript_Probability_Legend() {}

///////////////////////////////////////////////////////////////////////////////
// Set the colors used in the legend.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Probability_Legend::setColors() {

  colors[0] = RED;
  colors[1] = ORANGE;
  colors[2] = DARK_YELLOW;
  colors[3] = GREEN;
  colors[4] = BRIGHT_GREEN;
  colors[5] = LIGHT_BLUE;
  colors[6] = BLUE;
  colors[7] = DARK_PINK;

}

///////////////////////////////////////////////////////////////////////////////
// Set the texts used in the legend.
///////////////////////////////////////////////////////////////////////////////
void Postscript_Probability_Legend::setTexts() {

  texts[0] = "      BP Probability >= 99%";
  texts[1] = "99% > BP Probability >= 95%";
  texts[2] = "95% > BP Probability >= 90%";
  texts[3] = "90% > BP Probability >= 80%";
  texts[4] = "80% > BP Probability >= 70%";
  texts[5] = "70% > BP Probability >= 60%";
  texts[6] = "60% > BP Probability >= 50%";
  texts[7] = "50% > BP Probability";

}

