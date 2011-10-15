/*
 * An implementation file for a class that creates Postscript files for images
 * in the RNAstructure program. This class encapsulates common functions of
 * the Postscript library in utility methods, like plot and structure drawing.
 * Note that each of these utility methods is self-contained and manages its
 * own memory, so they should be used only as miniature stand-alone utilities.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Edited for use with Postscript Library, April 2010
 * Written by Jessica Reuter
 */

#include "Postscript_Wrapper.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Wrapper::Postscript_Wrapper() {}

///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Wrapper::~Postscript_Wrapper() {}

///////////////////////////////////////////////////////////////////////////////
// Write a Postscript Dynalign dot plot, for sequence 1 of the Dynalign file.
///////////////////////////////////////////////////////////////////////////////
bool Postscript_Wrapper::plotDynalign1( DotPlotHandler* handler ) {

  // Set the plot type, just in case it hasn't been set already.
  handler->setPlotType( TYPE_DYNALIGN1 );

  // Parse data for the plot.
  handler->readDynalignSeq1Data();
  if( handler->isError() ) { return false; }

  // Create the plot object.
  Postscript_Dynalign_Plot* data = new Postscript_Dynalign_Plot( handler );
  if( data->isError() ) { return false; }

  // Write the plot object to a file.
  data->writeOutput();
  if( data->isError() ) {
    delete data;
    return false;
  }

  // Delete the plot object.
  delete data;
  return true;  

}

///////////////////////////////////////////////////////////////////////////////
// Write a Postscript Dynalign dot plot, for sequence 2 of the Dynalign file.
///////////////////////////////////////////////////////////////////////////////
bool Postscript_Wrapper::plotDynalign2( DotPlotHandler* handler ) {

  // Set the plot type, just in case it hasn't been set already.
  handler->setPlotType( TYPE_DYNALIGN2 );

  // Parse data for the plot.
  handler->readDynalignSeq2Data();

  // Create the plot object.
  Postscript_Dynalign_Plot* data = new Postscript_Dynalign_Plot( handler );
  if( data->isError() ) { return false; }

  // Write the plot object to a file.
  data->writeOutput();
  if( data->isError() ) {
    delete data;
    return false;
  }

  // Delete the plot object.
  delete data;
  return true;

}

///////////////////////////////////////////////////////////////////////////////
// Write a Postscript free energy dot plot.
///////////////////////////////////////////////////////////////////////////////
bool Postscript_Wrapper::plotEnergy( DotPlotHandler* handler ) {

  // Set the plot type, just in case it hasn't been set already.
  handler->setPlotType( TYPE_ENERGY );

  // Parse data for the plot.
  handler->readFoldingData();

  // Create the plot object.
  Postscript_Energy_Plot* data = new Postscript_Energy_Plot( handler );
  if( data->isError() ) { return false; }

  // Write the plot object to a file.
  data->writeOutput();
  if( data->isError() ) {
    delete data;
    return false;
  }

  // Delete the plot object.
  delete data;
  return true;

}

///////////////////////////////////////////////////////////////////////////////
// Write a Postscript base pair probabilities dot plot.
///////////////////////////////////////////////////////////////////////////////
bool Postscript_Wrapper::plotProbability( DotPlotHandler* handler ) {

  // Set the plot type, just in case it hasn't been set already.
  handler->setPlotType( TYPE_PROBABILITY );

  // Parse data for the plot.
  handler->readPartitionData();

  // Create the plot object.
  Postscript_Probability_Plot* data =
    new Postscript_Probability_Plot( handler );
  if( data->isError() ) { return false; }

  // Write the plot object to a file.
  data->writeOutput();
  if( data->isError() ) {
    delete data;
    return false;
  }

  // Delete the plot object.
  delete data;
  return true;

}

///////////////////////////////////////////////////////////////////////////////
// Write Postscript circular structures without annotation.
///////////////////////////////////////////////////////////////////////////////
bool Postscript_Wrapper::structureCircular( string in, string out ) {

  // Create the structure(s).
  Postscript_Circle_Structure* data =
    new Postscript_Circle_Structure( in, out );
  if( data->isError() ) { return false; }

  // Write the structure(s) to a file.
  data->writeAllStructures();
  if( data->isError() ) {
    delete data;
    return false;
  }

  // Delete the structure(s).
  delete data;
  return true;

}

///////////////////////////////////////////////////////////////////////////////
// Write Postscript circular structures with base pair probability annotation.
///////////////////////////////////////////////////////////////////////////////
bool Postscript_Wrapper::structureCircular_Probability( string in, string out,
							string file ) {

  // Create the structure(s).
  Postscript_Circle_Structure* data =
    new Postscript_Circle_Structure( in, out );
  if( data->isError() ) { return false; }

  // Add probability annotation.
  data->addProbabilityAnnotation( file );
  if( data->isError() ) {
    delete data;
    return false;
  }

  // Write the structure(s) to a file.
  data->writeAllStructures();
  if( data->isError() ) {
    delete data;
    return false;
  }

  // Delete the structure(s).
  delete data;
  return true;

}

///////////////////////////////////////////////////////////////////////////////
// Write Postscript circular structures with SHAPE annotation.
///////////////////////////////////////////////////////////////////////////////
bool Postscript_Wrapper::structureCircular_SHAPE( string in, string out,
						  string file ) {

  // Create the structure(s).
  Postscript_Circle_Structure* data =
    new Postscript_Circle_Structure( in, out );
  if( data->isError() ) { return false; }

  // Add SHAPE annotation.
  data->addSHAPEAnnotation( file );
  if( data->isError() ) {
    delete data;
    return false;
  }

  // Write the structure(s) to a file.
  data->writeAllStructures();
  if( data->isError() ) {
    delete data;
    return false;
  }

  // Delete the structure(s).
  delete data;
  return true;

}

///////////////////////////////////////////////////////////////////////////////
// Write Postscript radial structures without annotation.
///////////////////////////////////////////////////////////////////////////////
bool Postscript_Wrapper::structureRadial( string in, string out ) {

  // Create the structure(s).
  Postscript_Radial_Structure* data =
    new Postscript_Radial_Structure( in, out );
  if( data->isError() ) { return false; }

  // Write the structure(s) to a file.
  data->writeAllStructures();
  if( data->isError() ) {
    delete data;
    return false;
  }

  // Delete the structure(s).
  delete data;
  return true;

}

///////////////////////////////////////////////////////////////////////////////
// Write Postscript radial structures with base pair probability annotation.
///////////////////////////////////////////////////////////////////////////////
bool Postscript_Wrapper::structureRadial_Probability( string in, string out,
						      string file ) {

  // Create the structure(s).
  Postscript_Radial_Structure* data =
    new Postscript_Radial_Structure( in, out );
  if( data->isError() ) { return false; }

  // Add probability annotation.
  data->addProbabilityAnnotation( file );
  if( data->isError() ) {
    delete data;
    return false;
  }

  // Write the structure(s) to a file.
  data->writeAllStructures();
  if( data->isError() ) {
    delete data;
    return false;
  }

  // Delete the structure(s).
  delete data;
  return true;

}

///////////////////////////////////////////////////////////////////////////////
// Write Postscript radial structures with SHAPE annotation.
///////////////////////////////////////////////////////////////////////////////
bool Postscript_Wrapper::structureRadial_SHAPE( string in, string out,
						string file ) {
  // Create the structure(s).
  Postscript_Radial_Structure* data =
    new Postscript_Radial_Structure( in, out );
  if( data->isError() ) { return false; }

  // Add SHAPE annotation.
  data->addSHAPEAnnotation( file );
  if( data->isError() ) {
    delete data;
    return false;
  }

  // Write the structure(s) to a file.
  data->writeAllStructures();
  if( data->isError() ) {
    delete data;
    return false;
  }

  // Delete the structure(s).
  delete data;
  return true;

}
