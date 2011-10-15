/*
 * A header file for a class that creates Postscript files for images in the
 * RNAstructure program. This class encapsulates common functions of the
 * Postscript library in utility methods, like plot and structure drawing.
 * Note that each of these utility methods is self-contained and manages its
 * own memory, so they should be used only as miniature stand-alone utilities.
 *
 * Note that this class is the Postscript Library's entry point with SWIG, when
 * used with the RNAstructure java interface. Compliance with SWIG is why the
 * constructor and destructor are explicitly defined, even though they don't
 * need to be from a C++ perspective because they're empty.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Edited for use with Postscript Library, April 2010
 * Written by Jessica Reuter
 */

#ifndef POSTSCRIPT_WRAPPER_H
#define POSTSCRIPT_WRAPPER_H

#include "Postscript_Circle_Structure.h"
#include "Postscript_Dynalign_Plot.h"
#include "Postscript_Energy_Plot.h"
#include "Postscript_Probability_Plot.h"
#include "Postscript_Radial_Structure.h"

class Postscript_Wrapper {
 public:
  // Public constuctor, destructor, and utility methods

  /*
   * Name:        Constructor
   * Description: Does nothing since utility methods are self-contained, here
   *              for class completeness.
   */
  Postscript_Wrapper();

  /*
   * Name:        Destructor
   * Description: Does nothing since utility methods are self-contained, here
   *              for class completeness.
   */
  ~Postscript_Wrapper();

  /*
   * Name:        plotDynalign1
   * Description: Writes a Postscript file of a Dynalign dot plot.
   *              This plot handles sequence 1 from the Dynalign calculation.
   * Arguments:
   *     1. in
   *        The input file to make the Postscript image from.
   *     2. out
   *        The output file the Postscript image is written to.
   *     3. entries
   *        The number of entries in the plot's legend.
   *     4. minimum
   *        The minimum bound of the plot.
   *     5. maximum
   *        The maximum bound of the plot.
   * Returns:
   *     True if no error occurred during image file creation, false if an
   *     error did occur.
   */
  bool plotDynalign1( DotPlotHandler* handler );

  /*
   * Name:        plotDynalign2
   * Description: Writes a Postscript file of a Dynalign dot plot.
   *              This plot handles sequence 2 from the Dynalign calculation.
   * Arguments:
   *     1. in
   *        The input file to make the Postscript image from.
   *     2. out
   *        The output file the Postscript image is written to.
   *     3. entries
   *        The number of entries in the plot's legend.
   *     4. minimum
   *        The minimum bound of the plot.
   *     5. maximum
   *        The maximum bound of the plot.
   * Returns:
   *     True if no error occurred during image file creation, false if an
   *     error did occur.
   */
  bool plotDynalign2( DotPlotHandler* handler );

  /*
   * Name:        plotEnergy
   * Description: Writes a Postscript file of a folding free energy dot plot.
   * Arguments:
   *     1. in
   *        The input file to make the Postscript image from.
   *     2. out
   *        The output file the Postscript image is written to.
   *     3. entries
   *        The number of entries in the plot's legend.
   *     4. minimum
   *        The minimum bound of the plot.
   *     5. maximum
   *        The maximum bound of the plot.
   * Returns:
   *     True if no error occurred during image file creation, false if an
   *     error did occur.
   */
  bool plotEnergy( DotPlotHandler* handler );

  /*
   * Name:        plotProbability
   * Description: Writes a Postscript file of a base pair probability dot plot.
   * Arguments:
   *     1. in
   *        The input file to make the Postscript image from.
   *     2. out
   *        The output file the Postscript image is written to.
   *     3. entries
   *        The number of entries in the plot's legend.
   *        Default is 5 entries.
   *     4. minimum
   *        The minimum bound of the plot.
   *     5. maximum
   *        The maximum bound of the plot.
   * Returns:
   *     True if no error occurred during image file creation, false if an
   *     error did occur.
   */
  bool plotProbability( DotPlotHandler* handler );

  /*
   * Name:        structureCircular
   * Description: Creates a Postscript file of all structures in a CT file,
   *              using a circular layout for each structure, with no
   *              annotation on bases or base pairs.
   * Arguments:
   *     1. in
   *        The input file to make the Postscript image from.
   *     2. out
   *        The output file the Postscript image is written to.
   * Returns:
   *     True if no error occurred during image file creation, false if an
   *     error did occur.
   */
  bool structureCircular( string in, string out );

  /*
   * Name:        structureCircular_Probability
   * Description: Creates a Postscript file of all structures in a CT file,
   *              using a circular layout for each structure, with base pairing
   *              probability data annotating base pairs.
   * Arguments:
   *     1. in
   *        The input file to make the Postscript image from.
   *     2. out
   *        The output file the Postscript image is written to.
   *     3. file
   *        The data file which contains base pair probabilities.
   * Returns:
   *     True if no error occurred during image file creation, false if an
   *     error did occur.
   */
  bool structureCircular_Probability( string in, string out, string file );

  /*
   * Name:        structureCircular_SHAPE
   * Description: Creates a Postscript file of all structures in a CT file,
   *              using a circular layout for each structure, with SHAPE data
   *              annotating bases.
   * Arguments:
   *     1. in
   *        The input file to make the Postscript image from.
   *     2. out
   *        The output file the Postscript image is written to.
   *     3. file
   *        The SHAPE data file used for annotation.
   * Returns:
   *     True if no error occurred during image file creation, false if an
   *     error did occur.
   */
  bool structureCircular_SHAPE( string in, string out, string file );

  /*
   * Name:        structureRadial
   * Description: Creates a Postscript file of all structures in a CT file,
   *              using a radial layout for each structure, with no annotation
   *              on bases or base pairs.
   * Arguments:
   *     1. in
   *        The input file to make the Postscript image from.
   *     2. out
   *        The output file the Postscript image is written to.
   * Returns:
   *     True if no error occurred during image file creation, false if an
   *     error did occur.
   */
  bool structureRadial( string in, string out );

  /*
   * Name:        structureRadial_Probability
   * Description: Creates a Postscript file of all structures in a CT file,
   *              using a radial layout for each structure, with base pairing
   *              probability data annotating base pairs.
   * Arguments:
   *     1. in
   *        The input file to make the Postscript image from.
   *     2. out
   *        The output file the Postscript image is written to.
   *     3. file
   *        The data file which contains base pair probabilities.
   * Returns:
   *     True if no error occurred during image file creation, false if an
   *     error did occur.
   */
  bool structureRadial_Probability( string in, string out, string file );

  /*
   * Name:        structureRadial_SHAPE
   * Description: Creates a Postscript file of all structures in a CT file,
   *              using a radial layout for each structure, with SHAPE data
   *              annotating bases.
   * Arguments:
   *     1. in
   *        The input file to make the Postscript image from.
   *     2. out
   *        The output file the Postscript image is written to.
   *     3. file
   *        The SHAPE data file used for annotation.
   * Returns:
   *     True if no error occurred during image file creation, false if an
   *     error did occur.
   */
  bool structureRadial_SHAPE( string in, string out, string file );
};

#endif /* POSTSCRIPT_WRAPPER_H */
