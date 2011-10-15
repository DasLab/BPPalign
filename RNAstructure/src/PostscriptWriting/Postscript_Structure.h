/*
 * A header file for a class that writes a Postscript RNA secondary structure.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#ifndef POSTSCRIPT_STRUCTURE_H
#define POSTSCRIPT_STRUCTURE_H

#include "Postscript_Image.h"
#include "Postscript_Annotation_Handler.h"
#include "Postscript_Probability_Legend.h"
#include "Postscript_SHAPE_Legend.h"

class Postscript_Structure : public Postscript_Image {
 public:
  // Public methods

  /*
   * Name:        addProbabilityAnnotation
   * Description: Add base pair probability annotation to the structure.
   * Arguments:
   *     1. file
   *        The name of the file which holds the necessary annotation data.
   */
  void addProbabilityAnnotation( string file );

  /*
   * Name:        addSHAPEAnnotation
   * Description: Add SHAPE annotation to the structure.
   * Arguments:
   *     1. file
   *        The name of the file which holds the necessary annotation data.
   */
  void addSHAPEAnnotation( string file );

  /*
   * Name:        calculatePseudoknottedPairs
   * Description: Calculate the number of pseudoknotted pairs in a structure.
   */
  void calculatePseudoknottedPairs();

  /*
   * Name:        getNumPairs
   * Description: Get the number of pairs in the structure.
   * Returns:
   *     The number of pairs in the structure.
   */
  virtual int getNumPairs();

  /*
   * Name:        getNumPseudoknottedPairs
   * Description: Get the number of pseudoknotted pairs in the structure.
   * Returns:
   *     The number of pseudoknotted pairs in the structure.
   */
  int getNumPseudoknottedPairs();

  /*
   * Name:        getPair
   * Description: Get a particular pairing in the structure.
   * Arguments:
   *     1. pair
   *        The nucleotide, one-indexed, whose pair to get.
   * Returns:
   *        The nucleotide, one-indexed, the nucleotide is paired to.
   */
  int getPair( int pair );

  /*
   * Name:        getStrand
   * Description: Get the RNA strand that backs this object.
   * Returns:
   *      The RNA strand.
   */
  RNA* getStrand();

  /*
   * Name:        parsePairs
   * Description: Saves the pairs in a particular structure into an easily
   *              accessible array for future use.
   * Arguments:
   *     1. structure
   *        The structure from which to parse pairs.
   *        Default is 1 (the first structure).
   */
  void parsePairs( int structure = 1 );

  /*
   * Name:        setNoAnnotation
   * Description: Set the type color code so no color annotation is allowed.
   */
  void setNoAnnotation();

  /*
   * Name:        setNucleotidesAnnotated
   * Description: Set the type color code so nucleotides can be color annotated
   *              if desired.
   */
  void setNucleotidesAnnotated();

  /*
   * Name:        setPairingsAnnotated
   * Description: Set the type color code so nucleotide pairs can be color
   *              annotated, if desired.
   */
  void setPairingsAnnotated();

  /*
   * Name:        setPairingsAndNucleotidesAnnotated
   * Description: Set the type color code so both nucleotides and nucleotide
   *              pairs can be color annotated, if desired.
   */
  void setPairingsAndNucleotidesAnnotated();

  /*
   * Name:        writeAllStructures
   * Description: Write all the structures in a particular CT file, not just
   *              one, as parsePairs does.
   */
  void writeAllStructures();

 protected:
  // Protected constructor, destructor, and methods

  /*
   * Name:        Constructor
   * Description: Constructor which initializes the RNA data structure for this
   *              class directly, creating the RNA object inside itself.
   *              Also, initializes all private variables.
   * Arguments:
   *     1. in
   *        The name of the file to create an internal RNA object with.
   *     2. out
   *        The name of the output file to write a structure to.
   */
  Postscript_Structure( string in, string out );

  /*
   * Name:        Destructor
   * Description: Deletes all dynamically allocated variables.
   */
  ~Postscript_Structure();

  /*
   * Name:        placeLegend
   * Description: Set the location that an annotation legend will appear for
   *              a particular structure.
   * Arguments:
   *     1. legend
   *        The legend to place.
   */
  virtual void placeLegend( Postscript_Legend* legend );

  /*
   * Name:        writeBackbone
   * Description: Write the backbone of the structure.
   *              The backbone is usually a sequence of lines, but not
   *              necessarily. The nucleotides themselves may function as a
   *              backbone, if so they and their labels will be written in this
   *              method, and the writeNucleotidesAndLabels method will not be
   *              implemented.
   * Arguments:
   *     1. out
   *        The output stream to the output file currently being written.
   */
  virtual void writeBackbone( ofstream &out );

  /*
   * Name:        writeBaseAnnotations
   * Description: Write an array of Postscript colors that correspond to each
   *              base in the structure, if the structure is annotated.
   *              Also writes the legend corresponding to these colors.
   * Arguments:
   *     1. out
   *        The output stream to the output file currently being written.
   */
  void writeBaseAnnotations( ofstream &out );

  /*
   * Name:        writeBasesArray
   * Description: Write the array of bases.
   */
  void writeBasesArray( ofstream &out );

  // Method inherited from Postscript_Image.h; see that file for description.
  virtual void writeDescriptor( ofstream &out );

  /*
   * Name:        writeNucleotidesAndLabels
   * Description: Write the nucleotides and labels onto the backkbone.
   *              The nucleotide labels themselves may function as a backbone,
   *              if so they will be written in the writeBackbone method, and
   *              this method will not be implemented.
   * Arguments:
   *     1. out
   *        The output stream to the output file currently being written.
   */
  virtual void writeNucleotidesAndLabels( ofstream &out );

  /*
   * Name:        writePair
   * Description: Write structure-specific Postscript for a particular pair.
   * Arguments:
   *     1. out
   *        The output stream to the output file currently being written.
   */
  virtual void writePair( ofstream &out );

  /*
   * Name:        writePairings
   * Description: Write the pairings between bases of the structure.
   * Arguments:
   *     1. out
   *        The output stream to the output file currently being written.
   */
  void writePairings( ofstream &out );

  /*
   * Name:        writePairingsArray
   * Description: Write the array of pairings as a variable.
   */
  virtual void writePairingsArray( ofstream &out );

  /*
   * Name:        writePairingsVariables
   * Description: Write specific variables holding numbers for particular types
   *              of pairings.
   */
  void writePairingsVariables( ofstream &out );

  // Method inherited from Postscript_Object.h; see that file for description
  virtual void writeSpecificImageType( ofstream &out );

 protected:
  // Protected variables

  // Boolean flags that determine what type of annotation, if any, is applied.
  bool isProbabilityAnnotated;
  bool isSHAPEAnnotated;

  // Vector of pairings.
  vector<int> pairings;

  // Vector that holds the number of pairs and number of pseudoknotted pairs.
  vector<int> pairTypes;

  // The structure number that's being drawn.
  int structureNumber;

  // Boolean that tells whether pairs should be drawn as thick lines (true) or
  // thin lines (false)
  bool useThickPairs;

 private:
  // Object that handles base annotation, if necessary.
  Postscript_Annotation_Handler* annotator;

  // Object that holds a probability annotation legend, if necessary.
  Postscript_Probability_Legend* probabilityLegend;

  // Object that holds a SHAPE annotation legend, if necessary.
  Postscript_SHAPE_Legend* shapeLegend;

  // Code that tells what colorization, if any, can be applied to a structure.
  // 0 = No colorization
  // 1 = Nucleotide colorization only
  // 2 = Pair colorization only
  // 3 = Both nucleotide and pair colorization
  int typeColorCode;
};

#endif /* POSTSCRIPT_STRUCTURE_H */
