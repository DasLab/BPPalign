/*
 * A header file for a class that writes annotations on a Postscript RNA
 * secondary structure.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#ifndef POSTSCRIPT_ANNOTATION_HANDLER_H
#define POSTSCRIPT_ANNOTATION_HANDLER_H

#include "Postscript_Piece.h"

class Postscript_Annotation_Handler : public Postscript_Piece {
 public:
  // Public constructor, destructor, and methods

  /*
   * Name:        Constructor
   * Description: Initializes the private data arrays.
   * Arguments:
   *     1. len
   *        The length of the sequence being annotated.
   *     2. structs
   *        The number of structures being annotated.
   *        Default is 0.
   */
  Postscript_Annotation_Handler( int len, int structs = 0 );

  /*
   * Name:        getProbabilityColor
   * Description: Get the Postscript color constant assigned to a particular
   *              nucleotide that shows what color it's annotated, for base
   *              pair probability data. Probability data is structure
   *              dependent, so a two-dimensional vector is necessary here.
   * Arguments:
   *     1. structure
   *        The structure number, one-indexed, to get annotation from.
   *     2. nucleotide
   *        The nucleotide, one-indexed, whose annotation to get
   */
  string getProbabilityColor( int structure, int nucleotide );

  /*
   * Name:        getSHAPEColor
   * Description: Get the Postscript color constant assigned to a particular
   *              nucleotide that shows what color it's annotated, for SHAPE
   *              data. SHAPE data is not structure dependent, so only a one-
   *              dimensional vector is necessary here.
   * Arguments:
   *     1. index
   *        The nucleotide (one-indexed) whose annotation color to get.
   * Returns:
   *     The Postscript color constant assigned, as a string.
   */
  string getSHAPEColor( int index );

  /*
   * Name:        readPartition
   * Description: Read base pairing probability data used to determine color
   *              annotation, from a partition function save file. Base pairing
   *              probabilities are structure dependent.
   * Arguments:
   *     1. file
   *        The name of the partition function save file to read.
   *     2. strand
   *        The strand being annotated (specific structures are necessary). The
   *        strand MUST be dervied from a CT file, so it has structures in it
   *        to analyze base pairing of.
   */
  void readPartition( string file, RNA* strand );

  /*
   * Name:        readSHAPE
   * Description: Read SHAPE data used to determine color annotation, from a
   *              SHAPE text data file. SHAPE data is not structure dependent,
   *              so a structure number is not necessary.
   * Arguments:
   *     1. file
   *        The name of the SHAPE data file to read.
   */
  void readSHAPE( string file );

 private:
  // Array of annotation colors for probability data
  vector< vector<char> > probabilityAnnotations;

  // Array of annotation colors for SHAPE data
  vector<char> shapeAnnotations;

  // The length of the sequence being annotated.
  int length;

  // The number of structures being annotated.
  int structures;
};

#endif /* POSTSCRIPT_ANNOTATION_HANDLER_H */
