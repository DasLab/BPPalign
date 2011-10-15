/*
 * A header file for a class that writes a single circular Postscript structure
 * whose pairs show a comparison between two distinct structures.
 *
 * Copyright 2010 Mathews Lab, University of Rochester
 * Written by Jessica S. Reuter
 */

#ifndef POSTSCRIPT_CIRCLE_STRUCTURES_H
#define POSTSCRIPT_CIRCLE_STRUCTURES_H

#include "Postscript_Circle_Structure.h"

class Postscript_Circle_Structures : public Postscript_Circle_Structure {
 public:
  // Public constructor and destructor

  /*
   * Name:        Constructor
   * Description: Constructor which initializes the RNA data structure for this
   *              class directly, creating the RNA object inside itself.
   * Arguments:
   *     1. in
   *        The name of the file to create an internal RNA object with.
   *     2. out
   *        The name of the Postscript output file to write.
   */
  Postscript_Circle_Structures( string in, string out );

  /*
   * Name:        Destructor
   * Description: Actually has nothing to do, here as a placeholder.
   */
  ~Postscript_Circle_Structures();

  /*
   * Name:        addSecondStructure
   * Description: Creates data structures necessary for the second structure
   *              which is overlaid over the first (base).
   * Arguments:
   *     1. second
   *        The name of the file to create a second internal RNA object with.
   *     2. number
   *        The structure in the second strand to be compared.
   */
  void addSecondStructure( string second, int number );

  /*
   * Name:        getNumPairs
   * Description: Get the number of unique pairs among the two structures.
   *              This method overrides getNumPairs() in
   *              Postscript_Structure.h.
   * Returns:
   *     The number of pairs among the structures.
   */
  int getNumPairs();

  /*
   * Name:        setColors
   * Description: Set the colors that distinguish the types of pairs in this
   *              combined structure.
   * Arguments:
   *     1. both
   *        The color used for pairs found in both structures.
   *     2. bothN
   *        The English name of the color used for pairs in both structures.
   *     3. one
   *        The color used for pairs found only in structure 1 (strand)
   *     4. oneN
   *        The English name of the color used for pairs in structure 1 only.
   *     5. two
   *        The color used for pairs found only in structure 2 (strand2)
   *     6. twoN
   *        The English name of the color used for pairs in structure 2 only.
   */
  void setColors( const string& both, const string& bothN,
		  const string& one, const string& oneN,
		  const string& two, const string& twoN );

  /*
   * Name:        setIdentifiers
   * Description: Set specialized names that identify structures 1 and 2, if
   *              they should not be called "Structure 1" and/or "Structure 2"
   *              in the final file.
   * Arguments:
   *     1. id1
   *        The new ID for structure 1
   *     2. id2
   *        The new ID for structure 2
   */
  void setIdentifiers( const string id1, const string id2 );

 protected:
  // Protected methods

  // Method inherited from Postscript_Object.h; see file for description
  virtual void writeDescriptor( ofstream &out );

  // Method inherited from Postscript_Structure.h; see file for description
  void writePairingsArray( ofstream &out );

 private:
  // Private variables

  // The comment for the original strand before a new strand was added.
  string comment2;

  // The file name of the original RNA strand before a new strand was added.
  string filename2;

  // The vector of pairings in the original RNA strand before a new strand was
  // added.
  vector<int> pairings2;

  // The number of unique pairs in the two structures.
  int uniquePairs;

  // References to the Postscript colors from Postscript_Object.h that will be
  // used for pairs in both structures, pairs only in structure 1, and pairs
  // only in structure 2.
  string bothColor;
  string oneColor;
  string twoColor;

  // The English names of the colors used in this combined structure.
  string bothName;
  string oneName;
  string twoName;

  // The general English identifiers for structures 1 and 2, not file names.
  string structureID1;
  string structureID2;
};

#endif /* POSTSCRIPT_CIRCLE_STRUCTURES_H */
