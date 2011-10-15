/*
 * An implementation file for a class that handles an object written in
 * Postscript.
 *
 * (c) 2010  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "Postscript_Object.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Object::Postscript_Object() {

  // Set the error to false.
  error = false;

}

///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
Postscript_Object::~Postscript_Object() {}

///////////////////////////////////////////////////////////////////////////////
// Escape "(" and ")" in strings where necessary so strings with these
// characters will display correctly in Postscript. Also, trim whitespace.
// If the user requests it, truncate the string.
///////////////////////////////////////////////////////////////////////////////
string Postscript_Object::escapeAndTrim( string data, int mode ) {

  // Create a string from the string given.
  string dataString( data );

  // Replace all "(" characters.
  size_t position1 = dataString.find( "(" );
  while( position1 != string::npos ) {
    dataString.replace( position1, 1, "\\(" );
    position1 = dataString.find( "(", position1 + 3 );
  }

  // Replace all ")" characters.
  size_t position2 = dataString.find( ")" );
  while( position2 != string::npos ) {
    dataString.replace( position2, 1, "\\)" );
    position2 = dataString.find( ")", position2 + 3 );
  }

  // Remove all leading and trailing whitespace.
  const char whites[] = { ' ', '\n', '\t', '\r' };
  int whiteLength = 4;

  for( int i = 0; i < whiteLength; i++ ) {
    string::size_type position3 = dataString.find_last_not_of( whites[i] );
    if( position3 != string::npos ) {
      dataString.erase( position3 + 1 );
      position3 = dataString.find_first_not_of( whites[i] );
      if( position3 != string::npos ) { dataString.erase( 0, position3 ); }
    } else { dataString.erase( dataString.begin(), dataString.end() ); }
  }

  // If truncation has been specified, truncate the string.
  if( ( mode != NO_TRUNC ) && ( dataString.length() > 40 ) ) {

    // If the mode is to truncate the beginning, treat it as a file name.
    if( mode == TRUNC_BEGIN ) {
      string::size_type forwardSlash = dataString.find_last_of( '/' );
      string::size_type backSlash = dataString.find_last_of( '\\' );
      string file = "";

      if( forwardSlash != string::npos ) {
	string forwardString = dataString.substr( forwardSlash );
	dataString = "FILEPATH" + forwardString;
      }

      else if( backSlash != string::npos ) {
	string backString = dataString.substr( backSlash );
	dataString = "FILEPATH" + backString;
      }
    }

    // If the mode is to truncate the end, treat it as a simple string.
    else if( mode == TRUNC_END ) {
      string shortString = dataString.substr( 0, 40 );
      dataString = shortString + "...";
    }
  }

  // Return the escaped and trimmed (and possibly truncated) string.
  return dataString;

}

///////////////////////////////////////////////////////////////////////////////
// Check whether this object has encountered an error.
///////////////////////////////////////////////////////////////////////////////
bool Postscript_Object::isError() {

  return error;

}
