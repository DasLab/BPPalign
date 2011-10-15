/*
 * Implementation for a program that parses Unix-like command lines to determine
 * the data being given for input
 * (c) 2008 Jessica Reuter
 */

#include "ParseCommandLine.h"

// Constructor
ParseCommandLine::ParseCommandLine( int num, char* args[], int params,
				    int noneFlags, const char* noneParams[],
				    int paramFlags, const char* flagParams[] ) {
  argc = num;
  argv = args;

  parameters = params;
  numFlagsNoParams = noneFlags;
  numFlagsParams = paramFlags;

  flagsNoParams = noneParams;
  flagsParams = flagParams;

  paramIndices = 0;
}

// Destructor
ParseCommandLine::~ParseCommandLine() {
  if( paramIndices != 0 ) { delete[] paramIndices; }
}

// Check that the command line is structured correctly
bool ParseCommandLine::checkLine( const char* helpFlags ) {
  // Check for a help request, if one exists, no parsing is necessary.
  char helpFlagsArray[100];
  strcpy( helpFlagsArray, helpFlags );

  char* pointer = strtok( helpFlagsArray, " " );
  while( pointer != NULL ) {
    for( int i = 1; i < argc; i++ ) {
      if( strcmp( pointer, argv[i] ) == 0 ) {
        usage();
        return false;
      }
    }
    pointer = strtok( NULL, " " );
  }

  // Create and initialize data structures
  paramIndices = new int[parameters];

  for( int i = 0; i < numFlagsNoParams; i++ ) {
    flagIndices[flagsNoParams[i]] = -1;
  }

  for( int i = 0; i < numFlagsParams; i++ ) {
    flagIndices[flagsParams[i]] = -1;
  }

  // Create a numerical model of the command line
  int* model = new int[argc-1];

  // Distinguish between general argument types
  // 0: param or option, 1: flag
  for( int i = 1; i < argc; i++ ) {
    if( strncmp( argv[i], "-", 1 ) == 0 ) {
      if( argv[i][1] == '.' || isdigit( argv[i][1] ) ) { model[i-1] = 0; }
      else { model[i-1] = 1; }
    } else { model[i-1] = 0; }
  }

  // Determine specific argument types
  // 0: parameter, 1: flag without parameter, 2: flag with parameter, 3: option
  for( int i = 0; i < argc - 1; i++ ) {
    if( model[i] == 0 ) {
      for( int j = 0; j < numFlagsParams; j++ ) {
	if( strcmp( flagsParams[j], argv[i] ) == 0 ) {
	  model[i] = 3;
	  break;
	}
      }
    } else {
      for( int j = 0; j < numFlagsParams; j++ ) {
	if( strcmp( flagsParams[j], argv[i+1] ) == 0 ) {
	  model[i] = 2;
	  break;
	}
      }
    }
  }

  // Scan the model for completeness and correctness
  int paramCounter = 0;

  for( int i = 0; i < argc - 1; i++ ) {
    if( model[i] == 0 ) {
      if( paramCounter < parameters ) {
	paramIndices[paramCounter] = i + 1;
      }
      paramCounter++;
    }
    else if( model[i] == 1 || model[i] == 2 ) {
      bool ok = false;

      for( int j = 0; j < numFlagsNoParams; j++ ) {
	if( strcmp( argv[i+1], flagsNoParams[j] ) == 0 ) {
	  ok = true;
	  break;
	}
      }

      if( !ok ) {
	for( int j = 0; j < numFlagsParams; j++ ) {
	  if( strcmp( argv[i+1], flagsParams[j] ) == 0 ) {
	    ok = true;
	    break;
	  }
	}
      }

      if( !ok ) {
	cerr << "Flag " << argv[i+1] << " does not exist." << endl;
	delete[] model;
	return false;
      } else if( model[i] == 2 && model[i+1] != 3 ) {
	cerr << "Option missing for flag: " << argv[i+1] << endl;
	delete[] model;
	return false;
      }

      if( ok ) {
        const char* nextFlag = referencePointer( argv[i+1] );
        flagIndices[nextFlag] = i + 1;
      }
    }
  }

  if( paramCounter != parameters ) {
    cerr << "Incorrect number of required parameters given." << endl;
    cerr << mainUsage << endl;
    cerr << "Use any of the following options to get a help message: "
	 << helpFlags << "." << endl;
    delete[] model;
    return false;
  }

  delete[] model;
  return true;
}

// Check if the command line contains a particular flag.
bool ParseCommandLine::contains( const char* flag ) {
  map<const char*, int>::iterator iterator = flagIndices.find( flag );
  if( iterator != flagIndices.end() ) { return iterator->second != -1; }
  else return false;
}

// Check if the command line contains one (or more) of a group of flags
bool ParseCommandLine::containsInGroup( const char* flagGroup ) {
  char group[100];
  strcpy( group, flagGroup );

  char* pointer = strtok( group, " " );
  while( pointer != NULL ) {
    const char* safePointer = referencePointer( pointer );
    if( contains( safePointer ) ) { return true; }
    pointer = strtok( NULL, " " );
  }

  return false;
}

// Get an option that immediately follows a flag
char* ParseCommandLine::getOption( const char* flag ) {
  if( contains( flag ) ) {
    return argv[ flagIndices.find( flag )->second + 1 ];
  } else return NULL;
}

// Get a required parameter
char* ParseCommandLine::getParameter( int index ) {
  if( index >= 1 && index <= parameters ) {
    return argv[ paramIndices[index - 1] ];
  } else return NULL;
}

// Ensure that a reference to a command line element is consistent.
const char* ParseCommandLine::referencePointer( const char* pointer ) {
  const char* reference = 0;
  for( int i = 0; i < numFlagsNoParams; i++ ) {
    if( strcmp( pointer, flagsNoParams[i] ) == 0 ) {
      reference = flagsNoParams[i];
    }
  }

  if( reference == 0 ) {
    for( int i = 0; i < numFlagsParams; i++ ) {
      if( strcmp( pointer, flagsParams[i] ) == 0 ) {
        reference = flagsParams[i];
      }
    }
  }

  return reference;
}

// Set an option from the command line. This function is only used for char*
// strings; a template is used for all number types.
bool ParseCommandLine::setOption( const char*& variable, const char* flagGroup,
                                  const char* message, bool file ) {
  char group[100];
  strcpy( group, flagGroup );

  char* flag = strtok( group, " " );
  while( flag != NULL ) {
    const char* safeFlag = referencePointer( flag );
    char* option = getOption( safeFlag );

    if( option != NULL ) {
      variable = option;

      // If string is meant to be an existing file name, check if it exists  
      if( file == true && ifstream( option ) == NULL ) {
        cerr << message << endl;
        return false;
      } else return true;
    }

    flag = strtok( NULL, " " );
  }

  return true;
}

// Set the usage strings
void ParseCommandLine::setUsageStrings( const char* newUsageString,
					const char* newParamString,
					const char* flagNoParamString,
					const char* flagParamString ) {
  mainUsage = newUsageString;
  paramUsage = newParamString;
  flagUsage1 = flagNoParamString;
  flagUsage2 = flagParamString;
}

// Print a usage message
void ParseCommandLine::usage() {
  cout << mainUsage << endl
       << "Short flags are case-insensitive, "
       << "and grouping of flags is not allowed." << endl
       << "Long options are case-sensitive." << endl << endl;

  cout << "=============================" << endl
       << "==== Required Parameters ====" << endl
       << "=============================" << endl
       << paramUsage << endl << endl;

  if( numFlagsParams != 0 ) {
    cout << "======================================" << endl
	 << "==== Option Flags With Parameters ====" << endl
	 << "======================================" << endl
	 << "All parameters must follow their associated flag directly."
	 << endl
	 << "Failure to do so may result in aberrant program behavior."
	 << endl << endl << flagUsage2 << endl << endl;
  }

  if( numFlagsNoParams != 0 || ( strlen( flagUsage1.c_str() ) > 0 ) ) {
    cout << "=========================================" << endl
	 << "==== Option Flags Without Parameters ====" << endl
	 << "=========================================" << endl
	 << flagUsage1 << endl;
  }
}
