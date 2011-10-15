t_structure.h anbd t_structure.cpp are totally redundant files to the structure.h and structure.cpp
in this directory.  These were added so that Visual Studio could compile PARTS.  PARTS also uses
RNAstructure/src/structure.h and RNAstructure/src/structure.cpp and Visual Syudio does not allow
a program to have two objects with same name even if they are in different directories!
