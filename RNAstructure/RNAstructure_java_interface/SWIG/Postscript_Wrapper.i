%module PostscriptProxy;

using namespace std;
%include "std_string.i"

%import "DotPlotHandler.i"

%{
#include "../../src/PostscriptWriting/Postscript_Wrapper.h"
%}

%include "../../src/PostscriptWriting/Postscript_Wrapper.h"