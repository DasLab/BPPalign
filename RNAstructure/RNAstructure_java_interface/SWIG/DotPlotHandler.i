%module DotPlotHandlerProxy;

%{
#include "../../src/DotPlotHandler.h"
%}

%include "std_string.i"

%ignore DotPlot_Constants::DEFAULT_ENTRIES;
%ignore DotPlot_Constants::MIN_COLORS;
%ignore DotPlot_Constants::MAX_COLORS;

%ignore DotPlot_Constants::TYPE_DYNALIGN1;
%ignore DotPlot_Constants::TYPE_DYNALIGN2;
%ignore DotPlot_Constants::TYPE_ENERGY;
%ignore DotPlot_Constants::TYPE_PROBABILITY;
%ignore DotPlot_Constants::TYPE_UNDEFINED;

%include "../../src/DotPlotHandler.h"