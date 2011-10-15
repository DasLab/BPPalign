%module TwoRNAProxy;

%{
#include "../../RNA_class/TwoRNA.h"
%}

%import "../../RNA_class/thermodynamics.h"
%import "../../RNA_class/RNA.h"

%include "std_string.i"
%include "../../RNA_class/TwoRNA.h"