%module ThermodynamicsProxy;

%{
#include "../../RNA_class/thermodynamics.h"
%}

%ignore Thermodynamics::GetDatatable();
%ignore Thermodynamics::GetEnthalpyTable();
%include "../../RNA_class/thermodynamics.h"