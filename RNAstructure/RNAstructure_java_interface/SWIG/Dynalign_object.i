%module Dynalign_objectProxy;

// Double typemaps for PFPRECISION
%typemap( jni ) PFPRECISION "jdouble"
%typemap( jtype ) PFPRECISION "double"
%typemap( jstype ) PFPRECISION "double"
%typemap( javain ) PFPRECISION "$javainput"
%typemap( in ) PFPRECISION %{ if ($1) delete $1; $1 = new double($input); %}
%typemap( out ) PFPRECISION %{ if ($1) $1; %}
%typemap( javaout ) PFPRECISION { return $jnicall; }

%{
#include "../../RNA_class/Dynalign_object.h"
%}

%import "TProgressDialog.h"
%import "../../RNA_class/TwoRNA.h"

%include "std_string.i"
%include "../../RNA_class/Dynalign_object.h"
