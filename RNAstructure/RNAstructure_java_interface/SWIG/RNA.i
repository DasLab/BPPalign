%module RNAProxy;

%typemap(javain) TProgressDialog *e "getCPtrAndAddReference($javainput)"

%typemap(javacode) RNA %{
  private TProgressDialog dialogReference;
  private long getCPtrAndAddReference(TProgressDialog dialog) {
    dialogReference = dialog;
    return TProgressDialog.getCPtr(dialog);
  }
%}

%{
#include "../../RNA_class/RNA.h"
%}

%import "TProgressDialog.h"
%import "../../RNA_class/thermodynamics.h"

%include "std_string.i"

%ignore RNA::GetStructure();
%include "../../RNA_class/RNA.h"