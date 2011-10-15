##########
## Set Library Names for Text Interfaces
## Options for Windows, Mac, and Linux (default)
##########

ifeq (${OPSYSTEM},Windows)
	RNA_LIBRARY = RNA_Library.dll
	HYBRID_RNA_LIBRARY = HybridRNA_Library.dll
	DYNALIGN_LIBRARY = Dynalign_Library.dll
	DYNALIGN_SMP_LIBRARY = Dynalign_SMP_Library.dll
	OLIGO_LIBRARY = Oligo_Library.dll
	PARTS_LIBRARY = PARTS_Library.dll
	POSTSCRIPT_LIBRARY = Postscript_Library.dll
else ifeq (${OPSYSTEM},Mac)
	RNA_LIBRARY = libRNA.dylib
	HYBRID_RNA_LIBRARY = libHybridRNA.dylib
	DYNALIGN_LIBRARY = libDynalign.dylib
	DYNALIGN_SMP_LIBRARY = libDynalign_SMP.dylib
	OLIGO_LIBRARY = libOligo.dylib
	PARTS_LIBRARY = libPARTS.dylib
	POSTSCRIPT_LIBRARY = libPS.dylib
else
	RNA_LIBRARY = libRNA.so
	HYBRID_RNA_LIBRARY = libHybridRNA.so
	DYNALIGN_LIBRARY = libDynalign.so
	DYNALIGN_SMP_LIBRARY = libDynalign_SMP.so
	OLIGO_LIBRARY = libOligo.so
	PARTS_LIBRARY = libPARTS.so
	POSTSCRIPT_LIBRARY = libPS.so
endif

##########
## Set Library Names for Usage With Java
## Options for Windows, Mac, and Linux (default)
##########

ifeq (${JAVA},yes)
ifeq (${OPSYSTEM},Windows)
	HYBRID_RNA_LIBRARY = HybridRNA_Library_GUI.dll
	DYNALIGN_LIBRARY = Dynalign_Library_GUI.dll
	OLIGO_LIBRARY = Oligo_Library_GUI.dll
	RNASTRUCTURE_LIBRARY = RNAstructure_GUI.dll
else ifeq (${OPSYSTEM},Mac)
	HYBRID_RNA_LIBRARY = libHybridRNA_GUI.dylib
	DYNALIGN_LIBRARY = libDynalign_GUI.dylib
	OLIGO_LIBRARY = libOligo_GUI.dylib
	RNASTRUCTURE_LIBRARY = libRNAstructure_GUI.dylib
else
	HYBRID_RNA_LIBRARY = libHybridRNA_GUI.so
	DYNALIGN_LIBRARY = libDynalign_GUI.so
	OLIGO_LIBRARY = libOligo_GUI.so
	RNASTRUCTURE_LIBRARY = libRNAstructure_GUI.so
endif
endif

##########
## Set any architecture constraints, as necessary.
##########

ifeq (${OPSYSTEM},Mac)
	ARCHITECTURE = -arch x86_64
else ifeq (${OPSYSTEM},Windows)
	ARCHITECTURE = ${LIBFLAGS}
endif

