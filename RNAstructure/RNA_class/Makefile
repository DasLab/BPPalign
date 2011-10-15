include ../compiler.h
include ../library_defines.h

##########
## Use the Java TProgressDialog and ProgressMonitor if libraries
## are to be integrated with Java
##########

TPROGRESSDIR = ../src/
OBSERVINGBAR = ../src/observingtextprogressbar.o
PROGRESSMONITOR =

Set_Java =
ifeq (${JAVA},yes)
	${CXXFLAGS} := -D_JAVA_GUI
	TPROGRESSDIR = ../RNAstructure_java_interface/SWIG/
	OBSERVINGBAR =
	PROGRESSMONITOR = ../RNAstructure_java_interface/SWIG/ProgressMonitor.o
endif

##########
## Define Internal Macros
##########

# Common files for the RNA library.
RNA_FILES = \
	RNA.o \
	../src/algorithm.o \
	../src/alltrace.o \
	../src/arrayclass.o \
	../src/dotarray.o \
	../src/draw.o \
	../src/forceclass.o \
	../src/MaxExpect.o \
	../src/MaxExpectStack.o \
	../src/outputconstraints.o \
	../src/pfunction.o \
	../src/probknot.o \
	../src/random.o \
	../src/rna_library.o \
	../src/stackclass.o \
	../src/stackstruct.o \
	../src/stochastic.o \
	../src/structure.o \
	thermodynamics.o \
	${TPROGRESSDIR}TProgressDialog.o \
	${PROGRESSMONITOR}

# Common files for the TwoRNA library.
TWORNA_FILES = \
	${RNA_FILES} \
	TwoRNA.o \
	../src/bimol.o

# Common files for the HybridRNA library.
HYBRID_FILES = \
	${TWORNA_FILES} \
	HybridRNA.o

# Common files for both versions of the Dynalign library.
DYNALIGN_FILES = \
	${TWORNA_FILES} \
	Dynalign_object.o \
	../src/align_services.o \
	../src/alignment_hmm_model.o \
	../src/aln_env_utils.o \
	../src/configfile.o \
	../src/dynalignarray.o \
	../src/dynalignheap.o \
	../src/dynalignstackclass.o \
	../src/hmm_arrays.o \
	../src/hmm_probs.o \
	../src/misc_decoding_services.o \
	../src/observable.o \
	../src/observer.o \
	${OBSERVINGBAR} \
	../src/process_output.o \
	../src/varray.o \
	../src/wendarray.o \
	../src/xlog_math.o

# Common files for the Dynalign library when compiled serially.
DYNALIGN_SERIAL_FILES = \
	${DYNALIGN_FILES} \
	../src/dynalign.o

# Common files for the Dynalign library when compiling for SMP.
DYNALIGN_SMP_FILES = \
	${DYNALIGN_FILES} \
	../src/dynalign-smp.o \
	../src/rank.o \
	../src/rankconsumer.o \
	../src/rankmanager.o \
	../src/rankproducer.o

# Common files for the Oligo library.
OLIGO_FILES = \
	${RNA_FILES} \
	Oligowalk_object.o \
	../src/alltrace_intermolecular.o \
	../src/intermolecular.o \
	../src/OligoScreenCalc.o \
	../src/pclass.o \
	../src/siRNAfilter.o \
	../src/thermo.o

##########
## Define Actions
##########

# Instructions on which make commands to use
all:
	@echo 'All executables and libraries are moved to the "exe" folder.'
	@echo 'Use "make RNA_class" to make "RNA_class."'
	@echo 'Use "make RNA_Library" to make the RNA library.'
	@echo 'Use "make HybridRNA_class" to make "HybridRNA_class."'
	@echo 'Use "make HybridRNA_Library" to make the HybridRNA library.'
	@echo 'Use "make Dynalign_class" to make "Dynalign_class."'
	@echo 'Use "make Dynalign_Library" to make the Dynalign library.'
	@echo 'Use "make DynalignSMP_class" to make "DynalignSMP_class."'
	@echo 'Use "make DynalignSMP_Library" to make the DynalignSMP library.'
	@echo 'Use "make OligoWalk_class" to make "OligoWalk_class."'
	@echo 'Use "make Oligo_Library" to make the Oligo library.'

# Build the RNA class shared library.
RNA_Library: ${RNA_FILES}
	${LIBCOMP} ${RNA_LIBRARY} ${RNA_FILES}
	mv ${RNA_LIBRARY} ../exe/

# Build the HybridRNA class shared library.
HybridRNA_Library: ${HYBRID_FILES}
	${LIBCOMP} ${HYBRID_RNA_LIBRARY} ${HYBRID_FILES}
	mv ${HYBRID_RNA_LIBRARY} ../exe/

# Build the Dynalign shared library for serial calculations.
Dynalign_Library: ${DYNALIGN_SERIAL_FILES}
	${LIBCOMP} ${DYNALIGN_LIBRARY} ${DYNALIGN_SERIAL_FILES}
	mv ${DYNALIGN_LIBRARY} ../exe/

# Build the Dynalign shared library for SMP calculations.
DynalignSMP_Library: ${DYNALIGN_SMP_FILES}
	${LIBCOMP} ${DYNALIGN_SMP_LIBRARY} ${DYNALIGN_SMP_FILES}
	mv ${DYNALIGN_SMP_LIBRARY} ../exe/

# Build the Oligo shared library.
Oligo_Library: $(OLIGO_FILES)
	${LIBCOMP} ${OLIGO_LIBRARY} ${OLIGO_FILES}
	mv $(OLIGO_LIBRARY) ../exe/

# Build the RNA class example program for RNA
RNA_class: ${RNA_FILES} RNA_class.o 
	${LINK} ${RNA_FILES} RNA_class.o
	mv RNA_class ../exe/

# Build the HybridRNA class example program for HybridRNA
HybridRNA_class: ${HYBRID_FILES} HybridRNA_class.o
	${LINK} ${HYBRID_FILES} HybridRNA_class.o
	mv HybridRNA_class ../exe/

# Build the Dynalign class example program for Dynalign_object
Dynalign_class: ${DYNALIGN_SERIAL_FILES} Dynalign_class.o
	${LINK} ${DYNALIGN_SERIAL_FILES} Dynalign_class.o
	mv Dynalign_class ../exe/

# Build the Dynalign class example program for Dynalign_object
DynalignSMP_class: ${DYNALIGN_SMP_FILES} DynalignSMP_class.o
	${LINK} ${DYNALIGN_SMP_FILES} DynalignSMP_class.o -lpthread
	mv DynalignSMP_class ../exe/

#Build the OligoWalk_class example program for Oligowalk_object
OligoWalk_class: ${OLIGO_FILES} OligoWalk_class.o
	${LINK} ${OLIGO_FILES} OligoWalk_class.o 
	mv Oligowalk_class ../exe/

##########
## Rules for Target-Specific Object Files
##########

RNA_class.o: RNA_class.cpp ../src/structure.h 

HybridRNA_class.o: HybridRNA_class.cpp ../src/structure.h

Dynalign_class.o: Dynalign_class.cpp

OligoWalk_class.o: OligoWalk_class.cpp

##########
## Rules for RNA Class Object Files (RNA_FILES)
##########

RNA.o: RNA.cpp RNA.h \
	../src/algorithm.cpp ../src/algorithm.h \
	../src/alltrace.h \
	../src/arrayclass.h \
	../src/bimol.h \
	../src/defines.h \
	../src/dotarray.h \
	../src/draw.h \
	../src/forceclass.h \
	../src/MaxExpect.h \
	../src/pfunction.h \
	../src/probknot.h \
	../src/platform.h \
	../src/random.h \
	../src/rna_library.h \
	../src/stackclass.h \
	../src/stackstruct.h \
	../src/stochastic.h \
	../src/structure.h \
	${TPROGRESSDIR}TProgressDialog.h \
	thermodynamics.cpp thermodynamics.h

../src/algorithm.o: ../src/algorithm.cpp ../src/algorithm.h \
	../src/arrayclass.h \
	../src/defines.h \
	../src/dotarray.h \
	../src/forceclass.h \
	../src/platform.h \
	../src/rna_library.h \
	../src/stackclass.h \
	../src/stackstruct.h \
	../src/structure.h \
	${TPROGRESSDIR}TProgressDialog.h

../src/alltrace.o: ../src/alltrace.cpp ../src/alltrace.h \
	../src/defines.h \
	../src/structure.h

../src/arrayclass.o: ../src/arrayclass.cpp ../src/arrayclass.h \
	../src/defines.h

../src/dotarray.o: ../src/dotarray.cpp ../src/dotarray.h \
	../src/defines.h

../src/draw.o: ../src/draw.cpp ../src/draw.h \
	../src/substructure.cpp

../src/forceclass.o: ../src/forceclass.cpp ../src/forceclass.h

../src/MaxExpect.o: ../src/MaxExpect.cpp ../src/MaxExpect.h ../src/defines.h

../src/MaxExpectStack.o:  ../src/MaxExpectStack.cpp ../src/MaxExpectStack.h

../src/outputconstraints.o: ../src/outputconstraints.cpp \
	../src/outputconstraints.h

../src/probknot.o: ../src/probknot.cpp ../src/probknot.h

../src/pfunction.o: ../src/pfunction.cpp ../src/pfunction.h

../src/random.o: ../src/random.cpp ../src/random.h

../src/rna_library.o: ../src/rna_library.cpp ../src/rna_library.h \
	../src/defines.h \
	../src/platform.h \
	../src/structure.h

../src/stackclass.o: ../src/stackclass.cpp ../src/stackclass.h \
	../src/defines.h

../src/stackstruct.o: ../src/stackstruct.cpp ../src/stackstruct.h

../src/stochastic.o: ../src/stochastic.cpp ../src/stochastic.h

../src/structure.o: ../src/structure.cpp ../src/structure.h \
	../src/defines.h \
	../src/platform.h

thermodynamics.o: thermodynamics.cpp thermodynamics.h

../src/TProgressDialog.o: ../src/TProgressDialog.cpp ../src/TProgressDialog.h

# The TProgressDialog used when libraries are made for integration with Java
# Also, the progress monitor
../RNAstructure_java_interface/SWIG/TProgressDialog.o: \
	../RNAstructure_java_interface/SWIG/TProgressDialog.cpp \
	../RNAstructure_java_interface/SWIG/TProgressDialog.h

../RNAstructure_java_interface/SWIG/ProgressMonitor.o: \
	../RNAstructure_java_interface/SWIG/ProgressMonitor.cpp \
	../RNAstructure_java_interface/SWIG/ProgressMonitor.h


##########
## Rules for TwoRNA Class Object Files (TWORNA_FILES)
##########

TwoRNA.o: TwoRNA.cpp TwoRNA.h \
	../src/algorithm.cpp ../src/algorithm.h \
	../src/arrayclass.h \
	../src/bimol.h \
	../src/defines.h \
	../src/dotarray.h \
	../src/forceclass.h \
	../src/pfunction.h \
	../src/platform.h \
	RNA.cpp RNA.h \
	../src/rna_library.h \
	../src/stackclass.h \
	../src/stackstruct.h \
	../src/structure.h \
	thermodynamics.cpp thermodynamics.h \
	${TPROGRESSDIR}TProgressDialog.h

../src/bimol.o: ../src/bimol.cpp ../src/bimol.h

##########
## Rules for Hybrid RNA class Object Files (HYBRID_FILES)
##########

HybridRNA.o: HybridRNA.cpp HybridRNA.h \
	../src/algorithm.cpp ../src/algorithm.h \
	../src/arrayclass.h \
	../src/bimol.h  \
	../src/defines.h \
	../src/dotarray.h \
	../src/forceclass.h \
	../src/pfunction.h \
	../src/platform.h \
	RNA.cpp RNA.h \
	../src/rna_library.h \
	../src/stackclass.h \
	../src/stackstruct.h \
	../src/structure.h \
	thermodynamics.cpp thermodynamics.h \
	${TPROGRESSDIR}TProgressDialog.h \
	TwoRNA.cpp TwoRNA.h

#########
## Rules for Dynalign Object Files (DYNALIGN_FILES)
#########

Dynalign_object.o: Dynalign_object.cpp Dynalign_object.h \
	../src/align_services.h \
	../src/hmm_arrays.h \
	../src/platform.h

../src/dynalign.o: ../src/dynalign.cpp ../src/dynalign.h \
	../src/algorithm.h \
	../src/arrayclass.h \
	../src/defines.h \
	../src/dynalignarray.h \
	../src/dynalignheap.h \
	../src/dynalignstackclass.h \
	../src/forceclass.h \
	../src/platform.h \
	../src/rna_library.h \
	../src/structure.h \
	${TPROGRESSDIR}TProgressDialog.h \
	../src/varray.h \
	../src/wendarray.h

../src/align_services.o: ../src/align_services.cpp \
	../src/align_services.h

../src/alignment_hmm_model.o: ../src/alignment_hmm_model.cpp \
	../src/alignment_hmm_model.h

../src/aln_env_utils.o: ../src/aln_env_utils.cpp ../src/aln_env_utils.h

../src/configfile.o: ../src/configfile.cpp ../src/configfile.h

../src/dynalignarray.o: ../src/dynalignarray.cpp ../src/dynalignarray.h \
	../src/defines.h \
	../src/dynalign.h

../src/dynalignheap.o: ../src/dynalignheap.cpp ../src/dynalignheap.h \
	../src/defines.h

../src/dynalignstackclass.o: ../src/dynalignstackclass.cpp \
	../src/dynalignstackclass.h \
	../src/defines.h

../src/hmm_arrays.o: ../src/hmm_arrays.cpp ../src/hmm_arrays.h

../src/hmm_probs.o: ../src/hmm_probs.cpp ../src/hmm_probs.h

../src/misc_decoding_services.o: ../src/misc_decoding_services.cpp \
	../src/misc_decoding_services.h

../src/observable.o: ../src/observable.cpp ../src/observable.h \
	../src/observer.h

../src/observer.o: ../src/observer.cpp ../src/observer.h

# This text progress bar is only used when the Java GUI isn't required.
../src/observingtextprogressbar.o: ../src/observingtextprogressbar.cpp \
	../src/observingtextprogressbar.h \
	../src/observer.h \
	../src/TProgressDialog.h

../src/process_output.o: ../src/process_output.cpp \
	../src/process_output.h

../src/varray.o: ../src/varray.cpp ../src/varray.h \
	../src/defines.h \
	../src/dynalign.h

../src/wendarray.o: ../src/wendarray.cpp ../src/wendarray.h \
	../src/defines.h

../src/xlog_math.o: ../src/xlog_math.cpp ../src/xlog_math.h

##########
## Rules for DynalignSMP Object Files (DYNALIGN_SMP_FILES)
##########

../src/dynalign-smp.o: ../src/dynalign.cpp ../src/dynalign.h \
	../src/algorithm.h \
	../src/arrayclass.h \
	../src/defines.h \
	../src/dynalignarray.h \
	../src/dynalignheap.h \
	../src/dynalignstackclass.h \
	../src/forceclass.h \
	../src/observingtextprogressbar.h \
	../src/platform.h \
	../src/rankconsumer.h \
	../src/rankmanager.h \
	../src/rankproducer.h \
	../src/rna_library.h \
	../src/structure.h \
	../src/TProgressDialog.h \
	../src/varray.h \
	../src/wendarray.h
	${CXX} ${CXXFLAGS} -DDYNALIGN_SMP -c -o $@ ../src/dynalign.cpp

../src/rank.o: ../src/rank.cpp ../src/rank.h \
	../src/workslice.h \
	../src/workunit.h

../src/rankconsumer.o: ../src/rankconsumer.cpp ../src/rankconsumer.h \
	../src/dynalign.h \
	../src/dynalignarray.h \
	../src/rankmanager.h \
	../src/rna_library.h \
	../src/structure.h \
	../src/varray.h \
	../src/wendarray.h \
	../src/workslice.h

../src/rankmanager.o: ../src/rankmanager.cpp ../src/rankmanager.h \
	../src/observable.h \
	../src/rank.h \
	../src/TProgressDialog.h \
	../src/workslice.h

../src/rankproducer.o: ../src/rankproducer.cpp ../src/rankproducer.h \
	../src/rank.h \
	../src/rankmanager.h \
	../src/workunit.h


##########
## Rules for Oligo Object Files (OLIGO_FILES)
##########

Oligowalk_object.o: Oligowalk_object.cpp Oligowalk_object.h

../src/alltrace_intermolecular.o: ../src/alltrace_intermolecular.cpp \
	../src/alltrace_intermolecular.h

../src/intermolecular.o: ../src/intermolecular.cpp ../src/intermolecular.h \
	../src/siRNAfilter.cpp ../src/siRNAfilter.h

../src/OligoScreenCalc.o: ../src/OligoScreenCalc.cpp ../src/OligoScreenCalc.h

../src/pclass.o: ../src/pclass.cpp ../src/pclass.h

../src/siRNAfilter.o: ../src/siRNAfilter.cpp ../src/siRNAfilter.h

../src/thermo.o: ../src/thermo.cpp ../src/thermo.h

##########
## Object and executable cleanup
## Executable cleanup removes all possible executables.
##########

clean:
	rm -f *.o
	rm -f ../src/*.o

realclean: clean
	rm -f ../exe/RNA_class
	rm -f ../exe/${RNA_LIBRARY}
	rm -f ../exe/HybridRNA_class
	rm -f ../exe/${HYBRID_RNA_LIBRARY}
	rm -f ../exe/Dynalign_class
	rm -f ../exe/${DYNALIGN_LIBRARY}
	rm -f ../exe/OligoWalk_class
	rm -f ../exe/${OLIGO_LIBRARY}
	rm -f ../exeDynalignSMP_class
	rm -f ../exe/$(DYNALIGN_SMP_LIBRARY)

