include ../../compiler.h

OBJ_FILES=../../src/algorithm.o ../../src/rna_library.o ../../src/structure.o globals.o \
	../../src/siRNAfilter.o ../../src/thermo.o ../../src/random.o ../../src/stochastic.o \
	../../src/intermolecular.o ../../src/alltrace.o ../../src/pclass.o ../../src/TProgressDialog.o \
	../../src/pfunction.o ../../src/forceclass.o ../../src/stackclass.o \
	../../src/alltrace_intermolecular.o ../../src/arrayclass.o ../../src/dotarray.o \
	../../src/stackstruct.o oligowalk.o 

# rules
#####################
all: 
	@echo 'Type make OligoWalk to make executable file in exe/'

OligoWalk: ../exe/OligoWalk
../exe/OligoWalk: ${OBJ_FILES}
	${LINK} ${OBJ_FILES}
	@echo "------"
	@echo "Installation complete."
	@echo "------"
	@echo 'A perl script, siRNAWalk.pl, is avaible to run OligoWalk and predict siRNA automatically'
	@echo 'Prerequisite: libsvm. (copy svm-predict and svm-scale to directory svm/)'


# Object Files 
##############################################
oligowalk.o: oligowalk.cpp ../../src/intermolecular.h ../../src/rna_library.h ../../src/algorithm.h ../../src/structure.h \
	../../src/thermo.h ../../src/pclass.h ../../src/alltrace.h ../../src/stochastic.h ../../src/siRNAfilter.h \
	../../src/TProgressDialog.h ../../src/platform.h ../../src/defines.h globals.h

../../src/intermolecular.o: ../../src/intermolecular.cpp ../../src/intermolecular.h ../../src/siRNAfilter.h \
	../../src/pclass.h ../../src/alltrace.h ../../src/stochastic.h ../../src/defines.h ../../src/TProgressDialog.h \
	../../src/platform.h ../../src/rna_library.h ../../src/structure.h ../../src/algorithm.h 

../../src/siRNAfilter.o: ../../src/siRNAfilter.cpp ../../src/siRNAfilter.h ../../src/algorithm.h ../../src/rna_library.h \
	../../src/structure.h ../../src/TProgressDialog.h 

../../src/pclass.o: ../../src/pclass.cpp ../../src/pclass.h ../../src/rna_library.h ../../src/algorithm.h \
	../../src/structure.h ../../src/TProgressDialog.h ../../src/defines.h ../../src/pfunction.h

../../src/pfunction.o: ../../src/pfunction.cpp ../../src/defines.h ../../src/pfunction.h ../../src/rna_library.h \
	../../src/structure.h ../../src/platform.h

../../src/stochastic.o: ../../src/stochastic.cpp ../../src/stochastic.h ../../src/random.h ../../src/pclass.h \
	../../src/platform.h ../../src/structure.h ../../src/rna_library.h ../../src/algorithm.h ../../src/TProgressDialog.h \
	../../src/defines.h ../../src/pclass.h 

../../src/alltrace.o: ../../src/alltrace.cpp ../../src/alltrace.h ../../src/defines.h ../../src/pclass.h \
	../../src/rna_library.h ../../src/algorithm.h ../../src/structure.h ../../src/TProgressDialog.h ../../src/platform.h 

../../src/algorithm.o:  ../../src/algorithm.cpp ../../src/algorithm.h ../../src/defines.h ../../src/structure.h \
	../../src/platform.h ../../src/TProgressDialog.h ../../src/rna_library.h

../../src/rna_library.o: ../../src/rna_library.cpp  ../../src/rna_library.h ../../src/platform.h ../../src/structure.h 

../../src/structure.o: ../../src/structure.cpp ../../src/structure.h ../../src/platform.h ../../src/defines.h 

globals.o: globals.cpp globals.h 

../../src/thermo.o: ../../src/thermo.cpp ../../src/thermo.h ../../src/algorithm.h ../../src/defines.h \
	../../src/TProgressDialog.h ../../src/rna_library.h ../../src/structure.h 

../../src/random.o: ../../src/random.cpp ../../src/random.h

../../src/TProgressDialog.o: ../../src/TProgressDialog.cpp ../../src/TProgressDialog.h

../../src/forceclass.o: ../../src/forceclass.cpp ../../src/forceclass.h

../../src/stackclass.o: ../../src/stackclass.cpp ../../src/stackclass.h

../../src/arrayclass.o: ../../src/arrayclass.cpp ../../src/arrayclass.h

../../src/alltrace_intermolecular.o: ../../src/alltrace_intermolecular.cpp ../../src/alltrace_intermolecular.h \
	../../src/defines.h ../../src/structure.h ../../src/rna_library.h ../../src/pclass.h

../../src/stackstruct.o: ../../src/stackstruct.cpp ../../src/stackstruct.h

../../src/dotarray.o: ../../src/dotarray.cpp ../../src/dotarray.h


##############################################

#clean the files
clean:
	rm -f *.o
	rm -f ../../src/*.o

realclean: clean
	rm -f ../exe/OligoWalk
