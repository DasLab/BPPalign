include ../compiler.h
include ../library_defines.h

##########
## Define any extensions to flags specific to the GUI
##########

Set_GUI_Flags =
ifeq (${OPSYSTEM},Linux)
	GUI_CFLAGS = ${CXXFLAGS} -I${INCLUDEPATH1} -I${INCLUDEPATH2}
	GUI_LFLAGS = ${LIBFLAGS}
else ifeq (${OPSYSTEM},Mac)
	GUI_CFLAGS = ${CXXFLAGS} -I${INCLUDEPATH1}
	GUI_LFLAGS = ${LIBFLAGS}
else
	GUI_CFLAGS = ${CXXFLAGS} -I${INCLUDEPATH1} -I${INCLUDEPATH2}
	GUI_LFLAGS = ${LIBFLAGS} -shared -I${INCLUDEPATH1} -I${INCLUDEPATH2}
endif

##########
## Define Repeated Commands
##########

SWIG_COMMAND = swig -java -c++ -package source -outdir .
OBJECT_COMMAND = ${CXX} ${GUI_CFLAGS} -c -o
SWIG_OBJECT_COMMAND = ${CXX} ${GUI_CFLAGS} -O -c -o

##########
## Define Actions
##########

# Rule to build the entire RNAstructure program, C++ and Java code
all:
	make buildC compileJava jar JAVA=yes

# Rule to create proxy files with SWIG -- not likely to be needed by end user
swigFiles: clean
	${SWIG_COMMAND} SWIG/Dynalign_object.i
	sed -i 's/Dynalign_objectProxyJNI.SWIGDynalign_objectUpcast(cPtr)/cPtr/g' Dynalign_object.java
	${SWIG_COMMAND} SWIG/Oligowalk_object.i
	sed -i 's/Oligowalk_objectProxyJNI.SWIGOligowalk_objectUpcast(cPtr)/cPtr/g' Oligowalk_object.java
	${SWIG_COMMAND} SWIG/ProgressMonitor.i
	${SWIG_COMMAND} SWIG/TProgressDialog.i
	${SWIG_COMMAND} SWIG/Thermodynamics.i
	${SWIG_COMMAND} SWIG/RNA.i
	sed -i 's/RNAProxyJNI.SWIGRNAUpcast(cPtr)/cPtr/g' RNA.java
	sed -i 's/RNAProxyJNI.RNA_SetProgress(swigCPtr, this, TProgressDialog.getCPtr(Progress), Progress/RNAProxyJNI.RNA_SetProgress(swigCPtr, this, getCPtrAndAddReference(Progress), Progress/g' RNA.java
	${SWIG_COMMAND} SWIG/TwoRNA.i
	${SWIG_COMMAND} SWIG/HybridRNA.i
	sed -i 's/HybridRNAProxyJNI.SWIGHybridRNAUpcast(cPtr)/cPtr/g' HybridRNA.java
	${SWIG_COMMAND} SWIG/Postscript_Wrapper.i
	${SWIG_COMMAND} SWIG/DotPlotHandler.i
	rm *Proxy.java
	mv *.java source/

# Remove unneeded object files and any executables
clean:
	rm -f *.o
	rm -f ../src/*.o
	rm -f ../RNA_class/*.o

# Build C++ code into object files and shared library
buildC:
	${SWIG_OBJECT_COMMAND} DotPlotHandler_wrap.o SWIG/DotPlotHandler_wrap.cxx
	${SWIG_OBJECT_COMMAND} DotPlotHandler.o ../src/DotPlotHandler.cpp
	${SWIG_OBJECT_COMMAND} Postscript_Wrapper_wrap.o SWIG/Postscript_Wrapper_wrap.cxx
	${SWIG_OBJECT_COMMAND} Dynalign_object_wrap.o SWIG/Dynalign_object_wrap.cxx
	${SWIG_OBJECT_COMMAND} Oligowalk_object_wrap.o SWIG/Oligowalk_object_wrap.cxx
	${SWIG_OBJECT_COMMAND} ProgressMonitor_wrap.o SWIG/ProgressMonitor_wrap.cxx
	${SWIG_OBJECT_COMMAND} TProgressDialog_wrap.o SWIG/TProgressDialog_wrap.cxx
	${SWIG_OBJECT_COMMAND} Thermodynamics_wrap.o SWIG/Thermodynamics_wrap.cxx
	${SWIG_OBJECT_COMMAND} RNA_wrap.o SWIG/RNA_wrap.cxx
	${SWIG_OBJECT_COMMAND} TwoRNA_wrap.o SWIG/TwoRNA_wrap.cxx
	${SWIG_OBJECT_COMMAND} HybridRNA_wrap.o SWIG/HybridRNA_wrap.cxx
	cd ../RNA_class; make HybridRNA_Library Dynalign_Library Oligo_Library JAVA=yes;
	cd ../src/PostscriptWriting; make Postscript_Library;
	${CXX} ${GUI_LFLAGS} -o ../exe/SharedLibrary.exe ../exe/*_GUI* \
	../exe/${POSTSCRIPT_LIBRARY} *.o
	mv ../exe/SharedLibrary.exe ../exe/${RNASTRUCTURE_LIBRARY}
	rm *.o SWIG/ProgressMonitor.o SWIG/TProgressDialog.o

# Compile Java code
compileJava:
	${JXX} ${JXXFLAGS} source/*.java

# Make the JAR file of GUI components
jar:
	cp -r ../manual/GUI/ docs/
	jar cvmf Manifest.txt RNAstructure.jar source/*.class sounds/*.wav docs/html/*.html docs/images/*.bmp images/*.*
	rm -rf docs/
	rm -f source/*.class
	mv RNAstructure.jar ../exe/RNAstructure.jar
