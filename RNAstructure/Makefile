###############################################################################
## Included definitions of outside variables.
###############################################################################

include compiler.h
include library_defines.h

###############################################################################
## General Build Actions
###############################################################################

# Echo message on command usage to build repository.
instructions:
	@echo 'Use "make all" to create all executables.'
	@echo 'Use "make AllSub" to create executable "AllSub."'
	@echo 'Use "make bifold" to create executable "bifold."'
	@echo 'Use "make bipartition" to create executable "bipartition."'
	@echo 'Use "make CircleCompare" to create executable "CircleCompare."'
	@echo 'Use "make ct2dot" to create executable "ct2dot."'
	@echo 'Use "make dot2ct" to create executable "dot2ct."'
	@echo 'Use "make draw" to create executable "draw."'
	@echo 'Use "make dynalign" to create exectuable "dynalign."'
	@echo 'Use "make dynalign-smp" to create executable "dynalign-smp."'
	@echo 'Use "make DuplexFold" to create executable "DuplexFold."'
	@echo 'Use "make DynalignDotPlot" to create executable' \
		'"DynalignDotPlot."'
	@echo 'Use "make EnergyPlot" to create executable "EnergyPlot."'
	@echo 'Use "make efn2" to create executable "efn2."'
	@echo 'Use "make Fold" to create executable "Fold."'
	@echo 'Use "make MaxExpect" to create executable "MaxExpect."'
	@echo 'Use "make NAPSS" to create the executable "NAPSS."'
	@echo 'Use "make oligoscreen" to create executable "oligoscreen."'
	@echo 'Use "make OligoWalk" to create executable "OligoWalk."'
	@echo 'Use "make partition" to create executable "partition."'
	@echo 'Use "make PARTS" to create executable "PARTS."'
	@echo 'Use "make ProbabilityPlot" to create executable' \
		'"ProbabilityPlot."'
	@echo 'Use "make ProbablePair" to create executable "ProbablePair."'
	@echo 'Use "make ProbKnot" to create executable "ProbKnot."'
	@echo 'Use "make refold" to create executable "refold."'
	@echo 'Use "make RemovePseudoknots" to create executable' \
		'"RemovePseudoknots."'
	@echo 'Use "make scorer" to create executable "scorer."'
	@echo 'Use "make stochastic" to create executable "stochastic."'

# Make all the executables at once.
all:
	@echo "Building of entire RNAstructure repository started."
	@echo
	make AllSub;
	make bifold;
	make bipartition;
	make CircleCompare;
	make ct2dot;
	make dot2ct;
	make draw;
	make DuplexFold;
	make dynalign;
	make dynalign-smp;
	make DynalignDotPlot;
	make efn2;
	make EnergyPlot;
	make Fold;
	make MaxExpect;
	make NAPSS;
	make oligoscreen;
	make OligoWalk;
	make partition;
	make PARTS;
	make ProbabilityPlot;
	make ProbablePair;
	make ProbKnot;
	make refold;
	make RemovePseudoknots;
	make scorer;
	make stochastic;
	@echo
	@echo "Building of entire RNAstructure repository finished."

# Copy the executables to the /usr/local directory.
install:
	cp -r exe/ /usr/local/RNAstructure

###############################################################################
## Rules and objects used in all text interfaces.
###############################################################################

# Common files list.
COMMON_FILES = \
	src/ParseCommandLine.o

# Unix-like command line parser.
src/ParseCommandLine.o: src/ParseCommandLine.cpp src/ParseCommandLine.h

###############################################################################
## Rules and objects used in some groups of interfaces, but not others.
###############################################################################

# Used in all interfaces where scoring is done.
src/score.o: src/score.cpp src/score.h

# Used in dot plot interfaces.
DOT_PLOT_GROUP = \
	DotPlots/DotPlots.o \
	${POSTSCRIPT_GROUP}

DotPlots/DotPlots.o: DotPlots/DotPlots.cpp DotPlots/DotPlots.h
src/DotPlotHandler.o: src/DotPlotHandler.cpp src/DotPlotHandler.h

###############################################################################
## Dependency commands and variables for libraries.
###############################################################################

# RNA library data.
exe/${RNA_LIBRARY}:
	cd RNA_class; make RNA_Library;
RNA_GROUP = ${COMMON_FILES} exe/${RNA_LIBRARY}

# HybridRNA library data.
exe/${HYBRID_RNA_LIBRARY}:
	cd RNA_class; make HybridRNA_Library;
HYBRID_RNA_GROUP = ${COMMON_FILES} exe/${HYBRID_RNA_LIBRARY}

# Dynalign library data.
exe/${DYNALIGN_LIBRARY}:
	cd RNA_class; make Dynalign_Library;
exe/${DYNALIGN_SMP_LIBRARY}:
	cd RNA_class; make DynalignSMP_Library;

# Oligo library data.
exe/${OLIGO_LIBRARY}:
	cd RNA_class; make Oligo_Library;
OLIGO_GROUP = ${COMMON_FILES} exe/${OLIGO_LIBRARY}

# Postscript library data.
exe/${POSTSCRIPT_LIBRARY}:
	cd src/PostscriptWriting; make Postscript_Library;
POSTSCRIPT_GROUP = ${COMMON_FILES} src/DotPlotHandler.o \
	exe/${POSTSCRIPT_LIBRARY} exe/${DYNALIGN_LIBRARY}

###############################################################################
## Commands and dependencies specific to the AllSub interface.
###############################################################################

# AllSub executable.
AllSub: exe/AllSub
exe/AllSub: AllSub/AllSub.o ${RNA_GROUP}
	${LINK} AllSub/AllSub.o ${RNA_GROUP}

# AllSub object.
AllSub/AllSub.o: AllSub/AllSub.cpp AllSub/AllSub.h

###############################################################################
## Commands and dependencies specific to the bifold interface.
###############################################################################

# bifold executable.
bifold: exe/bifold
exe/bifold: bifold/bifold.o ${HYBRID_RNA_GROUP}
	${LINK} bifold/bifold.o ${HYBRID_RNA_GROUP}

# bifold object.
bifold/bifold.o: bifold/bifold.cpp bifold/bifold.h

###############################################################################
## Commands and dependencies specific to the bipartition interface.
###############################################################################

# bipartition executable.
bipartition: exe/bipartition
exe/bipartition: bipartition/bipartition.o ${HYBRID_RNA_GROUP}
	${LINK} bipartition/bipartition.o ${HYBRID_RNA_GROUP}

# bipartition object.
bipartition/bipartition.o: \
	bipartition/bipartition.cpp bipartition/bipartition.h

###############################################################################
## Commands and dependencies specific to the CircleCompare interface.
###############################################################################

# Define CircleCompare file group.
CIRCLECOMPARE_GROUP = \
	CircleCompare/CircleCompare.o \
	CircleCompare/CircleCompare_Interface.o \
	src/score.o \
	${POSTSCRIPT_GROUP}

# CircleCompare executable.
CircleCompare: exe/CircleCompare
exe/CircleCompare: ${CIRCLECOMPARE_GROUP}
	${LINK} ${CIRCLECOMPARE_GROUP}

# Rules for unique files in CircleCompare file group.
CircleCompare/CircleCompare.o: \
	CircleCompare/CircleCompare.cpp CircleCompare/CircleCompare.h

CircleCompare/CircleCompare_Interface.o: \
	CircleCompare/CircleCompare_Interface.cpp \
	CircleCompare/CircleCompare_Interface.h \
	CircleCompare/CircleCompare.cpp \
	CircleCompare/CircleCompare.h

###############################################################################
## Commands and dependencies specific to the ct2dot interface.
###############################################################################

# ct2dot executable.
ct2dot: exe/ct2dot
exe/ct2dot: ct2dot/ct2dot.o ${RNA_GROUP}
	${LINK} ct2dot/ct2dot.o ${RNA_GROUP}

# ct2dot object.
ct2dot/ct2dot.o: ct2dot/ct2dot.cpp ct2dot/ct2dot.h

###############################################################################
## Commands and dependencies specific to the dot2ct interface.
###############################################################################

# dot2ct executable.
dot2ct: exe/dot2ct
exe/dot2ct: dot2ct/dot2ct.o ${RNA_GROUP}
	${LINK} dot2ct/dot2ct.o ${RNA_GROUP}

# dot2ct object.
dot2ct/dot2ct.o: dot2ct/dot2ct.cpp dot2ct/dot2ct.h

###############################################################################
## Commands and dependencies specific to the draw interface.
###############################################################################

# draw executable.
draw: exe/draw
exe/draw: draw/DrawStructure.o ${POSTSCRIPT_GROUP}
	${LINK} draw/DrawStructure.o ${POSTSCRIPT_GROUP}

# draw object.
draw/DrawStructure.o: draw/DrawStructure.cpp draw/DrawStructure.h

###############################################################################
## Commands and dependencies specific to the DuplexFold interface.
###############################################################################

# DuplexFold executable.
DuplexFold: exe/DuplexFold
exe/DuplexFold: DuplexFold/DuplexFold.o ${HYBRID_RNA_GROUP}
	${LINK} DuplexFold/DuplexFold.o ${HYBRID_RNA_GROUP}

# DuplexFold object.
DuplexFold/DuplexFold.o: DuplexFold/DuplexFold.cpp DuplexFold/DuplexFold.h

###############################################################################
## Commands and dependencies specific to the dynalign interfaces.
## This section contains targets for both dynalign and dynalign-smp.
###############################################################################

dynalign: exe/dynalign
exe/dynalign: dynalign/dynaligninterface.o exe/${DYNALIGN_LIBRARY}
	${LINK} dynalign/dynaligninterface.o exe/${DYNALIGN_LIBRARY}

dynalign-smp: exe/dynalign-smp
exe/dynalign-smp: dynalign/dynaligninterface-smp.o exe/${DYNALIGN_SMP_LIBRARY}
	${LINK} dynalign/dynaligninterface-smp.o exe/${DYNALIGN_SMP_LIBRARY} -lpthread

dynalign/dynaligninterface.o: dynalign/dynaligninterface.cpp
dynalign/dynaligninterface-smp.o: dynalign/dynaligninterface.cpp
	${CXX} ${CXXFLAGS} -DDYNALIGN_SMP -c -o $@ \
	dynalign/dynaligninterface.cpp

###############################################################################
## Commands and dependencies specific to the EnergyPlot interface.
###############################################################################

# EnergyPlot executable.
DynalignDotPlot: exe/DynalignDotPlot
exe/DynalignDotPlot: ${DOT_PLOT_GROUP}
	${LINK} ${DOT_PLOT_GROUP}

###############################################################################
## Commands and dependencies specific to the efn2 interface.
###############################################################################

# efn2 executable.
efn2: exe/efn2
exe/efn2: efn2/efn2.o ${RNA_GROUP}
	${LINK} efn2/efn2.o ${RNA_GROUP}

# efn2 object.
efn2/efn2.o: efn2/efn2.cpp efn2/efn2.h

###############################################################################
## Commands and dependencies specific to the EnergyPlot interface.
###############################################################################

# EnergyPlot executable.
EnergyPlot: exe/EnergyPlot
exe/EnergyPlot: ${DOT_PLOT_GROUP}
	${LINK} ${DOT_PLOT_GROUP}

###############################################################################
## Commands and dependencies specific to the Fold inteerface.
###############################################################################

# Fold executable.
Fold: exe/Fold
exe/Fold: fold/Fold.o ${RNA_GROUP}
	${LINK} fold/Fold.o ${RNA_GROUP}

# Fold object.
fold/Fold.o: fold/Fold.cpp fold/Fold.h

###############################################################################
## Commands and dependencies specific to the MaxExpect interface.
###############################################################################

# MaxExpect executable.
MaxExpect: exe/MaxExpect
exe/MaxExpect: MaxExpect/MaxExpectInterface.o ${RNA_GROUP}
	${LINK} MaxExpect/MaxExpectInterface.o ${RNA_GROUP}

# MaxExpect object.
MaxExpect/MaxExpectInterface.o: \
	MaxExpect/MaxExpectInterface.cpp MaxExpect/MaxExpect.h

###############################################################################
## Commands and dependencies specific to the NAPSS interface.
###############################################################################

# NAPSS executable.
NAPSS: exe/NAPSS
exe/NAPSS:
	cd napss; make
	mv napss/NAPSS exe/

###############################################################################
## Commands and dependencies specific to the oligoscreen interface.
###############################################################################

# oligoscreen executable.
oligoscreen: exe/oligoscreen
exe/oligoscreen: oligoscreen/oligoscreen.o ${OLIGO_GROUP}
	${LINK} oligoscreen/oligoscreen.o ${OLIGO_GROUP}

# oligoscreen object.
oligoscreen/oligoscreen.o: \
	oligoscreen/oligoscreen.cpp oligoscreen/oligoscreen.h

###############################################################################
## Commands and dependencies specific to the OligoWalk interface.
###############################################################################

# OligoWalk executable.
OligoWalk:
	cd oligowalk/src; make OligoWalk
	mv oligowalk/exe/OligoWalk exe/

###############################################################################
## Commands and dependencies specific to the partition interface.
###############################################################################

# partition executable.
partition: exe/partition
exe/partition: pfunction/partition.o ${RNA_GROUP}
	${LINK} pfunction/partition.o ${RNA_GROUP}

# partition object.
pfunction/partition.o: pfunction/partition.cpp pfunction/partition.h

###############################################################################
## Commands and dependencies specific to the PARTS interface.
###############################################################################

# Build PARTS using the PARTS Makefile.
PARTS:  exe/PARTS
exe/PARTS:
	cd PARTS; make PARTS

###############################################################################
## Commands and dependencies specific to the ProbKnot interface.
###############################################################################

# ProbKnot executable.
ProbKnot: exe/ProbKnot
exe/ProbKnot: ProbKnot/ProbKnot_Interface.o ${RNA_GROUP}
	${LINK} ProbKnot/ProbKnot_Interface.o ${RNA_GROUP}

# ProbKnot object.
ProbKnot/ProbKnot_Interface.o: \
	ProbKnot/ProbKnot_Interface.cpp ProbKnot/ProbKnot_Interface.h

###############################################################################
## Commands and dependencies specific to the ProbabilityPlot interface.
###############################################################################

# ProbabilityPlot executable.
ProbabilityPlot: exe/ProbabilityPlot
exe/ProbabilityPlot: ${DOT_PLOT_GROUP}
	${LINK} ${DOT_PLOT_GROUP}

###############################################################################
## Commands and dependencies specific to the ProbablePair interface.
###############################################################################

# ProbablePair executable.
ProbablePair: exe/ProbablePair
exe/ProbablePair: ProbablePair/ProbablePair.o ${RNA_GROUP}
	${LINK} ProbablePair/ProbablePair.o ${RNA_GROUP}

# ProbablePair object.
ProbablePair/ProbablePair.o: \
	ProbablePair/ProbablePair.cpp ProbablePair/ProbablePair.h

###############################################################################
## Commands and dependencies specific to the RemovePseudoknots interface.
###############################################################################

# RemovePseudoknots executable.
RemovePseudoknots: exe/RemovePseudoknots
exe/RemovePseudoknots: RemovePseudoknots/RemovePseudoknots.o ${RNA_GROUP}
	${LINK} RemovePseudoknots/RemovePseudoknots.o ${RNA_GROUP}

# RemovePseudoknots object.
RemovePseudoknots/RemovePseudoknots.o: \
	RemovePseudoknots/RemovePseudoknots.cpp \
	RemovePseudoknots/RemovePseudoknots.h

###############################################################################
## Commands and dependencies specific to the refold interface.
###############################################################################

# refold executable.
refold: exe/refold
exe/refold: refold/refold.o ${RNA_GROUP}
	${LINK} refold/refold.o ${RNA_GROUP}

# refold object.
refold/refold.o: refold/refold.cpp refold/refold.h

###############################################################################
## Commands and dependencies specific to the scorer interface.
###############################################################################

# scorer executable.
scorer: exe/scorer
exe/scorer: scorer/Scorer_Interface.o src/score.o ${RNA_GROUP}
	${LINK} scorer/Scorer_Interface.o src/score.o ${RNA_GROUP}

# scorer object.
scorer/Scorer_Interface.o: \
	scorer/Scorer_Interface.cpp scorer/Scorer_Interface.h

###############################################################################
## Commands and dependencies specific to the stochastic interface.
###############################################################################

# stochastic executable.
stochastic: exe/stochastic
exe/stochastic: stochastic/stochastic.o ${RNA_GROUP}
	${LINK} stochastic/stochastic.o ${RNA_GROUP}

# stochastic object.
stochastic/stochastic.o: stochastic/stochastic.cpp stochastic/stochastic.h

###############################################################################
## Object and Executable Cleanup
###############################################################################

# Remove object files.
clean:
	find . -depth -name '*.o' -delete

# Remove object files and executables.
realclean: clean
	find exe -maxdepth 1 -type f -perm ${PERM} \
	! -name RNAstructureScript \
	! -name RNAstructure.bat \
	-delete
	find exe -name *${LIBEXT} -delete
	rm -f exe/RNAstructure.jar