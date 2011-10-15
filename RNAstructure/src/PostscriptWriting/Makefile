###############################################################################
## Include external variables.
## Note that there is an external dependency not included here; this library
## depends on the Dynalign library.
## However, since this library is always compiled into an archive, compiling
## without the Dynalign library will not cause errors. Linking without the
## Dynalign library, however, will cause linking errors.
###############################################################################

include ../../compiler.h
include ../../library_defines.h

###############################################################################
## Define File Lists Within the Postscript Library.
###############################################################################

POSTSCRIPT_FILES = \
	Postscript_Object.o \
	Postscript_Piece.o \
	Postscript_Annotation_Handler.o \
	Postscript_Legend.o \
	Postscript_Probability_Legend.o \
	Postscript_SHAPE_Legend.o \
	Postscript_RGB_Legend.o \
	Postscript_Energy_Legend.o \
	Postscript_Log10_Legend.o \
	Postscript_Image.o \
	Postscript_Structure.o \
	Postscript_Circle_Structure.o \
	Postscript_Circle_Structures.o \
	Postscript_Radial_Structure.o \
	Postscript_Dot_Plot.o \
	Postscript_Single_Plot.o \
	Postscript_Energy_Plot.o \
	Postscript_Probability_Plot.o \
	Postscript_Dynalign_Plot.o

LIBRARY_FILES = ${POSTSCRIPT_FILES} Postscript_Wrapper.o

###############################################################################
## Define Actions
###############################################################################

all: Postscript_Library
Postscript_Library: ${LIBRARY_FILES}
	ar rc ../../exe/${POSTSCRIPT_LIBRARY} ${LIBRARY_FILES}

###############################################################################
## Rules for Postscript Object Files (POSTSCRIPT_FILES)
###############################################################################

# General Postscript objects: general object, general large image, and general
# small piece
Postscript_Object.o: \
	Postscript_Object.cpp Postscript_Object.h

Postscript_Image.o: \
	Postscript_Image.cpp Postscript_Image.h \
	Postscript_Object.cpp Postscript_Object.h

Postscript_Piece.o: \
	Postscript_Piece.cpp Postscript_Piece.h \
	Postscript_Object.cpp Postscript_Object.h

# Postscript annotation handler
Postscript_Annotation_Handler.o: \
	Postscript_Annotation_Handler.cpp Postscript_Annotation_Handler.h \
	Postscript_Piece.cpp Postscript_Piece.h \
	Postscript_Object.cpp Postscript_Object.cpp

# Postscript legends
Postscript_Legend.o: \
	Postscript_Legend.cpp Postscript_Legend.h \
	Postscript_Piece.cpp Postscript_Piece.h \
	Postscript_Object.cpp Postscript_Object.h

Postscript_Probability_Legend.o: \
	Postscript_Probability_Legend.cpp Postscript_Probability_Legend.h \
	Postscript_Legend.cpp Postscript_Legend.h \
	Postscript_Piece.cpp Postscript_Piece.h \
	Postscript_Object.cpp Postscript_Object.h

Postscript_SHAPE_Legend.o: \
	Postscript_SHAPE_Legend.cpp Postscript_SHAPE_Legend.h \
	Postscript_Legend.cpp Postscript_Legend.h \
	Postscript_Piece.cpp Postscript_Piece.h \
	Postscript_Object.cpp Postscript_Object.h

Postscript_RGB_Legend.o: \
	Postscript_RGB_Legend.cpp Postscript_RGB_Legend.h \
	Postscript_Legend.cpp Postscript_Legend.h \
	Postscript_Piece.cpp Postscript_Piece.h \
	Postscript_Object.cpp Postscript_Object.h

Postscript_Energy_Legend.o: \
	Postscript_Energy_Legend.cpp Postscript_Energy_Legend.h \
	Postscript_RGB_Legend.cpp Postscript_RGB_Legend.h \
	Postscript_Legend.cpp Postscript_Legend.h \
	Postscript_Piece.cpp Postscript_Piece.h \
	Postscript_Object.cpp Postscript_Object.h

Postscript_Log10_Legend.o: \
	Postscript_Log10_Legend.cpp Postscript_Log10_Legend.h \
	Postscript_RGB_Legend.cpp Postscript_RGB_Legend.h \
	Postscript_Legend.cpp Postscript_Legend.h \
	Postscript_Piece.cpp Postscript_Piece.h \
	Postscript_Object.cpp Postscript_Object.h

# Postscript structures
Postscript_Structure.o: \
	Postscript_Structure.cpp Postscript_Structure.h \
	Postscript_Image.cpp Postscript_Image.h \
	Postscript_Object.cpp Postscript_Object.h

Postscript_Circle_Structure.o: \
	Postscript_Circle_Structure.cpp Postscript_Circle_Structure.h \
	Postscript_Structure.cpp Postscript_Structure.h \
	Postscript_Image.cpp Postscript_Image.h \
	Postscript_Object.cpp Postscript_Object.h

Postscript_Circle_Structures.o: \
	Postscript_Circle_Structures.cpp Postscript_Circle_Structures.h \
	Postscript_Circle_Structure.cpp Postscript_Circle_Structure.h \
	Postscript_Structure.cpp Postscript_Structure.h \
	Postscript_Image.cpp Postscript_Image.h \
	Postscript_Object.cpp Postscript_Object.h

Postscript_Radial_Structure.o: \
	Postscript_Radial_Structure.cpp Postscript_Radial_Structure.h \
	Postscript_Structure.cpp Postscript_Structure.h \
	Postscript_Image.cpp Postscript_Image.h \
	Postscript_Object.cpp Postscript_Object.h

# Postscript dot plots
Postscript_Dot_Plot.o: \
	Postscript_Dot_Plot.cpp Postscript_Dot_Plot.h \
	Postscript_Image.cpp Postscript_Image.h \
	Postscript_Object.cpp Postscript_Object.h

Postscript_Single_Plot.o: \
	Postscript_Single_Plot.cpp Postscript_Single_Plot.h \
	Postscript_Dot_Plot.cpp Postscript_Dot_Plot.h \
	Postscript_Image.cpp Postscript_Image.h \
	Postscript_Object.cpp Postscript_Object.h

Postscript_Energy_Plot.o: \
	Postscript_Energy_Plot.cpp Postscript_Energy_Plot.h \
	Postscript_Single_Plot.cpp Postscript_Single_Plot.h \
	Postscript_Dot_Plot.cpp Postscript_Dot_Plot.h \
	Postscript_Image.cpp Postscript_Image.h \
	Postscript_Object.cpp Postscript_Object.h

Postscript_Probability_Plot.o: \
	Postscript_Probability_Plot.cpp Postscript_Probability_Plot.h \
	Postscript_Single_Plot.cpp Postscript_Single_Plot.h \
	Postscript_Dot_Plot.cpp Postscript_Dot_Plot.h \
	Postscript_Image.cpp Postscript_Image.h \
	Postscript_Object.cpp Postscript_Object.h

Postscript_Dynalign_Plot.o: \
	Postscript_Dynalign_Plot.cpp Postscript_Dynalign_Plot.h \
	Postscript_Dot_Plot.cpp Postscript_Dot_Plot.h \
	Postscript_Image.cpp Postscript_Image.h \
	Postscript_Object.cpp Postscript_Object.h

###############################################################################
## Rule for Postscript utility wrapper -- basically, always recompile it.
###############################################################################

Postscript_Wrapper.o: \
	${POSTSCRIPT_FILES}