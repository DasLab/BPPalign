include ../compiler.h

COMMON_OBJ_FILES=../src/algorithm.o ../src/arrayclass.o \
	../src/configfile.o ../src/dotarray.o ../src/forceclass.o \
	../src/rna_library.o ../src/stackclass.o \
	../src/stackstruct.o \
	../src/structure.o ../src/TProgressDialog.o 

SEQ_OBJ_FILES=napss.o

NAPSS:	${COMMON_OBJ_FILES} ${SEQ_OBJ_FILES} 
	${LINK} ${COMMON_OBJ_FILES} ${SEQ_OBJ_FILES}

# 	Sequential Object Files (SEQ_OBJ_FILES)

napss.o: napss.cpp ../src/algorithm.h \
	../src/arrayclass.h ../src/configfile.h ../src/defines.h ../src/dotarray.h \
	../src/forceclass.h \
	../src/platform.h ../src/rna_library.h ../src/stackclass.h \
	../src/stackstruct.h \
	../src/structure.h ../src/TProgressDialog.h

# 	Common Object Files (COMMON_OBJ_FILES)

algorithm.o: ../src/algorithm.cpp ../src/algorithm.h ../src/arrayclass.h \ 
	../src/defines.h	\
	../src/dotarray.h ../src/forceclass.h ../src/platform.h \
	../src/rna_library.h \
	../src/stackclass.h ../src/stackstruct.h ../src/structure.h \ 
	../src/TProgressDialog.h

../src/arrayclass.o: ../src/arrayclass.cpp ../src/arrayclass.h ../src/defines.h

../src/configfile.o: ../src/configfile.h

../src/dotarray.o: ../src/dotarray.cpp ..//src/dotarray.h ../src/defines.h

../src/forceclass.o: ../src/forceclass.cpp ../src/forceclass.h

#../src/platform.o: ../src/platform.cpp ../src/platform.h

rna_library.o: ../src/crna_library.cpp ../src/rna_library.h ../src/defines.h \
	../src/platform.h ../src/structure.h \
	../src/stackclass.h ../src/stackstruct.h

../src/stackclass.o: ../src/stackclass.cpp ../src/stackclass.h ../src/defines.h

../src/stackstruct.o: ../src/stackstruct.cpp ../src/stackstruct.h

../src/structure.o: ../src/structure.cpp ../src/structure.h ../src/defines.h \
	../src/platform.h

../src/TProgressDialog.o: ../src/TProgressDialog.cpp ../src/TProgressDialog.h

clean:
	rm -f *.o
	rm -f ../src/*.o

