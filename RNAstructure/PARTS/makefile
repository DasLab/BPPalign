###############################################################################
## Include external rules and variable definitions.
###############################################################################

include ../compiler.h
include ../library_defines.h

###############################################################################
## Variables which define working directories.
###############################################################################

LIB_DIR = src/lib
PHMM_LIB_DIR = src/lib/phmm
PARTS_DIR = src/parts
PARTS_MAIN = src/parts/main.o

###############################################################################
## Rules for object groups.
###############################################################################

# Utility objects
util_objs = \
	${LIB_DIR}/utils/file/utils.o \
	${LIB_DIR}/structure/structure.o \
	${LIB_DIR}/utils/xmath/log/xlog_math.o \
	${LIB_DIR}/utils/xmath/linear/linear_math.o \
	${LIB_DIR}/utils/rng/rng.o \
	${LIB_DIR}/structure/folding_constraints.o

# phmm library objects
phmm_lib_objs = \
	${PHMM_LIB_DIR}/aln_env_utils.o \
	${PHMM_LIB_DIR}/p_alignment.o \
	${PHMM_LIB_DIR}/phmm_aln.o \
	${PHMM_LIB_DIR}/phmm_array.o \
	${PHMM_LIB_DIR}/phmm.o \
	${PHMM_LIB_DIR}/phmm_ml_loops.o \
	${PHMM_LIB_DIR}/phmm_pp_loops.o

# PARTS objects
parts_objs = \
	${PARTS_DIR}/ppf_ss.o \
	${PARTS_DIR}/ppf_v_mhe.o \
	${PARTS_DIR}/ppf_w_mhi.o \
	${PARTS_DIR}/ppf_w_l.o \
	${PARTS_DIR}/ppf_w.o \
	${PARTS_DIR}/ppf_w_mbl.o \
	${PARTS_DIR}/ppf_w_mb.o \
	${PARTS_DIR}/stoch_tb/stoch_sampled_structures.o \
	${PARTS_DIR}/stoch_tb/stoch_sampled_alignment.o \
	${PARTS_DIR}/stoch_tb/stoch_sampled_str_aln_sample_set.o \
	${PARTS_DIR}/stoch_tb/stoch_sampling_math.o \
	${PARTS_DIR}/ppf_w_ext.o \
	${PARTS_DIR}/ppf_operators.o \
	${PARTS_DIR}/process_sequences.o \
	${PARTS_DIR}/ppf_loops.o \
	${PARTS_DIR}/single_pf_array.o \
	${PARTS_DIR}/alignment_priors.o \
	${PARTS_DIR}/pf_alignment.o \
	${PARTS_DIR}/template_pf_array.o \
	${PARTS_DIR}/ppf_tb_stack.o \
	${PARTS_DIR}/map_alignment.o \
	${PARTS_DIR}/map_structures.o \
	${PARTS_DIR}/ppf_cli.o \
	${PARTS_DIR}/ppf_scale.o \
	${PARTS_DIR}/ppf_timer.o \
	${PARTS_DIR}/ppf_progress_bar.o \
	${PARTS_DIR}/map_mhr_info.o \
	${PARTS_DIR}/map_results.o \
	${PARTS_DIR}/pp_results.o \
	${PARTS_DIR}/array_mem_manager.o \
	${PARTS_DIR}/array_file_manager.o \
	${LIB_DIR}/utils/ansi_string/ansi_string.o \
	${LIB_DIR}/utils/rng/seed_manager.o

###############################################################################
## PARTS Library and command line executable building rules.
###############################################################################

all: PARTS

PARTS: ${PARTS_MAIN}
	cd ../RNA_class; make RNA_Library;
	make PARTS_Library;
	${LINK} ${PARTS_MAIN} ../exe/${PARTS_LIBRARY} ../exe/${RNA_LIBRARY}
	mv PARTS ../exe/

PARTS_Library: ${parts_objs} ${phmm_lib_objs} ${util_objs}
	ar rc ${PARTS_LIBRARY} ${parts_objs} ${util_objs} ${phmm_lib_objs} 
	mv ${PARTS_LIBRARY} ../exe/

###############################################################################
## Object and executable cleanup.
###############################################################################

# Remove object files.
clean:
	find . -depth -name '*.o' -delete

# Remove object files and executables.
realclean: clean
	rm -f ../exe/PARTS
	rm -f ../exe/${PARTS_LIBRARY}
