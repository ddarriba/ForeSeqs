bin_PROGRAMS = foreseqs

foreseqs_SOURCES = ForeSeqs.cpp \
									 Utils.cpp \
									 Alignment.cpp \
									 Predictor.cpp \
									 DnaModel.cpp \
									 ProteinModel.cpp \
								   Model.cpp \
									 Phylo.cpp


foreseqs_CFLAGS = -g -O3 -Wall $(SIMD_FLAGS)
foreseqs_CPPFLAGS = -g -O3 -Wall $(SIMD_FLAGS)
foreseqs_LDFLAGS = -g -O3 -Wall
foreseqs_LDADD = -l:libpll_algorithm.a -l:libpll_binary.a -l:libpll_optimize.a -l:libpll_msa.a -l:libpll_tree.a -l:libpll_util.a -l:libpll.a

pkginclude_HEADERS = Predictor.h DnaModel.h ProteinModel.h Model.h Utils.h
