bin_PROGRAMS = foreseqs

foreseqs_SOURCES = \
	Predictor.cpp \
	DnaModel.cpp \
	ProteinModel.cpp \
	Model.cpp \
	Utils.cpp \
	ForeSeqs.cpp
foreseqs_CFLAGS = -g -O3 -Wall $(SIMD_FLAGS)
foreseqs_CPPFLAGS = -g -O3 -Wall $(SIMD_FLAGS)
foreseqs_LDFLAGS = -g -O3 -Wall

pkginclude_HEADERS = Predictor.h DnaModel.h ProteinModel.h Model.h Utils.h
