bin_PROGRAMS = seqpred

seqpred_SOURCES = \
	SeqPred.cpp \
	Predictor.cpp \
	DnaModel.cpp \
	ProteinModel.cpp \
	Model.cpp \
	Utils.cpp
seqpred_CFLAGS = -g -O3 -Wall $(SIMD_FLAGS)
seqpred_CPPFLAGS = -g -O3 -Wall $(SIMD_FLAGS)
seqpred_LDFLAGS = -g -O3 -Wall

pkginclude_HEADERS = Predictor.h DnaModel.h ProteinModel.h Model.h Utils.h
