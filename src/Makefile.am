bin_PROGRAMS = seqpred

seqpred_SOURCES = \
	SeqPred.cpp \
	Predictor.cpp \
	DnaModel.cpp \
	ProteinModel.cpp \
	Model.cpp \
	Utils.cpp
seqpred_CFLAGS = -g -O3 -Wall
seqpred_CPPFLAGS = -g -O3 -Wall
# seqpred_LDFLAGS = -Llibpll/lib -lpll-sse3-pthreads -lm -lpthread
seqpred_LDFLAGS = -g -O3 -Wall

#if HAVE_PTHREADS
#  seqpred_CFLAGS += -DHAVE_PTHREADS
#  seqpred_LDFLAGS += -DHAVE_PTHREADS
#endif

#if HAVE_AVX
#  seqpred_CFLAGS += -DHAVE_AVX
#  seqpred_LDFLAGS += -DHAVE_AVX
#endif

#if HAVE_SSE
#  seqpred_CFLAGS += -DHAVE_SSE
#  seqpred_LDFLAGS += -DHAVE_SSE
#endif

