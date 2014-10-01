CC=g++
CFLAGS=-c -Wall -Ilibpll/include
LDFLAGS=-Llibpll/lib
LIBS=-lpll-sse3 -lm

SOURCES=src/SeqPred.cpp src/Predictor.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=seqpred

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

oneline:
	g++ -g -O2 -Ilibpll/include -Llibpll/lib -o seqpred src/SeqPred.cpp src/Predictor.cpp -lpll-sse3 -lm
