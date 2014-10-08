CC=g++
CFLAGS=-g -Ilibpll/include -O4 -Wall
LDFLAGS=-g -Llibpll/lib -O4 -Wall
LIBS=-lpll-sse3 -lm

DBGFLAGS= # -DPRINT_ANCESTRAL -DPRINT_TRACE
SOURCES=src/SeqPred.cpp src/Predictor.cpp src/DnaModel.cpp src/ProteinModel.cpp src/Model.cpp src/Utils.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=seqpred

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

.cpp.o:
	$(CC) $(DBGFLAGS) -c $(CFLAGS) $< -o $@

clean:
	rm -f src/*.o $(EXECUTABLE)

oneline:
	g++ -g -O2 -Ilibpll/include -Llibpll/lib -o seqpred src/SeqPred.cpp src/Predictor.cpp -lpll-sse3 -lm
