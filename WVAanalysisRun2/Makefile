CXXFLAGS=-Wall # -pedantic -ansi
ROOT_LIB:=`root-config --libs --glibs`
ROOT_FLAGS:=`root-config --cflags --ldflags`
ROOT_INCLUDE:=`root-config --incdir`

DEPS= interface/setOutputTree.h interface/METzCalculator.h interface/analysisUtils.h interface/setInputTree.h interface/METzCalculator_Run2.h interface/pseudodataNtuple.h interface/SF.h
DEPS_OBJ= lib/setOutputTree.o lib/METzCalculator.o lib/analysisUtils.o lib/setInputTree.o lib/METzCalculator_Run2.o lib/pseudodataNtuple.o lib/SF.o

CC = g++
CFLAGS = -Wall

lib/%.o: src/%.cc $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< $(ROOT_LIB) $(ROOT_FLAGS)

all: produceWWNtuples.exe

produceWWNtuples.exe: bin/produceWWNtuples.cpp $(DEPS_OBJ)
	g++ $(CFLAGS) -o $@ $^ $(ROOT_LIB) $(ROOT_FLAGS)

clean:
	rm -f lib/*.o
	rm -f lib/*.d
	rm -f *.exe
