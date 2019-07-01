#============================================================
CXX = g++
CXXFLAGS = -O3 -g -std=c++11 
#============================================================

all: effMass

jkEffMass: effMass.o
	${CXX} ${CXXFLAGS} -o effMass effMass.o

effMass.o: effMass.cpp Makefile 
	${CXX} ${CXXFLAGS} -c effMass.cpp

#============================================================

all: creutzRatio

creutzRatio: creutzRatio.o
	${CXX} ${CXXFLAGS} -o creutzRatio creutzRatio.o

creutzRatio.o: creutzRatio.cpp Makefile 
	${CXX} ${CXXFLAGS} -c creutzRatio.cpp

#============================================================

ALL_SOURCES = Makefile creutzRatio.cpp effMass.cpp

clean:
	rm -f effMass effMass.o creutzRatio creutzRatio.o core*
