OBJS = jkEffMass.o
CXX = g++
CXXFLAGS = -O3 -g -std=c++11 
#============================================================

all: jkEffMass

jkEffMass: jkEffMass.o
	${CXX} ${CXXFLAGS} -o jkEffMass jkEffMass.o

jkEffMass.o: jkEffMass.cpp Makefile 
	${CXX} ${CXXFLAGS} -c jkEffMass.cpp

#============================================================

all: creutzRatio

creutzRatio: creutzRatio.o
	${CXX} ${CXXFLAGS} -o creutzRatio creutzRatio.o

creutzRatio.o: creutzRatio.cpp Makefile 
	${CXX} ${CXXFLAGS} -c creutzRatio.cpp

#============================================================

ALL_SOURCES = Makefile creutzRatio.cpp jkEffMass.cpp

clean:
	rm -f jkEffMass jkEffMass.o creutzRatio creutzRatio.o core*
