#Your path to...
ARPACK=/usr/local/Cellar/arpack/3.6.2

ARPACK_LIB= -L${ARPACK}/lib/ -larpack
#ARPACK_FLAGS= -DUSE_ARPACK

#============================================================

TARGET = u1

SOURCES = main.cpp ran2s.cpp rangs.cpp
OBJS    = main.o ran2s.o rangs.o
INCLUDES = ran2s.h arpack_interface.h dOpHelpers.h

CXX=g++
CXXFLAGS = -O3 -g -std=c++11 ${ARPACK_FLAGS}

#============================================================

all: $(TARGET)

u1: ${OBJS}
	$(CXX) $(CXXFLAGS) -o u1 $(OBJS) ${ARPACK_LIB}

main.o: main.cpp ${INCLUDES} Makefile
	${CXX} ${CXXFLAGS} -c main.cpp

#============================================================

ALL_SOURCES = Makefile $(SOURCES) $(INCLUDES)

clean:
	rm -f $(TARGET) $(OBJS) core*

tar: $(ALL_SOURCES) 
	tar cvf $(TARGET).tar $(ALL_SOURCES) 

$(TARGET).ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)
