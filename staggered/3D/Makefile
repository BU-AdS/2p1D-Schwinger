#Your path to...
ARPACK=/usr/local/Cellar/arpack/3.6.2

ARPACK_LIB= -L${ARPACK}/lib/ -larpack
ARPACK_FLAGS= -DUSE_ARPACK

#============================================================

TARGET = u1

SOURCES = main.cpp
OBJS    = main.o
INC_PATH=-I/Users/deanhowarth/2p1D-Schwinger/include

CXX=g++
CXXFLAGS = -O3 -g -std=c++11 ${INC_PATH} ${ARPACK_FLAGS} 

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
