TARGET = jk_error
OBJS   = jk_error.o
CXX = g++
CXXFLAGS = -O3 -g -std=c++11 
#============================================================

all: ${TARGET}

${TARGET}: ${OBJS}
	${CXX} ${CXXFLAGS} -o ${TARGET} ${OBJS}

jk_error.o: jk_error.cpp Makefile 
	${CXX} ${CXXFLAGS} -c jk_error.cpp

#============================================================

ALL_SOURCES = Makefile $(SOURCES) 

clean:
	rm -f $(TARGET) $(OBJS) core*
