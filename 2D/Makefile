.SUFFIXES:
.SUFFIXES: .o .cpp

QDPHOME = /home/brower/arch/qdp++/scalar
LIBNAME = qdp

LIBDIRS = -L/usr/local/Cellar/arpack/3.6.2/lib
#LIBDIRS = -L/usr/local/opt/libomp/lib
#INC = -I/usr/local/opt/libomp/include
INC = -I/usr/local/Cellar/arpack/3.6.2/libexec/include/arpack


#============================================================
TARGET = u1

C_SOURCES = main.cpp ran2s.cpp rangs.cpp 
C_OBJS    = main.o ran2s.o rangs.o 
MY_INCLUDES = ran2s.h arpack_interface.h

#CXX = pgc++ -fast -Minfo=all -Mprof=ccff
CXX=g++
CXXFLAGS = -g -O2 -std=c++11
#$(LIBDIRS) -larpack $(INC)
#-Xpreprocessor -fopenmp -lomp -DUSE_OMP

#============================================================
all: $(TARGET)

.o:.cpp	$(MY_INCLUDES)
	$(CXX)  -c  $(CXXFLAGS) $<  

$(TARGET) :   $(C_OBJS)
	$(CXX) $(CXXFLAGS)  $^ -o $@

# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites 
#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES)

NOTES = u1StringTension.nb

clean:
	rm -f $(TARGET) $(C_OBJS) core 

TAGS:	$(ALL_SOURCES)
	etags $(ALL_SOURCES)

tar: $(ALL_SOURCES) $(NOTES)
	tar cvf $(TARGET).tar $(ALL_SOURCES)  $(NOTES)

$(TARGET).ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)
