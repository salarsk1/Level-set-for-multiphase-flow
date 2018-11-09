# C++ compiler (and linker)
CXX= g++
# Compiler flags
CFLAGS=-Dnullptr=0
STRICT_WARNINGS=#Wall -Werror -Wno-unused-function -Wno-sign-compare
CFLAGS+=${STRICT_WARNINGS} -std=c++11
# Include libraries
CFLAGS+= -I/usr/local/include -I/opt/mpich3/include \
			-I/home/salar/developments/cfd_561/final/
CFLAGS+= -DNO_GETARG
LFLAGS+=-L/opt/openmpi/lib/ -L/home/salar/developments/cfd_561/final/ \
			-fopenmp
SRC =  $(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)
EXE=exe

all: $(SRC) $(EXE) $(OBJ)
all: CFLAGS+= $(STRICT_WARNINGS)
opt: $(SRC) $(EXE) $(OBJ)
opt: CFLAGS+=-O3
debug: $(SRC) $(EXE) $(OBJ)
debug: CFLAGS+=-g -DDEBUG_MODE
debug: LFLAGS+=-g

profile: $(SRC) $(EXE) $(OBJ)
profile: CFLAGS+=-DGPROF
profile: LFLAGS+=-lprofiler

#tags: $(SRC)
#	@ctags -R .

$(EXE): $(OBJ)
	$(CXX) $(OBJ) $(LFLAGS) -o $@

.cpp.o:
	$(CXX) -c -std=c++0x $(CFLAGS) $< -o $@

clean:
	@rm -f $(OBJ) *.gch tags exe
