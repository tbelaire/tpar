####################
# Sources
vpath %.cpp src
vpath %.h   src

# Excludes main.o since tests don't want to link with that.
OBJS := partition.o util.o circuit.o xor_func.o oracle.o dotqc.o
####################

DEBUGFLAGS = -O0 -g

# Pick one or the other
# include compiler.gcc.mk
# Includes some clang++
include compiler.llvm.mk

############# Warnings
CXXFLAGS += -Wall -Werror
# Fixed
CXXFLAGS += -Wunused-function -Wunused-parameter -Wuninitialized -Wunused-variable
# Fixing
CXXFLAGS +=
# Not fixed yet
CXXFLAGS += -Wno-sign-conversion
##################################

LDFLAGS += -lboost_program_options


.PHONY: all clean check

all: tpar test

include test.mk
test : $(OBJS)

tpar: $(OBJS) main.o
	$(CXX) $(LDFLAGS) $(CPPFLAGS) $^ -o tpar

# If we are to put all the built .o files into a build directory,
# we either have to specify $(BUILD)/%.o everywhere,
# or we have to call make from inside the build folder, and
# use vpath to locate the sources.
#
%.o : %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $^ -o $@

clean:
	-rm *.o *.a
