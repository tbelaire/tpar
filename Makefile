####################
# Sources
vpath %.cpp src
vpath %.h   src

# Excludes main.o since tests don't want to link with that.
OBJS := partition.o util.o circuit.o
####################

include compiler.mk
CXXFLAGS += -Wall
# Fixed
CXXFLAGS += -Wheader-hygiene
# Fixing
CXXFLAGS +=
# Not fixed yet
CXXFLAGS += -Wno-sign-conversion -Wno-unused-variable
# Don't care about these
CXXFLAGS += -Wno-c99-extensions -Wno-c++98-compat

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
