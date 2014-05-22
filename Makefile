####################
# Sources
vpath %.cpp src

OBJS := partition.o util.o circuit.o
####################

include compiler.mk

.PHONY: all clean test check

all: tpar

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
