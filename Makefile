BUILD = build
####################
# Sources
vpath %.cpp src
vpath %.o   $(BUILD)
OBJS = partition.o util.o circuit.o main.o

#####################
# Compiler Clang
CXX = /usr/bin/clang++
# CXX = g++-4.8
LDFLAGS =
CPPFLAGS = -O3 -std=c++11 -g


.PHONY: all clean

all: tpar
	echo "Made all"



tpar: $(OBJS)
	$(CXX) $(LDFLAGS) $(CPPFLAGS) -o tpar $(OBJS)

# partition.o: src/partition.cpp
# 	$(CXX) -c $(FLAGS) src/partition.cpp

# util.o: src/util.cpp
# 	$(CXX) -c $(FLAGS) src/util.cpp

# circuit.o: src/circuit.cpp
# 	$(CXX) -c $(FLAGS) src/circuit.cpp

# main.o: src/main.cpp
# 	$(CXX) -c $(FLAGS) src/main.cpp

clean:
	rm *.o
