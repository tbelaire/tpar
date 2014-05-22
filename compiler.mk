# Shiny new clang
CC = /usr/local/Cellar/llvm-all/3.4.1/bin/clang
CXX = /usr/local/Cellar/llvm-all/3.4.1/bin/clang++

# # Apple Clang
# CC = /usr/bin/clang
# CXX = /usr/bin/clang++

# # GCC 4.8
# CC = gcc-4.8
# CXX = g++-4.8


DEBUGFLAGS = -O0 -g

# Using C++11 and libc++ from llvm3.4.1
LDFLAGS = -lc++abi -lstdc++
CXXFLAGS = -std=c++11 -stdlib=libc++ $(DEBUGFLAGS)
CPPFLAGS = -I /usr/local/Cellar/llvm-all/3.4.1/include/c++/v1
