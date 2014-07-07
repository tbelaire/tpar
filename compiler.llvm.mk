CC = /usr/local/Cellar/llvm-all/3.4.1/bin/clang
CXX = /usr/local/Cellar/llvm-all/3.4.1/bin/clang++

# Apple Clang
# CC = /usr/bin/clang
# CXX = /usr/bin/clang++

# Using C++11 and libc++ from llvm3.4.1
LDFLAGS = -lc++abi -lstdc++
CXXFLAGS = -std=c++11 -stdlib=libc++ $(DEBUGFLAGS)
CPPFLAGS = -I /usr/local/Cellar/llvm-all/3.4.1/include/c++/v1


# Warnings that a llvm only
CXXFLAGS += -Wheader-hygiene

# Don't care about these
CXXFLAGS += -Wno-c99-extensions -Wno-c++98-compat

# Turn on runtime checks for memory errors
# CXXFLAGS += -fsanitize=address
