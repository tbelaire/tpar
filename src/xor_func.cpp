#include <ostream>
#include <iostream>

#include "types.h"
#include "util.h"

using namespace std;


xor_func::xor_func(const bool neg, const initializer_list<int> lst) :
    bitset(lst.size()),
    negated(neg)
{
    int index = 0;
    for (const int i : lst) {
        this->bitset.set(index++, (i > 0) ? 1 : 0);
    }
}
std::ostream& operator<<(std::ostream& out, const xor_func& f) {
    out << (f.is_negated() ? "~" : " ");
    for (int i = 0; i < f.size(); i++) {
        out << (int) f.test(i);
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, const std::vector<xor_func>& arr) {
    for (const xor_func& f : arr) {
        out << f << endl;
    }
    return out;
}

void output_with_names(std::ostream& out, xor_func f, std::map<int, int> val_map, std::vector<std::string> names) {

    bool following = false;
    out << "[";
    for (int i = 0; i < f.size(); i++) {
      if (f.test(i)) {
          if(following) { out << ", "; }
          out << names.at(val_map.at(i));
          following = true;
      }
    }
    out << "]";
}

// Returns a int in range [0, 2^offset -1]
// It will be different for each xor_func,
// beyond that there are no guarantees
unsigned long
xor_func::slice(size_t start, size_t offset) {
    const size_t left_edge = offset;
    const size_t right_edge = start;
    const unsigned long bits = bitset.to_ulong();
    /* cout << "slicing " <<  bits << "["; */
    /* cout << left_edge <<"..." << right_edge << "]" << endl; */
    // Line up bits[start+offset] with result[0]
    const unsigned long aligned = bits >> (right_edge);
    // Mask out all the higher bits
    const unsigned long mask = ~((~0) << left_edge);
    return aligned & mask;
}


void extend_row_length(vector<xor_func>& arr, int length) {
    for(xor_func& f : arr) {
        f.resize(length);
    }
}
