#include <ostream>

#include "types.h"
#include "util.h"

using namespace std;


xor_func::xor_func(bool neg, initializer_list<int> lst) :
    bitset(lst.size()),
    negated(neg)
{
    int index = 0;
    for (const int i : lst) {
        this->bitset.set(index++, (i > 0) ? 1 : 0);
    }
}
std::ostream& operator<<(std::ostream& out, const xor_func& f) {
    out << (f.is_negated() ? "~" : " ") << f.bitset;
    return out;
}
