#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <list>
#include <map>
#include <set>

class xor_func;

using exponent_val = unsigned char;
using exponent = std::pair<xor_func, exponent_val>;
using exponents_set = std::map<xor_func, exponent_val>;

// [(Str, [Str])]
using gatelist = std::list<std::pair<std::string, std::list<std::string>>>;

using partitioning = std::list<std::set<xor_func>>;
using path_iterator = std::list<std::pair<xor_func, partitioning::iterator>>::iterator;

enum synth_type { AD_HOC, GAUSS, PMH };
#endif // TYPES_H
