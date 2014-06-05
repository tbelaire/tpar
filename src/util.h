/*--------------------------------------------------------------------
  Tpar - T-gate optimization for quantum circuits
  Copyright (C) 2013  Matthew Amy and The University of Waterloo,
  Institute for Quantum Computing, Quantum Circuits Group

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

Author: Matthew Amy
---------------------------------------------------------------------*/
#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <list>
#include <map>
#include <set>
#include <boost/dynamic_bitset.hpp>
#include <boost/optional.hpp>

class xor_func;
using exponent_val = unsigned char;
using exponent = std::pair<xor_func, exponent_val>;
using exponents_set = std::map<xor_func, exponent_val>;

// [(Str, [Str])]
using gatelist = std::list<std::pair<std::string, std::list<std::string>>>;

using partitioning = std::list<std::set<xor_func>>;
using path_iterator = std::list<std::pair<xor_func, partitioning::iterator>>::iterator;

enum synth_type { AD_HOC, GAUSS, PMH };

extern bool disp_log;
extern synth_type synth_method;

// NOTE, the first n+h bits describe the combination of inputs,
// and the final bit represents if it's negated or not.
class xor_func {
    private:
        boost::dynamic_bitset<> bitset;
        bool negated;
    public:
        xor_func(bool neg, std::initializer_list<int> lst);
        xor_func(size_t size) : bitset(size), negated(false) {}
        /* xor_func() : bitset(), negated(false) {} */
        xor_func(bool negated, boost::dynamic_bitset<> bitset)
            : bitset(bitset), negated(negated) {}

        bool is_negated() const { return negated; }
        void negate() { negated = !negated; }

        xor_func operator^(const xor_func& b) const {
            return xor_func(negated ^ b.negated,
                    bitset ^ b.bitset);
        }

        xor_func& operator^=(const xor_func& b) {
            this->bitset ^= b.bitset;
            this->negated ^= b.negated; // double check this
            return *this;
        }
        bool operator==(const xor_func& b) const {
            return this->bitset == b.bitset
                && this->negated == b.negated;
        }
        bool operator<(const xor_func& b) const {
            if(this->negated == b.negated) {
                return this->bitset < b.bitset;
            } else {
                return this->negated < b.negated;
            }
        }
        friend std::ostream& operator<<(std::ostream& out, const xor_func& f);
        bool contains(const xor_func& b) const {
            // Explicitly don't care about `negated`.
            return ((~ this->bitset) & b.bitset).none();
        }
        // Pass through functions
        size_t size() const { return bitset.size(); }
        bool test(const size_t i) const { return bitset.test(i); }
        void set(const size_t i, const bool val = true) {
            bitset.set(i, val);
        }
        void reset(const size_t i) { bitset.reset(i); }
        void reset() { bitset.reset(); }
        void flip(const size_t i) { bitset.flip(i); }
        bool none() const { return bitset.none(); }
        bool any() const { return bitset.any(); }

        bool operator[](const size_t& i) const {
            return this->bitset[i];
        }
};

class ind_oracle {
  private:
    int num;
    int dim;
    int length;
  public:
    ind_oracle() { num = 0; dim = 0; length = 0; }
    ind_oracle(int numin, int dimin, int lengthin) { num = numin; dim = dimin; length = lengthin; }

    void set_dim(int newdim) { dim = newdim; }
    boost::optional<xor_func> retrieve_lin_dep(const std::set<xor_func> & lst) const;

    bool operator()(const std::set<xor_func> & lst) const;
};

void print_wires(const xor_func * wires, int num, int dim);
int compute_rank(int m, int n, const std::vector<xor_func> bits);
int compute_rank(const std::vector<xor_func> bits);
int compute_rank(int m, int n, const xor_func * bits);
int compute_rank(int n, const exponents_set & expnts, const std::set<xor_func> & lst);
int compute_rank(const std::set<xor_func> & lst);

gatelist construct_circuit(exponents_set & phase,
    const partitioning & part,
    const std::vector<xor_func> in,
    const std::vector<xor_func> out,
    const int num,
    const int dim,
    const std::vector<std::string> names);

xor_func init_xor_func(std::initializer_list<int> lst);
std::vector<xor_func> init_matrix(
        std::initializer_list<std::initializer_list<int>>);

// Note, requies the inputs to be sorted.
template<class InputIt1, class InputIt2>
size_t set_intersection_count(InputIt1 first1, InputIt1 last1,
                          InputIt2 first2, InputIt2 last2)
{
    size_t count = 0;
    while (first1 != last1 && first2 != last2) {
        if (*first1 < *first2) {
            ++first1;
        } else  {
            if (!(*first2 < *first1)) {
                count++;
                first1++;
            }
            ++first2;
        }
    }
    return count;
}

enum class list_compare_result { EQUAL, DISJOINT, OVERLAPPED };
list_compare_result
list_compare(const std::list<std::string> & a, const std::list<std::string> & b);
#endif // UTIL_H
