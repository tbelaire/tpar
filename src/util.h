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

#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "partition.h"

typedef boost::dynamic_bitset<>            xor_func;
typedef std::pair<char, xor_func >              exponent;
typedef std::list<std::pair<std::string, std::list<std::string> > > gatelist; // [(Str, [Str])]

enum synth_type { AD_HOC, GAUSS, PMH };

extern bool disp_log;
extern synth_type synth_method;

class ind_oracle {
  private:
    int num;
    int dim;
    int length;
  public:
    ind_oracle() { num = 0; dim = 0; length = 0; }
    ind_oracle(int numin, int dimin, int lengthin) { num = numin; dim = dimin; length = lengthin; }

    void set_dim(int newdim) { dim = newdim; }
    int retrieve_lin_dep(const std::vector<exponent> & expnts, const std::set<int> & lst) const;

    bool operator()(const std::vector<exponent> & expnts, const std::set<int> & lst) const;
};

void print_wires(const xor_func * wires, int num, int dim);
int compute_rank(int m, int n, const std::vector<xor_func> bits);
int compute_rank(int m, int n, const xor_func * bits);
int compute_rank(int n, const std::vector<exponent> & expnts, const std::set<int> & lst);

gatelist construct_circuit(const std::vector<exponent> & phase,
    const partitioning & part,
    const std::vector<xor_func> in,
    const std::vector<xor_func> out,
    const int num,
    const int dim,
    const std::vector<std::string> names);

xor_func init_xor_func(std::initializer_list<int> lst);
std::vector<xor_func> init_matrix_transpose(
        std::initializer_list<std::initializer_list<int>>);

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
