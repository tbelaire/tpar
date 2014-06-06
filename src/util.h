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

#include "types.h"
#include "xor_func.h"


extern bool disp_log;
extern synth_type synth_method;

void print_wires(const std::vector<xor_func>& wires);
int compute_rank(int m, int n, const std::vector<xor_func> bits);
int compute_rank(const std::vector<xor_func> bits);
int compute_rank(int m, int n, const xor_func * bits);
int compute_rank(const std::set<xor_func> & lst);

int to_upper_echelon(int m, int n,
        std::vector<xor_func>& arr,
        std::function<void(int)> do_negate,
        std::function<void(int, int)> do_swap,
        std::function<void(int, int)> do_xor);
gatelist to_upper_echelon(int m, int n,
        const std::vector<xor_func>& bits,
        const std::vector<std::string>& names);
// Pair of matrix versions
void to_upper_echelon(int m, int n,
        const std::vector<xor_func>& bits,
        std::vector<xor_func>& mat);
void to_upper_echelon(int m, int n,
        std::vector<xor_func>& bits,
        std::vector<xor_func>& mat);

void backfill_matrix(int m, int n,
        std::vector<xor_func>& bits,
        std::function<void(int, int)> do_xor);
gatelist to_lower_echelon(const int m, const int n,
        std::vector<xor_func>& bits,
        const std::vector<std::string> names);
void to_lower_echelon(const int m, const int n,
        std::vector<xor_func>& bits,
        std::vector<xor_func>& mat);

gatelist construct_circuit(exponents_set & phase,
    const partitioning & part,
    const std::vector<xor_func> in,
    const std::vector<xor_func> out,
    const int num,
    const int dim,
    const std::vector<std::string> names);

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
