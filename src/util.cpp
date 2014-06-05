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

#include <map>
#include <cmath>
#include <cassert>
#include <stdexcept>

#include <boost/lexical_cast.hpp>

#include "util.h"
#include "xor_func.h"
#include "oracle.h"

using namespace std;

bool disp_log = true;
synth_type synth_method = PMH;

void print_wires(const vector<xor_func>& wires) {
  for (auto f : wires) {
    for (int j = 0; j < f.size(); j++) {
      if (f.test(j)) cout << "1";
      else           cout << "0";
    }
    cout << "\n";
  }
}

// Commands for making certain circuits
gatelist xor_com(int a, int b, const vector<string> names) {
  return {{"tof", {names[a], names[b]}}};
}

gatelist swap_com(int a, int b, const vector<string> names) {
  return {
      {"tof", {names[a], names[b]}},
      {"tof", {names[b], names[a]}},
      {"tof", {names[a], names[b]}},
  };
}

gatelist x_com(int a, const vector<string> names) {
  return {{"tof", {names[a]}}};
}

// TODO, can this just use to_upper_echelon?
// Make triangular to determine the rank
int compute_rank_dest(vector<xor_func> tmp) {
  int rank = 0;

  if(tmp.size() == 0) { return 0; } // Empty vector has 0 rank
  int func_size = tmp.at(0).size() ;  // - 1?
  // Make triangular
  for (int i = 0; i < func_size; i++) {
    bool flg = false;
    for (int j = rank; j < tmp.size(); j++) {
      if (tmp[j].test(i)) {
        // If we haven't yet seen a vector with bit i set...
        if (!flg) {
          // If it wasn't the first vector we tried, swap to the front
          if (j != rank) swap(tmp[rank], tmp[j]);
          flg = true;
        } else {
          tmp[j] ^= tmp[rank];
        }
      }
    }
    if (flg) rank++;
  }

  return rank;
}


// If they're giving the info to me, might as well check it.
int compute_rank(int m, int n, const vector<xor_func> bits) {
    if( m != bits.size() ) {
        print_wires(bits);
        throw std::logic_error("Bad sizes in compute_rank: m="
                + boost::lexical_cast<string>(m) + ", bits.size()="
                + boost::lexical_cast<string>(bits.size()));
    }
    if( m > 0 && n != bits[0].size() ){
        print_wires(bits);
        throw std::logic_error("Bad sizes in compute_rank: n="
                + boost::lexical_cast<string>(n) + ", bits[0].size()="
                + boost::lexical_cast<string>(bits.size()));
    }
    return compute_rank_dest(bits);
}

int compute_rank(const vector<xor_func> bits) {
    return compute_rank_dest(bits);
}

int compute_rank(int m, int n, const xor_func * bits) {
  int ret;
  // TODO n is the number of bits used in the xor_funcs
  // perhaps do an assert?
  (void) n;

  // Make a copy of the bitset
  vector<xor_func> tmp;
  for(int i = 0; i < m; i++) {
    tmp.push_back(bits[i]);
  }
  ret = compute_rank_dest(tmp);
  return ret;
}

int compute_rank(const set<xor_func> & set) {
  // Make a copy of the bitset
  vector<xor_func> tmp;
  tmp.reserve(set.size());
  for(const xor_func& f : set) {
      tmp.push_back(f);
  }
  return compute_rank_dest(tmp);
}

// Why take a set too?
// TODO AUDIT
int compute_rank(int n, const exponents_set & expnts, const set<xor_func> & lst) {
  int ret;
  int m = lst.size();
  (void) n;

  vector<xor_func> tmp;
  for (const auto& xpt: expnts) {
      // TODO audit.
      if(xpt.second != 0) { // Nonzero coefficient
          tmp.push_back(xpt.first);
      }
  }
  ret = compute_rank_dest(tmp);
  return ret;
}

int to_upper_echelon(int m, int n, vector<xor_func> arr,
        std::function<void(int)> do_negate,
        std::function<void(int, int)> do_swap);
gatelist to_upper_echelon(int m, int n, vector<xor_func> bits,
        const vector<string> names);
void to_upper_echelon(int m, int n, vector<xor_func> bits,
        vector<xor_func> mat);

int to_upper_echelon(int m, int n,
        vector<xor_func> bits,
        std::function<void(int)> do_negate,
        std::function<void(int, int)> do_swap,
        std::function<void(int, int)> do_xor){

  // Clear out all the X gate caused negations
  for (int j = 0; j < m; j++) {
    if (bits[j].is_negated()) {
        bits[j].negate();
        do_negate(j);
    }
  }

  int rank = 0;

  // Make triangular
  for (int i = 0; i < n; i++) {
    bool flg = false;
    for (int j = rank; j < m; j++) {
      if (bits[j].test(i)) {
        // If we haven't yet seen a vector with bit i set...
        if (!flg) {
          // If it wasn't the first vector we tried, swap to the front
          if (j != rank) {
            swap(bits[rank], bits[j]);
            do_swap(j, rank);
          }
          flg = true;
        } else {
          bits[j] ^= bits[rank];
          do_xor(j, rank);
        }
      }
    }
    if (flg) rank++;
  }
  return rank;
}

// Make echelon form
gatelist to_upper_echelon(int m, int n, vector<xor_func> bits, const vector<string> names) {
  gatelist acc;
  to_upper_echelon(m, n, bits,
          [&acc, &names](int j){
            acc.splice(acc.end(), x_com(j, names));
          },
          [&acc, &names](int r1, int r2){
            acc.splice(acc.end(), swap_com(r1, r2, names));
          },
          [&acc, &names](int target, int i){
            acc.splice(acc.end(), xor_com(target, i, names));
          });
  return acc;
}


// Make echelon form
void to_upper_echelon(int m, int n, vector<xor_func> bits, vector<xor_func> mat) {
  to_upper_echelon(m, n, bits,
          [&mat, m](int j){
            mat[j].set(m);
          },
          [&mat](int r1, int r2){
            swap(mat[r1], mat[r2]);
          },
          [&mat](int target, int i){
            mat[target] ^= mat[i];
          });
}

// TODO audit
gatelist to_lower_echelon(int m, int n, vector<xor_func> bits, const vector<string> names) {
  gatelist acc;

  (void) m;
  for (int i = n-1; i > 0; i--) {
    for (int j = i - 1; j >= 0; j--) {
      if (bits[j].test(i)) {
        bits[j] ^= bits[i];
        acc.splice(acc.end(), xor_com(i, j, names));
      }
    }
  }
  return acc;
}

// TODO audit
void to_lower_echelon(int m, int n, vector<xor_func>& bits, vector<xor_func>& mat) {
  (void) m;
  for (int i = n-1; i > 0; i--) {
    for (int j = i - 1; j >= 0; j--) {
      if (bits[j].test(i)) {
        bits[j] ^= bits[i];
        mat[j] ^= mat[i];
      }
    }
  }
}

// Existing code version with bool instead of NULLs
gatelist fix_basis(int m, int n, int k, const xor_func * fst, xor_func * snd, xor_func * mat, const vector<string> names, bool has_mat);

// Fixed interface versions
gatelist fix_basis(int m, int n, int k, const vector<xor_func> fst, vector<xor_func> snd, const vector<string> names);
void fix_basis(int m, int n, int k, const vector<xor_func> fst, vector<xor_func> snd, vector<xor_func> mat);

gatelist fix_basis(int m, int n, int k, const vector<xor_func> fst, vector<xor_func> snd, const vector<string> names) {
    return fix_basis(m, n, k, &fst[0], &snd[0], nullptr, names, false);
}
void fix_basis(int m, int n, int k, const vector<xor_func> fst, vector<xor_func> snd, vector<xor_func> mat) {
    fix_basis(m, n, k, &fst[0], &snd[0], &mat[0], {}, true);
}
// Expects two matrices in echelon form, the second being a subset of the
//   rowspace of the first. It then morphs the second matrix into the first
gatelist fix_basis(int m, int n, int k, const xor_func * fst, xor_func * snd, xor_func * mat, const vector<string> names, bool has_mat) {
  gatelist acc;
  int j = 0;
  bool flg = false;
  map<int, int> pivots;  // mapping from columns to rows that have that column as pivot
  for (int i = 0; i < n; i++) pivots[i] = -1;

  // First pass makes sure tmp has the same pivots as fst
  for (int i = 0; i < m; i++) {
    // Find the next pivot
    while (j < n && !fst[i].test(j)) j++;
    if (j < n) {
      pivots[j] = i;
      flg = false;
      for (int h = i; !flg && h < k; h++) {
        // We found a vector with the same pivot
        if (snd[h].test(j)) {
          flg = true;
          if (h != i) {
            swap(snd[h], snd[i]);
            if (!has_mat)  acc.splice(acc.end(), swap_com(h, i, names));
            else           swap(mat[h], mat[i]);
          }
        }
      }
      // There was no vector with the same pivot
      if (!flg) {
        if (k >= m) {
          cout << "FATAL ERROR: second space not a subspace\n" << flush;
          exit(1);
        }
        snd[k] = fst[i];
        if (k != i) {
          swap(snd[k], snd[i]);
          if (!has_mat) acc.splice(acc.end(), swap_com(k, i, names));
          else          swap(mat[k], mat[i]);
        }
        k++;
      }
    }
  }

  // Second pass makes each row of tmp equal to that row of fst
  for (int i = 0; i < m; i++) {
    for (int j = i +1; j < n; j++) {
      if (fst[i][j] != snd[i][j]) {
        if (pivots[j] == -1) {
          cout << "FATAL ERROR: cannot fix basis\n" << flush;
          exit(1);
        } else {
          snd[i] ^= snd[pivots[j]];
          if (! has_mat) acc.splice(acc.end(), xor_com(pivots[j], i, names));
          else           mat[i] ^= mat[pivots[j]];
        }
      }
    }
    if (!(snd[i] == fst[i])) {
      cout << "FATAL ERROR: basis differs\n" << flush;
      exit(1);
    }
  }

  return acc;
}

// A := B^{-1} A
void compose(int num, vector<xor_func> A, const vector<xor_func> B) {
  vector<xor_func> tmp = B;
  to_upper_echelon(num, num, tmp, A);
  to_lower_echelon(num, num, tmp, A);
}

//------------------------- CNOT synthesis methods

// Gaussian elimination based CNOT synthesis
gatelist gauss_CNOT_synth(int n, int m, vector<xor_func> bits, const vector<string> names) {
  int i, j, k;
  gatelist lst;
  list<string> tmp_list1, tmp_list2;

  for (j = 0; j < n; j++) {
    if (bits[j].test(n)) {
      bits[j].reset(n);
      lst.splice(lst.begin(), x_com(j, names));
    }
  }

  // Make triangular
  for (i = 0; i < n; i++) {
    bool flg = false;
    for (j = i; j < n + m; j++) {
      if (bits[j].test(i)) {
        // If we haven't yet seen a vector with bit i set...
        if (!flg) {
          // If it wasn't the first vector we tried, swap to the front
          if (j != i) {
            swap(bits[i], bits[j]);
            lst.splice(lst.begin(), swap_com(i, j, names));
          }
          flg = true;
        } else {
          bits[j] ^= bits[i];
          lst.splice(lst.begin(), xor_com(i, j, names));
        }
      }
    }
    if (!flg) {
      cout << "ERROR: not full rank\n";
      exit(1);
    }
  }

  //Finish the job
  for (i = n-1; i > 0; i--) {
    for (j = i - 1; j >= 0; j--) {
      if (bits[j].test(i)) {
        bits[j] ^= bits[i];
        lst.splice(lst.begin(), xor_com(i, j, names));
      }
    }
  }

  return lst;
}

// Patel/Markov/Hayes CNOT synthesis
gatelist Lwr_CNOT_synth(int n, int m, vector<xor_func> bits, const vector<string> names, bool rev) {
  gatelist acc;
  int sec, tmp, row, col, i;
  int patt[1<<m];

  for (sec = 0; sec < ceil(n / m); sec++) {

    for (i = 0; i < (1<<m); i++) {
      patt[i] = -1;
    }
    for (row = sec*m; row < n; row++) {
      tmp = 0;
      for (i = 0; i < m; i++) {
        if (bits[row].test(sec*m + i)) tmp += (1 << i);
      }
      if (patt[tmp] == -1) {
        patt[tmp] = row;
      } else if (tmp != 0) {
        bits[row] ^= bits[patt[tmp]];
        if (rev) acc.splice(acc.begin(), xor_com(row, patt[tmp], names));
        else acc.splice(acc.end(), xor_com(patt[tmp], row, names));
      }
    }

    for (col = sec*m; col < (sec+1)*m; col++) {
      for (row=col + 1; row < n; row++) {
        if (bits[row].test(col)) {
          if (not(bits[col].test(col))) {
            bits[col] ^= bits[row];
            bits[row] ^= bits[col];
            bits[col] ^= bits[row];
            if (rev) {
              acc.splice(acc.begin(), xor_com(col, row, names));
              acc.splice(acc.begin(), xor_com(row, col, names));
              acc.splice(acc.begin(), xor_com(col, row, names));
            } else {
              acc.splice(acc.end(), xor_com(row, col, names));
              acc.splice(acc.end(), xor_com(col, row, names));
              acc.splice(acc.end(), xor_com(row, col, names));
            }
          } else {
            bits[row] ^= bits[col];
            if (rev) acc.splice(acc.begin(), xor_com(row, col, names));
            else acc.splice(acc.end(), xor_com(col, row, names));
          }
        }
      }
    }
  }

  return acc;
}

gatelist CNOT_synth(int n, vector<xor_func> bits, const vector<string> names) {
  gatelist acc, tmp;
  int i, j, m = (int)(log((double)n) / (log(2) * 2));

  for (j = 0; j < n; j++) {
    if (bits[j].test(n)) {
      bits[j].reset(n);
      acc.splice(acc.begin(), x_com(j, names));
    }
  }

  acc.splice(acc.end(), Lwr_CNOT_synth(n, m, bits, names, false));
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      bits[j].set(i, bits[i][j]); // Transpose
      bits[i].reset(j);
    }
  }
  acc.splice(acc.end(), Lwr_CNOT_synth(n, m, bits, names, true));
  acc.reverse();

  return acc;
}

// Construct a circuit for a given partition
gatelist construct_circuit(exponents_set & phase,
    const partitioning & part,
    const vector<xor_func> in,
    const vector<xor_func> out,
    const int num,
    const int dim,
    const vector<string> names) {
  gatelist ret, tmp, rev;
  vector<xor_func> bits, pre, post;
  bool flg = true;

  assert(in.size() == num);
  assert(out.size() == num);
  // flg = forall i. in[i] == out[i]
  for (int i = 0; i < num; i++) {
    flg &= (in[i] == out[i]);
  }
  for (int i = 0; i < num; i++) {
    /* flg &= (in[i] == out[i]); */
    if (synth_method != AD_HOC) {
      pre.emplace_back(xor_func(num));
      post.emplace_back(xor_func(num));
      pre[i].set(i);
      post[i].set(i);
    }
  }
  if (flg && (part.size() == 0)) return ret;

  // Reduce in to echelon form to decide on a basis
  if (synth_method == AD_HOC) ret.splice(ret.end(), to_upper_echelon(num, dim, in, names));
  else to_upper_echelon(num, dim, in, pre);

  // For each partition... Compute *it, apply T gates, uncompute
  for (partitioning::const_iterator it = part.begin(); it != part.end(); it++) {
    cerr << "List is of length: " << it->size() << endl;
    auto ti = it->begin();
    for (int i = 0; i < num;  i++) {
      if (ti != it->end()){
          it++;
          /* cerr << "Ran out of list" << endl; */
          /* exit(4); */
      }
      if (i < it->size()) { bits[i] = *ti; }
      else                { bits[i] = xor_func(dim); }
    }

    // prepare the bits
    if (synth_method == AD_HOC) {
      tmp = to_upper_echelon(it->size(), dim, bits, names);
      tmp.splice(tmp.end(), fix_basis(num, dim, it->size(), in, bits, names));
      rev = tmp;
      rev.reverse();
      ret.splice(ret.end(), rev);
    } else {
      to_upper_echelon(it->size(), dim, bits, post);
      fix_basis(num, dim, it->size(), in, bits, post);
      compose(num, pre, post);
      if (synth_method == GAUSS) ret.splice(ret.end(), gauss_CNOT_synth(num, 0, pre, names));
      else if (synth_method == PMH) ret.splice(ret.end(), CNOT_synth(num, pre, names));
    }

    // apply the T gates
    list<string> tmp_lst;
    ti = it-> begin();
    for (int i = 0; ti != it->end(); ti++, i++) {
      tmp_lst.clear();
      tmp_lst.push_back(names[i]);
      if (phase[*ti] <= 4) {
        if (phase[*ti] / 4 == 1) ret.push_back(make_pair("Z", tmp_lst));
        if (phase[*ti] / 2 == 1) ret.push_back(make_pair("P", tmp_lst));
        if (phase[*ti] % 2 == 1) ret.push_back(make_pair("T", tmp_lst));
      } else {
        if (phase[*ti] == 5 || phase[*ti] == 6) ret.push_back(make_pair("P*", tmp_lst));
        if (phase[*ti] % 2 == 1) ret.push_back(make_pair("T*", tmp_lst));
      }
    }

    // unprepare the bits
    if (synth_method == AD_HOC) ret.splice(ret.end(), tmp);
    else {
      pre = post;
      // re-initialize post
      for (int i = 0; i < num; i++) {
        post[i] = xor_func(num);
        post[i].set(i);
      }
    }
  }

  // Reduce out to the basis of in
  for (int i = 0; i < num; i++) {
    bits[i] = out[i];
  }
  if (synth_method == AD_HOC) {
    tmp = to_upper_echelon(num, dim, bits, names);
    tmp.splice(tmp.end(), fix_basis(num, dim, num, in, bits, names));
    tmp.reverse();
    ret.splice(ret.end(), tmp);
  } else {
    to_upper_echelon(num, dim, bits, post);
    fix_basis(num, dim, num, in, bits, post);
    compose(num, pre, post);
    if (synth_method == GAUSS) ret.splice(ret.end(), gauss_CNOT_synth(num, 0, pre, names));
    else if (synth_method == PMH) ret.splice(ret.end(), CNOT_synth(num, pre, names));
  }

  return ret;
}

vector<xor_func> init_matrix(
        initializer_list<initializer_list<int>> lst){
    vector<xor_func> cols{};
    for (const auto& l : lst)
     {
         cols.push_back({false, l});
     }
    return cols;
}

list_compare_result
list_compare(const list<string> & a, const list<string> & b) {
    const set<string> set_a{a.begin(), a.end()},
                      set_b{b.begin(), b.end()};
    const auto overlap = set_intersection_count(
            set_a.begin(), set_a.end(),
            set_b.begin(), set_b.end());
    if (overlap == 0) {
        return list_compare_result::DISJOINT;
    } else if (overlap == set_a.size() && overlap == set_b.size()) {
        return list_compare_result::EQUAL;
    } else {
        return list_compare_result::OVERLAPPED;
    }

}
