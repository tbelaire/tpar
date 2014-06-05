#include "oracle.h"
#include "util.h"
#include <boost/optional.hpp>

using namespace std;


bool ind_oracle::operator()(const set<xor_func> & lst) const {
  if (lst.size() > this->num) return false;
  if (lst.size() == 1 || (this->num - lst.size()) >= this->dim) return true;

  set<xor_func>::const_iterator it;
  int i, j, rank = 0;
  bool flg;
  vector<xor_func> tmp{};

  for (const xor_func& f : lst) {
    tmp.push_back(f);
  }

  rank = compute_rank(tmp);

  return (this->num - lst.size()) >= (this->dim - rank);
}

//TODO audit this
// Shortcut to find a linearly dependent element faster
boost::optional<xor_func>
ind_oracle::retrieve_lin_dep(const set<xor_func> & lst) const {
  set<xor_func>::const_iterator it;
  map<int, xor_func> mp;
  bool flg;
  vector<xor_func> tmp;

  int rank = 0;

  {
    int i = 0;
    for (const xor_func& f : lst) {
        tmp.push_back(f);
        mp.emplace(i, f);
        i++;
    }
  }

  for (int j = 0; j < lst.size(); j++) {
    if (tmp[j].test(length)) tmp[j].reset(length);
  }

  for (int i = 0; i < length; i++) {
    flg = false;
    for (int j = rank; j < lst.size(); j++) {
      if (tmp[j].test(i)) {
        // If we haven't yet seen a vector with bit i set...
        if (!flg) {
          // If it wasn't the first vector we tried, swap to the front
          if (j != rank) {
            swap(tmp[rank], tmp[j]);
            xor_func tmpr = mp.at(rank);
            mp.at(rank) = mp.at(j);
            mp.at(j) = tmpr;
          }
          flg = true;
        } else {
          tmp[j] ^= tmp[rank];
          if (tmp[j].none()) return mp.at(j);
        }
      }
    }
    if (flg) rank++;
  }

  assert((num - lst.size()) >= (dim - rank));
  return boost::optional<xor_func>();
}

