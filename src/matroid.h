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
#include <deque>
#include <assert.h>
#include <boost/optional.hpp>

#include "oracle.h"
#include "partition.h"


#ifndef MATROID
#define MATROID

struct path {
  std::list<std::pair <xor_func, partitioning::iterator> > lst;

  path() { }
  path(xor_func i, partitioning::iterator ref) { lst.push_front(make_pair(i, ref)); }
  path(const path & p) { lst = p.lst; }
  path(xor_func i, partitioning::iterator ref, path & p) {
    lst = p.lst;
    lst.push_front(make_pair(i, ref));
  }

  std::pair<xor_func, partitioning::iterator> head() { return lst.front(); }
  xor_func                               head_elem() { return lst.front().first; }
  partitioning::iterator                 head_part() { return lst.front().second; }

  path_iterator begin() { return lst.begin(); }
  path_iterator   end() { return lst.end(); }

  void insert(xor_func i, partitioning::iterator ref) { lst.push_front(make_pair(i, ref)); }
};

//-------------------------------------- Matroids

// Implements a matroid partitioning algorithm
template <class cont, typename oracle_type>
void add_to_partition(partitioning & ret, xor_func i, const cont & elts, const oracle_type & oracle) {
  partitioning::iterator Si;
  std::set<xor_func>::iterator yi, zi;

  // The node q contains a queue of paths and an iterator to each node's location.
  //    Each path's first element is the element we grow more paths from.
  //    If x->y is in the path, then we can replace x with y.
  std::deque<path> node_q;
  path t;
  path_iterator p;
  bool flag;
  std::map<xor_func, bool> marked;

  // Reset everything
  node_q.clear();
  for (const exponent& xpt : elts) {
    marked[xpt.first] = false;
  }
  flag = false;

  // Insert element to be partitioned
  node_q.push_back(path(i, ret.end()));
  marked[i] = true;

  // BFS loop
  while (!node_q.empty() && !flag) {
    // The head of the path is what we're currently considering
    t = node_q.front();
    node_q.pop_front();

    for (Si = ret.begin(); Si != ret.end() && !flag; Si++) {
      if (Si != t.head_part()) {
        // Add the head to Si. If Si is independent, leave it, otherwise we'll have to remove it
        Si->insert(t.head_elem());

        if (oracle(*Si)) {
          // We have the shortest path to a partition, so make the changes:
          //	For each x->y in the path, remove x from its partition and add y
          for (p = t.begin(); p != --(t.end()); ) {
            Si = p->second;
            (Si)->erase(p->first);
            (Si)->insert((++p)->first);
          }
          flag = true;
        } else {
          // For each element of Si, if removing it makes an independent set, add it to the queue
          for (yi = Si->begin(); yi != Si->end(); yi++) {
            if (!marked[*yi]) {
              // Generate an iterator to the position before yi
              zi = yi;
              if (zi != Si->begin()) zi--;
              // Take yi out
              auto tmp = *yi;
              Si->erase(yi);
              if (oracle(*Si)) {
                // Put yi back in
                yi = Si->insert(Si->begin(), tmp);
                // Add yi to the queue
                node_q.push_back(path(*yi, Si, t));
                marked[*yi] = true;
              } else {
                yi = Si->insert(Si->begin(), tmp);
              }
            }
          }
          // Remove CURRENT from Si
          Si->erase(t.head_elem());
        }
      }
    }
  }

  // We were unsuccessful trying to edit the current partitions
  if (!flag) {
    std::set<xor_func> newset;
    newset.insert(i);
    ret.push_front(newset);
  }

}

// Partition the matroid
template <class T, typename oracle_type>
partitioning partition_matroid(const std::vector<T> & elts, const oracle_type & oracle) {
  partitioning ret;

  // For each element of the matroid
  for (int i = 0; i < elts.size(); i++) {
    add_to_partition(ret, elts[i], elts, oracle);
  }
  return ret;
}

// Not so polymorphic, sorry
template <typename oracle_type>
void repartition(partitioning & partition, exponents_set & elts, const oracle_type & oracle ) {
    std::list<xor_func> acc;

    for (std::set<xor_func>& part : partition) {
        boost::optional<xor_func> dep = oracle.retrieve_lin_dep(part);
        if(dep) {
            part.erase(*dep);
            acc.push_back(*dep);
        }
    }
    for(const auto& dep : acc) {
        add_to_partition(partition, dep, elts, oracle);
    }
}

template <class T, typename oracle_type>
void repartition(partitioning & part, const std::vector<T> & elts, const oracle_type & oracle) {
  partitioning::iterator Si;
  std::set<xor_func>::iterator yi;

  std::list<xor_func> acc;

  for (Si = part.begin(); Si != part.end(); Si++) {
    boost::optional<xor_func> tmp = oracle.retrieve_lin_dep(*Si);
    if (tmp) {
      Si->erase(*tmp);
      acc.push_back(*tmp);
    }
    //assert(oracle(*Si));
  }

  for (std::list<xor_func>::iterator it = acc.begin(); it != acc.end(); it++) {
    add_to_partition(part, *it, elts, oracle);
  }
}

#endif
