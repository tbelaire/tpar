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

#include "partition.h"
#include <list>

using namespace std;

template<typename T>
bool is_disjoint(const set<T> & A, const set<T> & B) {
  return 0 == set_intersection_count(A.begin(), A.end(), B.begin(), B.end());
}

ostream& operator<<(ostream& output, const partitioning& part) {

  for (auto Si = part.begin(); Si != part.end(); Si++) {
    output << "{";
    for (auto yi = Si->begin(); yi != Si->end(); yi++) {
      output << *yi << ",";
    }
    output << "}";
  }

  return output;
}

// Take a partition and a set of xor_funcs, and return all partitions that are not
//   disjoint with the set, also removing them from the partition
partitioning freeze_partitions(partitioning & part, set<xor_func> & st) {
  partitioning ret;
  partitioning::iterator it, tmp;

  for (it = part.begin(); it != part.end();) {
    if (!is_disjoint(*it, st)) {
      tmp = it;
      it++;
      ret.splice(ret.begin(), part, tmp);
    } else {
      it++;
    }
  }
  return ret;
}

int num_elts(partitioning & part) {
  int tot = 0;
  partitioning::iterator it;
  for (it = part.begin(); it != part.end(); it++) {
    tot += it->size();
  }
  return tot;
}

partitioning create(const set<xor_func> & st) {
  return partitioning{st};
}

//-------------------------------------- Matroids

// Implements a matroid partitioning algorithm
void add_to_partition(partitioning & ret, xor_func i, const ind_oracle & oracle) {
  partitioning::iterator Si;
  set<xor_func>::iterator yi, zi;

  // The node q contains a queue of paths and an iterator to each node's location.
  //    Each path's first element is the element we grow more paths from.
  //    If x->y is in the path, then we can replace x with y.
  deque<path> node_q;
  path t;
  path_iterator p;
  bool flag;
  map<xor_func, bool> marked;

  // Reset everything
  node_q.clear();
  /* for (const exponent& xpt : elts) { */
  /*   marked[xpt.first] = false; */
  /* } */
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
    set<xor_func> newset;
    newset.insert(i);
    ret.push_front(newset);
  }

}

// Partition the matroid
partitioning partition_matroid(const vector<xor_func> & elts, const ind_oracle & oracle) {
  partitioning ret;

  // For each element of the matroid
  for (int i = 0; i < elts.size(); i++) {
    add_to_partition(ret, elts[i], oracle);
  }
  return ret;
}

void repartition(partitioning & partition, const ind_oracle & oracle ) {
    list<xor_func> acc;

    for (set<xor_func>& part : partition) {
        boost::optional<xor_func> dep = oracle.retrieve_lin_dep(part);
        if(dep) {
            part.erase(*dep);
            acc.push_back(*dep);
        }
    }
    for(const auto& dep : acc) {
        add_to_partition(partition, dep, oracle);
    }
}

