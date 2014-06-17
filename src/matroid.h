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

#ifndef MATROID_H
#define MATROID_H

#include <vector>
#include <deque>
#include <assert.h>
#include <boost/optional.hpp>

#include "oracle.h"
#include "partition.h"

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



#endif // MATROID_H
