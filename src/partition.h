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

#ifndef PARTITION_H
#define PARTITION_H

#include <list>
#include <set>
#include <iostream>

#include "util.h"
#include "matroid.h"

std::ostream& operator<<(std::ostream& output, const partitioning& part);
partitioning freeze_partitions(partitioning & part, std::set<xor_func> & st);

int num_elts(partitioning & part);
partitioning create(const std::set<xor_func> & st);
void add_to_partition(partitioning & ret, xor_func i, const ind_oracle & oracle);
void repartition(partitioning & partition, const ind_oracle & oracle );
partitioning partition_matroid(const std::vector<xor_func> & elts, const ind_oracle & oracle);
#endif
