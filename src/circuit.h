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

#include <string>
#include <list>
#include <iostream>
#include <set>
#include <map>
#include "matroid.h"
#include "util.h"

// Recognized gates are T, T*, P, P*, Z, Z*, Z a b c, tof a b, tof a, X, H

// Internal representation of a .qc circuit circuit
struct dotqc {
  int n;                   // number of unknown inputs
  int m;                   // number of known inputs (initialized to |0>)
  std::list<std::string> names;      // names of qubits
  std::map<std::string, bool> zero;  // mapping from qubits to 0 (non-zero) or 1 (zero)
  gatelist circ;           // Circuit

  void input(std::istream& in);
  void output(std::ostream& out) const;
  void print() const {output(std::cout);}
  void clear() {n = 0; m = 0; names.clear(); zero.clear(); circ.clear();}
  void append(std::pair<std::string, std::list<std::string> > gate);
  void remove_swaps();
  int count_t_depth() const;
  void print_stats() const;
  void remove_ids();
  bool operator==(const dotqc& other) const;
};
std::ostream& operator<<(std::ostream& out, const dotqc& circuit);

// ------------------------- Hadamard version
struct Hadamard {
  int qubit;        // Which qubit this hadamard is applied to
  int prep;         // Which "value" this hadamard prepares

  std::set<xor_func> in;      // exponent terms that must be prepared before the hadamard
  std::vector<xor_func> wires; // state of the wires when this hadamard is applied
};

// Characteristic of a circuit
struct character {
  const int n;                        // number of unknown inputs
  int m;                        // number of zero-initialized ancilla qubits
  const int h;                        // number of hadamards
  std::vector<std::string>   names;  // names of qubits
  std::vector<bool>     zero;        // Which qubits start as 0
  std::map<int, int>    val_map;     // which value corresponds to which qubit
  exponents_set         phase_expts; // a list of exponents of \omega in the mapping
  std::vector<xor_func> outputs;     // the xors computed into each qubit
  // TODO: make this a dependency graph instead
  std::list<Hadamard>   hadamards;   // a list of the hadamards in the order we saw them
  character(const dotqc &input);
  void output(std::ostream& out) const;
  void print() {output(std::cout);}
  void print_outputs() const;
  void add_ancillae(int num);
  dotqc synthesize();
  dotqc synthesize_unbounded();
};

// ------------------------- {CNOT, T} version

enum circuit_type { CNOTT, OTHER, UNKNOWN };

struct metacircuit {
  int n;                        // number of unknown inputs
  int m;                        // number of known inputs (initialized to |0>)
  std::list<std::string> names;           // names of qubits
  std::map<std::string, bool> zero;       // mapping from qubits to 0 (non-zero) or 1 (zero)
  std::list<std::pair<circuit_type, dotqc> > circuit_list;     // A list of subcircuits

  void partition_dotqc(dotqc & input);
  void output(std::ostream& out);
  void print() {output(std::cout);}
  void optimize();
  dotqc to_dotqc();
};
