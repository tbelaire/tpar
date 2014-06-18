#include <string>
#include <list>
#include <iostream>
#include <map>

#include "types.h"

// Recognized gates are T, T*, P, P*, Z, Z*, Z a b c, tof a b, tof a, X, H

// Internal representation of a .qc circuit circuit
struct dotqc {
  int n;                   // number of unknown inputs
  int m;                   // number of known inputs (initialized to |0>)
  std::list<std::string> names;      // names of qubits
  std::map<std::string, bool> zero;  // mapping from qubits to 0 (non-zero) or 1 (zero)
  std::vector<std::string> input_wires, output_wires; // the .i and .o lines
  gatelist circ;           // Circuit

  void input(std::istream& in);
  void output(std::ostream& out) const;
  void print() const {output(std::cout);}
  void clear() {n = 0; m = 0; names.clear(); zero.clear(); circ.clear();}
  void append(std::pair<std::string, std::list<std::string> > gate);
  void remove_swaps();
  int count_t_depth() const;
  int count_h() const;
  void print_stats() const;
  void remove_ids();
  bool operator==(const dotqc& other) const;
};
std::ostream& operator<<(std::ostream& out, const dotqc& circuit);

