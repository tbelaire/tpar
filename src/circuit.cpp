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

#include <algorithm>
#include <sstream>

#include <boost/optional.hpp>

#include "circuit.h"
#include "xor_func.h"
#include "oracle.h"

using namespace std;

// Forward decl to use in construct
void insert_phase (unsigned char c, xor_func f, map<xor_func, unsigned char> & phases);
// Parse a {CNOT, T} circuit
// NOTE: a qubit's number is NOT the same as the bit it's value represents
character::character(const dotqc &input) :
    n(input.n),
    m(input.m),
    h(input.count_h()),
    names(n + m + h),
    zero(n + m),
    val_map(),
    phase_expts(),
    outputs(),
    hadamards()
{
    int val_max = 0;
    map<string, int> name_map, gate_lookup;
    gate_lookup["T"] = 1;
    gate_lookup["T*"] = 7;
    gate_lookup["P"] = 2;
    gate_lookup["P*"] = 6;
    gate_lookup["Z"] = 4;
    gate_lookup["Y"] = 4; // TODO investigate

    // Initialize names and wires
    vector<xor_func> wires;
    wires.reserve( n + m );
    int name_max = 0;
    for (const auto& name : input.names) {
        // name_map maps a name to a wire
        name_map[name] = name_max;
        // names maps a wire to a name
        names[name_max] = name;
        // zero mapping
        zero[name_max]  = input.zero.at(name);
        // each wire has an initial value j, unless it starts in the 0 state
        wires.push_back(xor_func(n + h));
        if (!zero[name_max]) {
            wires[name_max].set(val_max);
            val_map[val_max++] = name_max;
        }
        name_max++;
    }

    // map<string, int> name_map
    // type gatelist = [(Str, [Str]]
    int gate_index = -1;
    for (const auto& gate : input.circ ) {
        bool flg = false;
        gate_index++;
        if (gate.first == "tof" && gate.second.size() == 2) {
            auto gate_inputs_it = gate.second.begin();
            int first_input = name_map[*gate_inputs_it];
            gate_inputs_it++;
            int second_input = name_map[*gate_inputs_it];
            wires[second_input] ^= wires[first_input];
        } else if ((gate.first == "tof" || gate.first == "X")
                && gate.second.size() == 1) {
            int first_input = name_map[*(gate.second.begin())];
            wires[first_input].flip(n + h);
        } else if (gate.first == "Y" && gate.second.size() == 1) {
            int first_input = name_map[*(gate.second.begin())];
            insert_phase(gate_lookup[gate.first], wires[first_input], phase_expts);
            wires[name_map[*(gate.second.begin())]].flip(n + h);
        } else if (gate.first == "T" || gate.first == "T*" ||
                gate.first == "P" || gate.first == "P*" ||
                (gate.first == "Z" && gate.second.size() == 1)) {
            int first_input = name_map[*(gate.second.begin())];
            insert_phase(gate_lookup[gate.first], wires[first_input], phase_expts);
        } else if (gate.first == "Z" && gate.second.size() == 3) {
            list<string>::const_iterator tmp_it = gate.second.begin();
            // The three inputs
            int a = name_map[*(tmp_it++)];
            int b = name_map[*(tmp_it++)];
            int c = name_map[*tmp_it];
            insert_phase(1, wires[a], phase_expts);
            insert_phase(1, wires[b], phase_expts);
            insert_phase(1, wires[c], phase_expts);
            insert_phase(7, wires[a] ^ wires[b], phase_expts);
            insert_phase(7, wires[a] ^ wires[c], phase_expts);
            insert_phase(7, wires[b] ^ wires[c], phase_expts);
            insert_phase(1, wires[a] ^ wires[b] ^ wires[c], phase_expts);
        } else if (gate.first == "H") {
            // This WILL confuse you later on you idiot
            //   You zero the "destroyed" qubit, compute the rank, then replace the
            //   value with each of the phase exponents to see if the rank increases
            //   i.e. the system is inconsistent. This is so you don't have to make
            //   a new matrix -- i.e. instead of preparing the new value and computing
            //   rank, then adding each phase exponent and checking the rank you do it
            //   in place
            Hadamard new_h;
            new_h.qubit = name_map[*(gate.second.begin())];
            new_h.prep  = val_max++;
            new_h.wires = wires;

            // Check previous exponents to see if they're inconsistent
            wires[new_h.qubit].reset();
            int rank = compute_rank(n + m, n + h, wires);
            for (const auto& expt : phase_expts) {
                if (expt.second != 0) {
                    wires[new_h.qubit] = expt.first;
                    if (compute_rank(n + m, n + h, wires) > rank) {
                        new_h.in.insert(expt.first);
                    }
                }
            }

            // Done creating the new hadamard
            hadamards.push_back(new_h);

            // Prepare the new value
            wires[new_h.qubit].reset();
            wires[new_h.qubit].set(new_h.prep);

            // Give this new value a name
            val_map[new_h.prep] = name_max;
            names[name_max] = names[new_h.qubit];
            names[name_max++].append(to_string(new_h.prep));

        } else {
            cout << "ERROR: not a {H, CNOT, X, Y, Z, P, T} circuit" << endl;
            cout << "on the " << gate_index << "th gate_index" << endl;
            cout << "gate.first is: '" << gate.first << "'" << endl;
            cout << "gate.second is: [";
            for (const auto& i : gate.second) {
                cout << i << ", ";
            }
            cout << "]" << endl;
            throw logic_error("Can't parse circuit");

            phase_expts.clear();
        }
    }

}

void character::output(ostream& out) const {
  bool flag;

  out << "U|";
  for (int i = 0; i < (n + m); i++) {
    if (i != 0)  out << " ";
    if (zero[i]) out << "()";
    else out << names[i];
  }

  out << "> --> w^(";

  // Print the phase exponents
  for (auto it = phase_expts.begin(); it != phase_expts.end(); it++) {
    if (it != phase_expts.begin()) out << "+";
    out << (int)(it->second) << "*";
    if (it->first.is_negated()) out << "~";
    for (int i = 0; i < (n + h); i++) {
      if (it->first.test(i)) out << names.at(val_map.at(i));
    }
  }
  out << ")|";

  // Print the output functions
  if (this->outputs.size() > 0) {
      for (int i = 0; i < (n + m); i++) {
          flag = false;
          out << "(";
          if (outputs[i].is_negated()) out << "~";
          for (int j = 0; j < (n + h); j++) {
              if (outputs[i].test(j)) {
                  if (flag) out << " ";
                  out << names.at(val_map.at(j));
                  flag = true;
              }
          }
          out << ")";
      }
  } else { out << "{}"; }
  out << ">\n";

  // Print the Hadamards:
  for (const auto& hadamard : hadamards) {
    out << "H:" << names[hadamard.qubit] << "-->" << hadamard.prep << "\n";
  }
}

void character::print_outputs() const {
    for(size_t i = 0; i < this->outputs.size(); i++) {
        cout << "[" << i << "] : size(" << this->outputs.size() << ")";
        cout << "{";
        for (size_t j = 0; j < this->outputs[i].size(); j++) {
            cout << this->outputs[i].test(j);
        }
        cout << "}" << endl;
    }
}

void insert_phase (unsigned char c, xor_func f, vector<pair<exponent_val, xor_func>> & phases) {
  bool flg = false;
  for (auto it = phases.begin(); it != phases.end() && !flg; it++) {
    if (it->second == f) {
      it->first = (it->first + c) % 8;
      flg = true;
    }
  }
  if (!flg) {
    phases.push_back(make_pair(c, xor_func(f)));
  }
}

void insert_phase (unsigned char c, xor_func f, map<xor_func, unsigned char> & phases) {
  phases[f] += c;
  phases[f] %= 8;
  /*  If we're pressed for space?
  if(phases[f] == 0) {
      phases.erase(f);
  }
  */
}


void character::add_ancillae(int num) {
  int num_qubits = this->n + this->m + num;
  stringstream ss;

  int new_m = this->m + num;
  vector<string> new_names(num_qubits);
  vector<bool> new_zero(num_qubits);
  vector<xor_func> new_out;
  new_out.reserve(num_qubits);

  for (auto it = hadamards.begin(); it != hadamards.end(); it++) {
    for (int i = n + m; i < num_qubits; i++) {
      it->wires.push_back(xor_func(n + h));
    }
  }

  for (int i = 0; i < num_qubits; i++) {
    if (i < (n + m)) {
      new_names[i] = names[i];
      new_zero[i]  = zero[i];
      new_out.push_back(outputs[i]);
    } else {
      ss.str("");
      ss << "__anc" << i - (n + m);
      new_names[i] = ss.str();
      new_zero[i] = true;
      new_out.push_back(xor_func(n + h));
    }
  }

  this->names = new_names;
  this->zero = new_zero;
  this->outputs = new_out;
  this->m = new_m;
}

//---------------------------- Synthesis

dotqc character::synthesize() {
  partitioning floats[2], frozen[2];
  dotqc ret{};
  // TODO investigate mask, is +1 needed?
  xor_func mask(n + h);      // Tells us what values we have prepared
  vector<xor_func> wires;        // Current state of the wires
  wires.reserve(n+m);
  list<xor_func> remaining[2];          // Which terms we still have to partition
  int dim = n;
  int tdepth = 0, h_count = 1, applied = 0;
  ind_oracle oracle(n + m, dim, n + h);

  assert(this->outputs.size() > 0);

  // initialize some stuff
  ret.n = n;
  ret.m = m;
  mask.negate();
  for (int i = 0, j = 0; i < n + m; i++) {
    ret.names.push_back(names[i]);
    ret.zero[names[i]] = zero[i];
    wires.push_back(xor_func(n + h));
    if (!zero[i]) {
      wires[i].set(j);
      mask.set(j++);
    }
  }

  // initialize the remaining list
  for (const auto& xpt : phase_expts) {
    if (xpt.second % 2 == 1) remaining[0].push_back(xpt.first);
    else if (xpt.second != 0) remaining[1].push_back(xpt.first);
  }

  // create an initial partition
  // cerr << "Adding new functions to the partition... " << flush;
  for (int j = 0; j < 2; j++) {
    for (auto it = remaining[j].begin(); it != remaining[j].end();) {
      if (mask.contains(*it)) {
        add_to_partition(floats[j], *it, oracle);
        it = remaining[j].erase(it);
      } else it++;
    }
  }
  if (disp_log) cerr << "  " << phase_expts.size() - (remaining[0].size() + remaining[1].size())
    << "/" << phase_expts.size() << " phase rotations partitioned\n" << flush;

  for (auto it = hadamards.begin(); it != hadamards.end(); it++, h_count++) {
    // 1. freeze partitions that are not disjoint from the hadamard input
    // 2. construct CNOT+T circuit
    // 3. apply the hadamard gate
    // 4. add new functions to the partition
    if (disp_log) cerr << "  Hadamard " << h_count << "/" << hadamards.size() << "\n" << flush;

    // determine frozen partitions
    for (int j = 0; j < 2; j++) {
      frozen[j] = freeze_partitions(floats[j], it->in);
      applied += num_elts(frozen[j]);
    }

    // Construct {CNOT, T} subcircuit for the frozen partitions
    {
        auto tmp = construct_circuit(phase_expts, frozen[0], wires, wires, n + m, n + h, names);
        ret.circ.splice(ret.circ.end(), tmp);
        tmp = construct_circuit(phase_expts, frozen[1], wires, it->wires, n + m, n + h, names);
        ret.circ.splice(ret.circ.end(), tmp);
    }
    for (int i = 0; i < n + m; i++) {
      wires[i] = it->wires[i];
    }
    if (disp_log) cerr << "    " << applied << "/" << phase_expts.size() << " phase rotations applied\n" << flush;

    // Apply Hadamard gate
    ret.circ.push_back(make_pair("H", list<string>(1, names[it->qubit])));
    wires[it->qubit].reset();
    wires[it->qubit].set(it->prep);
    mask.set(it->prep);

    // Check for increases in dimension
    int rank = compute_rank(n + m, n + h, &wires[0]);
    if (rank > dim) {
      if (disp_log) cerr << "    Dimension increased to " << rank << ", fixing partitions...\n" << flush;
      dim = rank;
      oracle.set_dim(dim);
      repartition(floats[0], oracle);
      repartition(floats[1], oracle);
    }

    // Add new functions to the partition
    for (int j = 0; j < 2; j++) {
      for (auto it = remaining[j].begin(); it != remaining[j].end();) {
        if (mask.contains(*it)) {
          add_to_partition(floats[j], *it, oracle);
          it = remaining[j].erase(it);
        } else it++;
      }
    }
    if (disp_log) {
        cerr << "    "
            << phase_expts.size() - (remaining[0].size() + remaining[1].size())
            << "/" << phase_expts.size() << " phase rotations partitioned\n"
            << flush;
    }
  }

  applied += num_elts(floats[0]) + num_elts(floats[1]);
  // Construct the final {CNOT, T} subcircuit
  {
      auto tmp = construct_circuit(this->phase_expts, floats[0],
              wires, wires, n + m, n + h, this->names);
      ret.circ.splice(ret.circ.end(), tmp);
      cout << "Printing outputs" << endl;
      this->print_outputs();
      cout << "Printed outputs" << endl;
      tmp = construct_circuit(this->phase_expts, floats[1],
              wires, this->outputs, n + m, n + h, this->names);
      ret.circ.splice(ret.circ.end(), tmp);
  }
  if (disp_log) cerr << "  " << applied << "/" << phase_expts.size() << " phase rotations applied\n" << flush;

  return ret;
}

dotqc character::synthesize_unbounded() {
  partitioning floats[2], frozen[2];
  dotqc ret;
  xor_func mask(n + h);      // Tells us what values we have prepared
  vector<xor_func> wires; // Current state of the wires
  wires.reserve(n + m);
  list<xor_func> remaining[2];          // Which terms we still have to partition
  int dim = n, tmp1, tmp2, tdepth = 0, h_count = 1, applied = 0, j;
  ind_oracle oracle(n + m, dim, n + h);
  gatelist circ;
  list<Hadamard>::iterator it;

  // initialize the wires
  for (int i = 0, j = 0; i < n + m; i++) {
    wires.push_back(xor_func(n + h));
    if (!zero[i]) {
      wires[i].set(j);
      mask.set(j++);
    }
  }

  // initialize the remaining list
  for (const auto& xpt : phase_expts) {
    if (xpt.second % 2 == 1) remaining[0].push_back(xpt.first);
    else if (xpt.second != 0) remaining[1].push_back(xpt.first);
  }

  // create an initial partition
  // cerr << "Adding new functions to the partition... " << flush;
  for (int j = 0; j < 2; j++) {
    for (auto it = remaining[j].begin(); it != remaining[j].end();) {
      if (mask.contains(*it)) {
        add_to_partition(floats[j], *it, oracle);
        it = remaining[j].erase(it);
      } else it++;
    }
  }
  if (disp_log) cerr << "  " << phase_expts.size() - (remaining[0].size() + remaining[1].size())
    << "/" << phase_expts.size() << " phase rotations partitioned\n" << flush;

  for (auto it = hadamards.begin(); it != hadamards.end(); it++, h_count++) {
    // 1. freeze partitions that are not disjoint from the hadamard input
    // 2. construct CNOT+T circuit
    // 3. apply the hadamard gate
    // 4. add new functions to the partition
    if (disp_log) cerr << "  Hadamard " << h_count << "/" << hadamards.size() << "\n" << flush;

    tmp1 = compute_rank(n + m, n + h, wires);
    // determine frozen partitions
    for (j = 0; j < 2; j++) {
      frozen[j] = freeze_partitions(floats[j], it->in);
      applied += num_elts(frozen[j]);
      // determine if we need to add ancillae
      if (frozen[j].size() != 0) {
          // TODO audit
        tmp2 = compute_rank(*(frozen[j].begin()));
        int etc = ((tmp1 - tmp2 < 0)?tmp1:tmp1 - tmp2) + num_elts(frozen[j]) - n - m;
        if (etc > 0) {
          for (int i = n + m; i < n + m + etc; i++) {
            wires.push_back(xor_func(n + h));
          }
          this->add_ancillae(etc);
        }
      }
    }

    if (disp_log) cerr << "    Synthesizing T-layer\n" << flush;
    // Construct {CNOT, T} subcircuit for the frozen partitions
    ret.circ.splice(ret.circ.end(),
        construct_circuit(phase_expts, frozen[0], wires, wires, n + m, n + h, names));
    ret.circ.splice(ret.circ.end(),
        construct_circuit(phase_expts, frozen[1], wires, it->wires, n + m, n + h, names));
    for (int i = 0; i < n + m; i++) {
      wires[i] = it->wires[i];
    }
    if (disp_log) cerr << "    " << applied << "/" << phase_expts.size() << " phase rotations applied\n" << flush;

    // Apply Hadamard gate
    ret.circ.push_back(make_pair("H", list<string>(1, names[it->qubit])));
    wires[it->qubit].reset();
    wires[it->qubit].set(it->prep);
    mask.set(it->prep);

    // Add new functions to the partition
    for (int j = 0; j < 2; j++) {
      for (auto it = remaining[j].begin(); it != remaining[j].end();) {
        if (mask.contains(*it)) {
          if (floats[j].size() == 0) {
              floats[j].push_back(set<xor_func>());
          }
          (floats[j].begin())->insert(*it);
          it = remaining[j].erase(it);
        } else it++;
      }
    }
    if (disp_log) cerr << "    " << phase_expts.size() - (remaining[0].size() + remaining[1].size())
      << "/" << phase_expts.size() << " phase rotations partitioned\n" << flush;
  }

  applied += num_elts(floats[0]) + num_elts(floats[1]);
  // Construct the final {CNOT, T} subcircuit

  // determine if we need to add ancillae
  tmp1 = compute_rank(n + m, n + h, wires);
  for (j = 0; j < 2; j++) {
    if (floats[j].size() != 0) {
      for (j = 0; j < 2; j++) {
          // TODO audit
        tmp2 = compute_rank(*(floats[j].begin()));
        int etc = tmp1 - tmp2 + num_elts(floats[j]) - n - m;
        if (etc > 0) {
          if (disp_log) cerr << "    " << "Adding " << etc << " ancilla(e)\n" << flush;
          for (int i = n + m; i < n + m + etc; i++) {
            wires.push_back(xor_func(n + h));
          }
          this->add_ancillae(etc);
        }
      }
    }
  }

  ret.circ.splice(ret.circ.end(),
      construct_circuit(phase_expts, floats[0],
          wires, wires, n + m, n + h, names));
  ret.circ.splice(ret.circ.end(),
      construct_circuit(phase_expts, floats[1],
          wires, this->outputs, n + m, n + h, names));
  if (disp_log) cerr << "  " << applied << "/" << phase_expts.size() << " phase rotations applied\n" << flush;

  ret.n = n;
  ret.m = m;
  for (int i = 0, j = 0; i < n + m; i++) {
    ret.names.push_back(names[i]);
    ret.zero[names[i]] = zero[i];
  }
  return ret;
}

// TODO bring back metacircuit
/*
//-------------------------------- old {CNOT, T} version code. Still used for the "no hadamards" option

void metacircuit::partition_dotqc(dotqc & input) {
  circuit_type current = UNKNOWN;

  n = input.n;
  m = input.m;
  circuit_list.clear();
  names = input.names;
  zero  = input.zero;

  dotqc acc;
  map<string, bool> zero_acc = input.zero;

  acc.zero = zero_acc;
  for (auto it = input.circ.begin(); it != input.circ.end(); it++) {
    if ((it->first == "T"   && it->second.size() == 1) ||
        (it->first == "T*"  && it->second.size() == 1) ||
        (it->first == "P"   && it->second.size() == 1) ||
        (it->first == "P*"  && it->second.size() == 1) ||
        (it->first == "X"   && it->second.size() == 1) ||
        (it->first == "Y"   && it->second.size() == 1) ||
        (it->first == "Z"   && (it->second.size() == 1 || it->second.size() == 3)) ||
        (it->first == "tof" && (it->second.size() == 1 || it->second.size() == 2))) {
      if (current == UNKNOWN) {
        current = CNOTT;
      } else if (current != CNOTT) {
        acc.n = acc.m = 0;
        for (auto ti = acc.zero.begin(); ti != acc.zero.end(); ti++) {
          if (ti->second) {
            acc.m += 1;
            if (find(acc.names.begin(), acc.names.end(), ti->first) == acc.names.end()) {
                acc.names.push_back(ti->first);
            }
          }
        }
        acc.n = acc.names.size() - acc.m;
        circuit_list.push_back(make_pair(current, acc));

        acc.clear();
        acc.zero = zero_acc;
        current = CNOTT;
      }
    } else {
      if (current == UNKNOWN) {
        current = OTHER;
      } else if (current != OTHER) {
        acc.n = acc.m = 0;
        for (auto ti = acc.zero.begin(); ti != acc.zero.end(); ti++) {
          if (ti->second) {
            acc.m += 1;
            if (find(acc.names.begin(), acc.names.end(), ti->first) == acc.names.end()) {
                acc.names.push_back(ti->first);
            }
          }
        }
        acc.n = acc.names.size() - acc.m;
        circuit_list.push_back(make_pair(current, acc));

        acc.clear();
        acc.zero = zero_acc;
        current = OTHER;
      }
    }
    for (auto iti = it->second.begin(); iti != it->second.end(); iti++) {
      zero_acc[*iti] = false;
    }
    acc.append(*it);
  }

  acc.n = acc.m = 0;
  for (auto ti = acc.zero.begin(); ti != acc.zero.end(); ti++) {
    if (ti->second) {
      acc.m += 1;
      if (find(acc.names.begin(), acc.names.end(), ti->first) == acc.names.end()) {
          acc.names.push_back(ti->first);
      }
    }
  }
  acc.n = acc.names.size() - acc.m;
  circuit_list.push_back(make_pair(current, acc));
}

void metacircuit::output(ostream& out) {
  for (auto it = circuit_list.begin(); it != circuit_list.end(); it++) {
    if (it->first == CNOTT) {
      character tmp{it->second};
      out << "CNOT, T circuit: " << tmp.n << " " << tmp.m << "\n";
      tmp.output(out);
    }
    else out << "Other: " << it->second.n << " " << it->second.m << "\n";
    it->second.output(out);
    out << "\n";
  }
}

dotqc metacircuit::to_dotqc() {
  dotqc ret;
  ret.n = n;
  ret.m = m;
  ret.names = names;
  ret.zero = zero;

  for (auto it = circuit_list.begin(); it != circuit_list.end(); it++) {
    for (auto ti = it->second.circ.begin(); ti != it->second.circ.end(); ti++) {
      ret.circ.push_back(*ti);
    }
  }

  return ret;
}

void metacircuit::optimize() {
  for (auto it =circuit_list.begin(); it != circuit_list.end(); it++) {
    if (it->first == CNOTT) {
      character tmp{it->second};
      it->second = tmp.synthesize();
    }
  }
}
*/
