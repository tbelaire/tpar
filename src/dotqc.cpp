#include "dotqc.h"
#include "util.h"
//----------------------------------------- DOTQC stuff

using namespace std;

void ignore_white(istream& in) {
  while (in.peek() == ' ' || in.peek() == ';') in.ignore();
}

void dotqc::input(istream& in) {
  int i, j;
  string buf, tmp;
  list<string> namelist;
  n = 0;

  // Inputs
  while (buf != ".v") in >> buf;
  ignore_white(in);
  while(in.peek() != '\n' && in.peek() != '\r') {
    in >> buf;
    names.push_back(buf);
    zero[buf] = 1;
    ignore_white(in);
  }

  // Primary inputs
  while (buf != ".i") in >> buf;
  ignore_white(in);
  while (in.peek() != '\n' && in.peek() != '\r') {
    n++;
    in >> buf;
    zero[buf] = 0;
    ignore_white(in);
  }

  m = names.size() - n;

  // Circuit
  while (buf != "BEGIN") in >> buf;
  in >> tmp;
  while (tmp != "END") {
    namelist.clear();
    // Build up a list of the applied qubits
    ignore_white(in);
    while (in.peek() != '\n' && in.peek() != '\r' && in.peek() != ';') {
      in >> buf;
      int pos = buf.find(';');
      if (pos != string::npos) {
        for (int i = buf.length() - 1; i > pos; i--) {
          in.putback(buf[i]);
        }
        in.putback('\n');
        buf.erase(pos, buf.length() - pos);
      }
      if (find(names.begin(), names.end(), buf) == names.end()) {
        cout << "ERROR: no such qubit \"" << buf << "\"\n" << flush;
        exit(1);
      } else {
        namelist.push_back(buf);
      }
      ignore_white(in);
    }
    if (tmp == "TOF") tmp = "tof";
    circ.push_back(make_pair(tmp, namelist));
    in >> tmp;
  }
}

void dotqc::output(ostream& out) const {
  // Inputs
  out << ".v";
  for (auto name_it = names.begin(); name_it != names.end(); name_it++) {
    out << " " << *name_it;
  }

  // Primary inputs
  out << "\n.i";
  for (auto name_it = names.begin(); name_it != names.end(); name_it++) {
    if (zero.at( *name_it ) == 0) out << " " << *name_it;
  }

  // Outputs
  out << "\n.o";
  for (auto name_it = names.begin(); name_it != names.end(); name_it++) {
    out << " " << *name_it;
  }

  // Circuit
  out << "\n\nBEGIN\n";
  for (auto it = circ.begin(); it != circ.end(); it++) {
    out << it->first;
    for (auto ti = (it->second).begin(); ti != (it->second).end(); ti++) {
      out << " " << *ti;
    }
    out << "\n";
  }
  out << "END\n";
}

int max_depth(const map<string, int> & depths, const list<string> & names) {
  int max = 0;

  for (const auto& name : names) {
      if (depths.at(name) > max) {
          max = depths.at(name);
      }
  }

  return max;
}

// Compute T-depth
int dotqc::count_t_depth() const {
  map<string, int> current_t_depth;
  int d;

  for (const auto& name : names) {
    current_t_depth[name] = 0;
  }

  for (auto it = circ.rbegin(); it != circ.rend(); it++) {
    d = max_depth(current_t_depth, it->second);
    if ((it->first == "T") || (it->first == "T*")) {
      d = d + 1;
    } else if ((it->first == "Z") && (it->second.size() >= 3)) {
      d = d + 3;
    }
    for (const auto wire : it->second) {
      current_t_depth[wire] = d;
    }
  }

  return max_depth(current_t_depth, names);
}

// Gather statistics and print
void dotqc::print_stats() const {
  int H = 0;
  int cnot = 0;
  int X = 0;
  int T = 0;
  int P = 0;
  int Z = 0;
  int tdepth = 0;
  bool tlayer = false;
  set<string> qubits;

  for (const auto& gate : this->circ) {
    // Add each on the inputs to the set of used qubits
    for (const string& qubit : gate.second) {
        qubits.insert(qubit);
    }
    if (gate.first == "T" || gate.first == "T*") {
      T++;
      if (!tlayer) {
        tlayer = true;
        tdepth++;
      }
    } else if (gate.first == "P" || gate.first == "P*") {
        P++;
    } else if (gate.first == "Z" && gate.second.size() == 3) {
      tdepth += 3;
      T += 7;
      cnot += 7;
    } else if (gate.first == "Z") {
        Z++;
    } else {
      tlayer = false;
      if (gate.first == "tof" && gate.second.size() == 2) {
          cnot++;
      } else if (gate.first == "tof" || gate.first == "X") {
          X++;
      } else if (gate.first == "H") {
          H++;
      }
    }
  }

  cout << "#   qubits: " << names.size() << "\n";
  cout << "#   qubits used: " << qubits.size() << "\n";
  cout << "#   H: " << H << "\n";
  cout << "#   cnot: " << cnot << "\n";
  cout << "#   X: " << X << "\n";
  cout << "#   T: " << T << "\n";
  cout << "#   P: " << P << "\n";
  cout << "#   Z: " << Z << "\n";
  cout << "#   tdepth (by partitions): " << tdepth << "\n";
  cout << "#   tdepth (by critical paths): " << count_t_depth() << "\n";

}

// Count the Hadamard gates
int dotqc::count_h() const {
  int ret = 0;

  for (const auto& gate : this->circ) {
    if (gate.first == "H") ret++;
  }

  return ret;
}

void dotqc::append(pair<string, list<string> > gate) {
  circ.push_back(gate);

  for (const string& wire : gate.second) {
    // If it's already not in, add it
    if (find(names.begin(), names.end(), wire) == names.end()) {
        names.push_back(wire);
    }
  }
}

// Optimizations
// TODO audit this
void dotqc::remove_swaps() {
  gatelist::iterator it, tt, ttt;
  list<string>::iterator iti;
  map<string, string>::iterator pit;
  string q1, q2, q1_map, q2_map, tmp;
  bool flg, found_q1, found_q2;
  int i;
  map<string, string> perm;

  for (it = circ.begin(), i = 0; i < (circ.size() - 3);) {
    flg = false;
    if (it->first == "tof" && it->second.size() == 2) {
      iti = it->second.begin();
      q1 = *(iti++);
      q2 = *iti;

      tt = it;
      tt++;
      ttt = tt;
      ttt++;
      if (tt->first == "tof" && tt->second.size() == 2 && ttt->first == "tof" && ttt->second.size() == 2) {
        flg = true;
        iti = tt->second.begin();
        flg &= *(iti++) == q2;
        flg &= *(iti) == q1;
        iti = ttt->second.begin();
        flg &= *(iti++) == q1;
        flg &= *(iti) == q2;

        if (flg) {
          circ.erase(it);
          circ.erase(tt);
          it = circ.erase(ttt);
          found_q1 = found_q2 = false;
          q1_map = (perm.find(q2) != perm.end()) ? perm[q2] : q2;
          q2_map = (perm.find(q1) != perm.end()) ? perm[q1] : q1;
          perm[q1] = q1_map;
          perm[q2] = q2_map;
        }
      }
    }
    if (!flg) {
      // Apply permutation
      for (iti = it->second.begin(); iti != it->second.end(); iti++) {
        if (perm.find(*iti) != perm.end()) *iti = perm[*iti];
      }
      i++;
      it++;
    }
  }
  for (; i < circ.size(); i++) {
    for (iti = it->second.begin(); iti != it->second.end(); iti++) {
      if (perm.find(*iti) != perm.end()) *iti = perm[*iti];
    }
    it++;
  }

  // fix outputs
  while(!perm.empty()) {
    pit = perm.begin();
    if (pit->first == pit->second) perm.erase(pit);
    else {
      list<string> tmp_list1, tmp_list2;
      q1 = pit->second;
      q2 = perm[pit->second];
      tmp_list1.push_back(q1);
      tmp_list1.push_back(q2);
      tmp_list2.push_back(q2);
      tmp_list2.push_back(q1);
      circ.push_back(make_pair("tof", tmp_list1));
      circ.push_back(make_pair("tof", tmp_list2));
      circ.push_back(make_pair("tof", tmp_list1));

      pit->second = q2;
      perm[q1] = q1;
    }
  }
}


void dotqc::remove_ids() {

  // Maps each gate to it's identity.
  map<string, string> ids;

  ids["tof"] = "tof";
  ids["Z"] = "Z";
  ids["H"] = "H";
  ids["P"] = "P*";
  ids["P*"] = "P";
  ids["T"] = "T*";
  ids["T*"] = "T";

  // Was an operation performed which might open up new chances to optimize.
  bool modified = true;
  const bool logging = false;

  while (modified) {
    modified = false;
    for (auto it = this->circ.begin(); it != this->circ.end(); it++) {
      if(logging) cout << "Looking at gate " << it->first << endl;
      bool flg = false;
      gatelist::iterator ti = it;
      ti++;
      // Peek forward, looking for
      // a gate with the same inputs
      // giving up when we see one with overlapping inputs.
      // Note that we don't care about the order of the inputs apart
      // from the first one.  That's because they are all control lines.
      for (; ti != circ.end() && !flg; ti++) {
        if(logging) cout << "Comparing with gate " << ti->first << endl;
        if(logging) cout << "Which has wires at:";
        for (const auto& w : it->second) {
            if(logging) cout << " " << w;
        }
        if(logging) cout << endl;
        const auto i = list_compare(it->second, ti->second);
        switch (i) {
            case list_compare_result::EQUAL:
                if (ids[it->first] == ti->first
                        && *(it->second.begin()) == *(ti->second.begin()) ) {
                    if(logging) cout << "Inside if" << endl;
                    ti = circ.erase(ti);
                    it = circ.erase(it);
                    modified = true;
                }
                flg = true;
                break;
            case list_compare_result::OVERLAPPED:
                flg = true;
                break;
            case list_compare_result::DISJOINT:
                break;
            default:
                cerr << "Unknown result from list_compare" <<endl;
                break;
        }
      }
    }
  }
  if(logging) cout << "DONE" << endl;
}

bool dotqc::operator==(const dotqc& other) const {
    return    (other.n == n
            && other.m == m
            && other.names == names
            && other.zero == zero
            && other.circ == circ);
}

ostream& operator<<(ostream& out, const dotqc& circuit) {
    circuit.output(out);
    return out;
}

//-------------------------------------- End of DOTQC stuff
