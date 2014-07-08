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

#include "circuit.h"
#include <cstdio>
#include <iomanip>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#ifndef CLOCK_MONOTONIC
#include <sys/time.h>
#define CLOCK_MONOTONIC 0


static int
clock_gettime(int foo, struct timespec *ts)
{
    struct timeval tv;
    (void) foo;

    gettimeofday(&tv, nullptr);
    ts->tv_sec = tv.tv_sec;
    ts->tv_nsec = tv.tv_usec * 1000;
    return (0);
}

#endif /* !CLOCK_MONOTONIC */

using namespace std;


int main(int argc, char *argv[]) {
  struct timespec start, end;
  dotqc circuit, synth;
  bool post_process = true;
  int anc = 0;


  po::options_description desc("Allowed Options");
  desc.add_options()
      ("help", "displays help")
      ("ancillae", po::value<string>(), "Number of Ancillae (Or 'n' or 'unbounded')")
      ("no-post-process", po::value<bool>(&post_process)->implicit_value(false)->default_value(true),
       "Remove identities in a post processing step")
      ("synth", po::value<string>())
      ("verbose,v", "Display additional logging")
      ;

  po::variables_map vm;
  try {
      po::store(po::command_line_parser(argc, argv)
              .options(desc).run(), vm);
      if (vm.count("help")) {
          cout << desc << endl;
          return 0;
      }
      po::notify(vm);
  }
  catch(std::exception& e)
  {
      cout << "Error: " << e.what() << endl;
      cout << desc << endl;
      return 1;
  }

  disp_log = vm.count("verbose");
  cout << "# Logging: " << (disp_log ? "Yes" : "No") << endl;

  if (vm.count("ancillae")) { // TODO unbounded and n
      string ancillae = vm["ancillae"].as<string>();
      if(ancillae == "n") {
          anc = -1;
      } else if (ancillae == "unbounded") {
          anc = -2;
      } else {
          try {
          anc = boost::lexical_cast<int>(ancillae);
          if (anc < 0) {
              cout << "Error: less than 0 ancillae" << endl;
              return 1;
          }
          } catch (const boost::bad_lexical_cast &)
          {
              cout << "Error: Malformed argument to --ancillae" << endl;
              cout << ancillae << endl;
              return 1;
          }
      }
      if(disp_log) {
          cout << "ancillae: " << anc << endl;
      }
  } else {
      if(disp_log) {
          cout << "ancillae: No" << endl;
      }
  }
  if (vm.count("synth")) {
      string syth_opt = vm["synth"].as<string>();
      if (syth_opt == "ADHOC") {
          cout << "Warning, untested syth type" << endl; // TODO
          synth_method = AD_HOC;
      } else if (syth_opt == "GAUSS") {
          cout << "Warning, untested syth type" << endl; // TODO
          synth_method = GAUSS;
      } else if (syth_opt == "PMH") {
          synth_method = PMH;
      } else {
          cout << "Error: Invalid argument to --synth" << endl; // TODO
      }
      if(disp_log) {
          cout << "synth method: " << syth_opt << endl;
      }
  } else {
      if(disp_log) {
          cout << "synth method: No" << endl;
      }
  }

  if (disp_log) cerr << "Reading circuit...\n" << flush;
  circuit.input(cin);
  cout << "# Original circuit\n" << flush;
  circuit.print_stats();
  cout << flush;

  circuit.remove_ids();
  if (disp_log) cerr << "With full character" << endl;
  if (disp_log) cerr << "Parsing circuit...\n" << flush;
  character c{circuit};
  if (disp_log) cerr << "anc is " << anc << endl;
  if (anc == -1) c.add_ancillae(c.n + c.m);
  else if (anc > 0) c.add_ancillae(anc);
  if (disp_log) cerr << "Resynthesizing circuit...\n" << flush;
  clock_gettime(CLOCK_MONOTONIC, &start);
  if (anc == -2) synth = c.synthesize_unbounded();
  else           synth = c.synthesize();
  clock_gettime(CLOCK_MONOTONIC, &end);

  if (post_process) {
      if (disp_log) cerr << "Applying post-processing...\n" << flush;
      synth.remove_swaps();
      synth.remove_ids();
  }
  cout << "# Optimized circuit\n";
  synth.print_stats();
  cout << fixed << setprecision(3);
  cout << "#   Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000)
      - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";
  synth.print();

  return 0;
}
