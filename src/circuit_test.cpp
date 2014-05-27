#include "circuit.h"

#include "gtest/gtest.h"

#include <sstream>

using namespace std;

// Dotqc testing
TEST(dotqc, superSimpleInput) {
    stringstream input;
    input << ".v 1" << endl;
    input << ".i 1" << endl;
    input << ".o 1" << endl
        << endl;
    input << "BEGIN" << endl
        << endl;
    input << "Z 1" << endl
        << endl;
    input << "END" << endl;

    dotqc input_dotqc;
    input_dotqc.input(input);

    EXPECT_EQ(input_dotqc.n, 1);
    EXPECT_EQ(input_dotqc.m, 0);
    list<string> names{"1"};
    EXPECT_EQ(input_dotqc.names, names);
    map<string, bool> zero{{"1", false}};
    EXPECT_EQ(input_dotqc.zero, zero);
    gatelist circ{{"Z", {"1"} } };
    EXPECT_EQ(input_dotqc.circ, circ);

    dotqc initalized_dotqc {.n = 1, .m = 0,
        .names = {"1"},
        .zero = {{"1", false}},
        .circ = {{"Z", {"1"} }}
    };
    EXPECT_EQ(input_dotqc, initalized_dotqc);
}

TEST(dotqc, somewhatSimpleInput) {
    stringstream input;
    input << ".v 1 2 3 4" << endl;
    input << ".i 1 2 3 4" << endl;
    input << ".o 1 2 3 4" << endl
        << endl;
    input << "BEGIN" << endl
        << endl;
    input << "Z 1" << endl;
    input << "Z 2 3" << endl
        << endl;
    input << "END" << endl;

    dotqc input_dotqc;
    input_dotqc.input(input);

    dotqc initalized_dotqc {.n = 4, .m = 0,
        .names = {"1", "2", "3", "4"},
        .zero = {{"1", false}, {"2", false}, {"3", false}, {"4", false}},
        .circ = {{"Z", {"1"}}, {"Z", {"2", "3"}}}
    };
    EXPECT_EQ(input_dotqc, initalized_dotqc);
}

// character testing
TEST(parsingFromDotQC, initalized) {
    character c;
    dotqc input_dotqc {.n = 1, .m = 0,
        .names = {"1"},
        .zero = {{"1", false}},
        .circ = {{"Z", {"1"} }}
    };
    c.parse_circuit(input_dotqc);
    EXPECT_EQ(1, c.n);
    EXPECT_EQ(0, c.m);
    EXPECT_EQ(0, c.h);
    EXPECT_EQ(1, c.names.size());
    EXPECT_EQ("1", c.names[0]);
    EXPECT_EQ(vector<bool>{0}, c.zero);

    EXPECT_EQ(1, c.val_map.size()); // TODO understand
    /*
    for (const auto& p : c.val_map) {
        cout << p.first << ": " << p.second << endl;
    }
    */
    // We don't have any outputs until we synthesise
    EXPECT_EQ(0, c.outputs.size());

    EXPECT_EQ(1, c.phase_expts.size()); // TODO understand

    // No H gates in this example
    EXPECT_EQ(0, c.hadamards.size());

    c.print_outputs();

    c.print();
}
