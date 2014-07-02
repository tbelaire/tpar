#include "circuit.h"
#include "util.h"

#include <gtest/gtest.h>

#include <sstream>

using namespace std;

// character testing
TEST(parsingFromDotQC, initalized) {
    dotqc input_dotqc {.n = 1, .m = 0,
        .names = {"1"},
        .zero = {{"1", false}},
        .input_wires = {},
        .output_wires = {},
        .circ = {{"Z", {"1"} }}
    };
    character c{input_dotqc};
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
    // That's not true anymore, we have names.size() outputs
    // TODO that might play poorly with synthesize unbounded.
    /* EXPECT_EQ(0, c.outputs.size()); */

    EXPECT_EQ(1, c.phase_expts.size()); // TODO understand

    // No H gates in this example
    EXPECT_EQ(0, c.hadamards.size());

    /* c.print_outputs(); */

    /* c.print(); */
}

//using exponent = pair<u8, xor_func>
void insert_phase (unsigned char c, xor_func f, vector<pair<exponent_val, xor_func>> & phases);
TEST(character, insertPhase) {
    const xor_func f0101 = {false, {0,1,0,1}};
    const xor_func f1010 = {false, {1,0,1,0}};
    const xor_func f1111 = {false, {1,1,1,1}};
    const xor_func f1000 = {false, {1,0,0,0}};
    vector<pair<exponent_val, xor_func>> xpt{ {1, f0101}, {3, f1010}, {1, f1111} };
    insert_phase(1, f0101, xpt);
    EXPECT_EQ(2, xpt[0].first);
    EXPECT_EQ(f0101, xpt[0].second);
    EXPECT_EQ(3, xpt[1].first);
    EXPECT_EQ(1, xpt[2].first);
    insert_phase(6, f0101, xpt);
    EXPECT_EQ(0, xpt[0].first);  // TODO is that expected?
}

void insert_phase (unsigned char c, xor_func f, exponents_set& phases);
TEST(character, insertPhaseMap) {
    const xor_func f0101 = {false, {0,1,0,1}};
    const xor_func f1010 = {false, {1,0,1,0}};
    const xor_func f1111 = {false, {1,1,1,1}};
    const xor_func f1000 = {false, {1,0,0,0}};
    exponents_set xpt{ {f0101, 1}, {f1010, 3}, {f1111, 1} };
    insert_phase(1, f0101, xpt);
    EXPECT_EQ(2, xpt[f0101]);
    EXPECT_EQ(3, xpt[f1010]);
    EXPECT_EQ(1, xpt[f1111]);
    insert_phase(6, f0101, xpt);
    EXPECT_EQ(0, xpt[f0101]);  // TODO is that expected?
}
