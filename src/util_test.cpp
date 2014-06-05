#include <boost/dynamic_bitset.hpp>
#include "gtest/gtest.h"

#include "util.h"

using namespace std;


TEST(utilTest, initMaxtrix) {
    xor_func f = {false, {1,0,0,1}};
    EXPECT_EQ(true, f.any());
    EXPECT_EQ(false, f.test(1));

    vector<xor_func> bits = init_matrix({
                {1,0,1,0},
                {0,1,1,0}
            });
    EXPECT_EQ(2, bits.size());
    EXPECT_EQ(true, bits.at(0).test(0));
    EXPECT_EQ(false, bits.at(0).test(1));
    EXPECT_EQ(false, bits.at(1).test(0));
}

int compute_rank(int m, int n, const std::vector<xor_func> bits);
TEST(computeRank, vectors) {

    EXPECT_EQ(1, compute_rank(2, 4, init_matrix({
                {1,0,1,0},
                {1,0,1,0},
            })));
    EXPECT_EQ(2, compute_rank(3, 4, init_matrix({
                {1,0,1,0},
                {1,1,1,0},
                {1,0,1,0},
            })));
    EXPECT_EQ(0, compute_rank(3, 4, init_matrix({
                {0,0,0,0},
                {0,0,0,0},
                {0,0,0,0},
            })));
    EXPECT_EQ(3, compute_rank(3, 4, init_matrix({
                {1,1,1,1},
                {1,0,0,0},
                {1,0,1,1},
            })));
    EXPECT_EQ(2, compute_rank(3, 4, init_matrix({
                {1,1,1,1},
                {1,0,1,0},
                {0,1,0,1},
            })));
}

// Same test data, just fed in differently.
int compute_rank(int m, int n, const xor_func * bits);
TEST(computeRank, pointers) {
    EXPECT_EQ(1, compute_rank(2, 4, &init_matrix({
                {1,0,1,0},
                {1,0,1,0},
            })[0]));
    EXPECT_EQ(2, compute_rank(3, 4, &init_matrix({
                {1,0,1,0},
                {1,1,1,0},
                {1,0,1,0},
            })[0]));
    EXPECT_EQ(0, compute_rank(3, 4, &init_matrix({
                {0,0,0,0},
                {0,0,0,0},
                {0,0,0,0},
            })[0]));
    EXPECT_EQ(3, compute_rank(3, 4, &init_matrix({
                {1,1,1,1},
                {1,0,0,0},
                {1,0,1,1},
            })[0]));
    EXPECT_EQ(2, compute_rank(3, 4, &init_matrix({
                {1,1,1,1},
                {1,0,1,0},
                {0,1,0,1},
            })[0]));
}
int compute_rank(int n, const exponents_set & expnts, const std::set<xor_func> & lst);
//TEST(computeRank, xptSets) {
    // I have no idea what this is going to do.
    // Looking at it's history, it only uses the set passed in for it's size
    // but it'd make more sense to compute the rank of the subset.
    // Which would involve reading the expoents_set and the set
    // in past versions, where the set is a set of indecies into expnts,
    // but now we can just compute directly.
//}
TEST(computeRank, sets) {
    EXPECT_EQ(3, compute_rank(set<xor_func>{
                {false, {0,0,0,0,1}},
                {false, {0,0,0,1,0}},
                {false, {0,0,0,1,1}},
                {false, {0,0,1,0,1}},
                {false, {0,0,1,1,1}},
                }));
    EXPECT_EQ(0, compute_rank(set<xor_func>{ }));
}

TEST(computeRank, negatedFuncs) {
    EXPECT_EQ(2, compute_rank(vector<xor_func>{
                {false, {0,0,0,0,1}},
                {false, {0,0,0,1,0}},
                {true,  {0,0,0,0,1}},
                {true,  {0,0,0,1,0}},
            }));
    EXPECT_EQ(1, compute_rank(vector<xor_func>{
                {false, {0,0,0,0,1}},
                {true,  {0,0,0,0,1}},
            }));
    EXPECT_EQ(2, compute_rank(vector<xor_func>{
                {false, {0,0,0,0,1}},
                {false, {0,0,0,1,0}},
                {true,  {0,0,0,1,1}},
            }));
}

gatelist xor_com(int a, int b, const vector<string> names);
TEST(components, xor) {
    const gatelist x = xor_com(1,2, {"A", "B", "C"});
    EXPECT_EQ(1, x.size());
    const pair<string, list<string>> gate = *x.begin();
    EXPECT_EQ("tof", gate.first);
    const auto inputs = gate.second;
    EXPECT_EQ(2, inputs.size());
    EXPECT_EQ("B", *(inputs.begin()));
    EXPECT_EQ("C", *(++inputs.begin()));
}

gatelist swap_com(int a, int b, const vector<string> names);
TEST(components, swap) {
    const gatelist swap_cir = swap_com(1,2, {"A", "B", "C"});
    EXPECT_EQ(3, swap_cir.size());
    auto it = swap_cir.begin();
    const auto a = *it; it++;
    const auto b = *it; it++;
    const auto c = *it; it++;

    EXPECT_EQ("tof", a.first);
    EXPECT_EQ("tof", b.first);
    EXPECT_EQ("tof", c.first);
    EXPECT_EQ(2, a.second.size());
    EXPECT_EQ(2, b.second.size());
    EXPECT_EQ(2, c.second.size());

    list<string> inputs{"B", "C"};
    list<string> inputs_r{"C", "B"};
    EXPECT_EQ(inputs, a.second);
    EXPECT_EQ(inputs_r, b.second);
    EXPECT_EQ(inputs, c.second);
}

gatelist x_com(int a, const vector<string> names);
TEST(components, x) {
    const gatelist x = x_com(2, {"A", "B", "C"});
    EXPECT_EQ(1, x.size());
    const pair<string, list<string>> gate = *x.begin();
    EXPECT_EQ("tof", gate.first);
    const auto inputs = gate.second;
    EXPECT_EQ(1, inputs.size());
    EXPECT_EQ("C", *(inputs.begin()));
}

int to_upper_echelon(int m, int n,
        vector<xor_func> bits,
        std::function<void(int)> do_negate,
        std::function<void(int, int)> do_swap,
        std::function<void(int, int)> do_xor);
TEST(echelon, upperCallCount) {
    const vector<xor_func>arr{
            {false, {0,1}},
            {true, {1,0}},
        };
    int num_negates = 0;
    int num_swaps = 0;
    int num_xors = 0;
    const int rank = to_upper_echelon(2,2, arr,
            [&num_negates](int j){
                (void)j;
                num_negates++;
            },
            [&num_swaps](int r1, int r2){
                (void)r1;
                (void)r2;
                num_swaps++;
            },
            [&num_xors](int r1, int r2){
                (void)r1;
                (void)r2;
                num_xors++;
            }
        );
    EXPECT_EQ(2, rank);
    EXPECT_EQ(1, num_swaps);
    EXPECT_EQ(1, num_negates);
    EXPECT_EQ(0, num_xors);
}

TEST(echelon, upperZeroRow) {
    const vector<xor_func>arr {
            {true, {0,1}},
            {true, {1,0}},
            {false, {1,1}},
        };
    int num_swaps = 0;
    int num_negates = 0;
    int num_xors = 0;
    const int rank = to_upper_echelon(3,2, arr,
            [&num_negates](int j){
                (void)j;
                num_negates++;
            },
            [&num_swaps](int r1, int r2){
                (void)r1;
                (void)r2;
                num_swaps++;
            },
            [&num_xors](int r1, int r2){
                (void)r1;
                (void)r2;
                num_xors++;
            }
        );
    EXPECT_EQ(2, rank);
    EXPECT_EQ(1, num_swaps);
    EXPECT_EQ(2, num_negates);
    EXPECT_EQ(2, num_xors);
}


gatelist to_upper_echelon(int m, int n, vector<xor_func> bits, const vector<string> names);
TEST(DISABLED_echelon, upperGates) {
    // TODO
    EXPECT_EQ(1,0);
}
void to_upper_echelon(int m, int n, vector<xor_func> bits, vector<xor_func> mat);
TEST(DISABLED_echelon, upperMat) {
    // TODO
    EXPECT_EQ(1,0);
}
gatelist to_lower_echelon(int m, int n, vector<xor_func> bits, const vector<string> names);
TEST(DISABLED_echelon, lower) {
    // TODO
    EXPECT_EQ(1,0);
}
void to_lower_echelon(int m, int n, vector<xor_func>& bits, vector<xor_func>& mat);

// TODO fix_basis
TEST(DISABLED_fixBasis, basic) {
    EXPECT_EQ(1,0);
}
// TODO compose
TEST(DISABLED_compose, basic) {
    EXPECT_EQ(1,0);
}

// MAYBE test *_CNOT_synth

// Future TODO construct_curcuit



TEST(listCompare, baseline) {
    EXPECT_EQ(list_compare_result::EQUAL,
              list_compare({"A", "B", "C"}, {"A", "B", "C"}));
    EXPECT_EQ(list_compare_result::OVERLAPPED,
              list_compare({"A", "D", "B"}, {"A", "B", "C"}));
    EXPECT_EQ(list_compare_result::DISJOINT,
              list_compare({"A"}, {"D"}));
    EXPECT_EQ(list_compare_result::DISJOINT,
              list_compare({"A", "B", "C"}, {"D", "E", "F"}));
}

TEST(listCompare, abABC) {
    EXPECT_EQ(list_compare_result::OVERLAPPED,
              list_compare({"A", "B"}, {"A", "B", "C"}));
}

TEST(listCompare, orderingTest) {
    EXPECT_EQ(list_compare_result::EQUAL, list_compare({"A", "C", "B"}, {"A", "B", "C"}));
}
/*  Needs to be manually inspected.
void print_wires(const vector<xor_func> wires);
TEST(utilTest, printing) {
    cout << "Should see:" <<endl;
    cout << "1111" << endl;
    cout << "1000" << endl;
    cout << "1011" << endl << endl;

    print_wires( init_matrix({
                {1,1,1,1},
                {1,0,0,0},
                {1,0,1,1},
            }));
}
*/
