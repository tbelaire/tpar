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
    ASSERT_EQ(2, bits.size());
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
    ASSERT_EQ(1, x.size());
    EXPECT_EQ("tof: B, C", stringify_gate(*x.begin()));
}

gatelist swap_com(int a, int b, const vector<string> names);
TEST(components, swap) {
    const gatelist swap_cir = swap_com(1,2, {"A", "B", "C"});
    EXPECT_EQ(3, swap_cir.size());
    auto it = swap_cir.begin();
    const auto a = *it; it++;
    const auto b = *it; it++;
    const auto c = *it; it++;
    EXPECT_EQ("tof: B, C", stringify_gate(a));
    EXPECT_EQ("tof: C, B", stringify_gate(b));
    EXPECT_EQ("tof: B, C", stringify_gate(c));

}

// TODO X instead of tof is ok, right?
gatelist x_com(int a, const vector<string> names);
TEST(components, x) {
    const gatelist x = x_com(2, {"A", "B", "C"});
    EXPECT_EQ(1, x.size());
    EXPECT_EQ("X: C", stringify_gate(*x.begin()));
}

TEST(echelon, upperCallCount) {
    vector<xor_func>arr{
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
    vector<xor_func>arr {
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


TEST(echelon, upperGates) {
    const vector<xor_func>arr {
            {true, {0,1}},
            {false, {1,1}},
        };
    gatelist gates = to_upper_echelon(2,2, arr, {"A", "B"});
    auto g = gates.begin();
    EXPECT_EQ(4, gates.size()); // 1 X, 3 Swap
    EXPECT_EQ("X: A", stringify_gate(*g));
    g++;
    // Swap
    EXPECT_EQ("tof: B, A", stringify_gate(*g));
    g++;
    EXPECT_EQ("tof: A, B", stringify_gate(*g));
    g++;
    EXPECT_EQ("tof: B, A", stringify_gate(*g));
    g++;
}
TEST(DISABLED_echelon, upperMat) {
    vector<xor_func>A {
            {false, {1,0,0}},
            {false, {0,0,1}},
            {false, {0,1,0}},
        };
    vector<xor_func>B {
            {false, {1,0,0}},
            {false, {0,1,0}},
            {false, {0,0,1}},
        };
    to_upper_echelon(3, 3, A, B);
    EXPECT_EQ(true,  B[0][0]);
    EXPECT_EQ(false, B[0][1]);
    EXPECT_EQ(false, B[0][2]);

    EXPECT_EQ(false, B[1][0]);
    EXPECT_EQ(false, B[1][1]);
    EXPECT_EQ(true,  B[1][2]);

    EXPECT_EQ(false, B[2][0]);
    EXPECT_EQ(true,  B[2][1]);
    EXPECT_EQ(false, B[2][2]);
    EXPECT_EQ(1,0);
}
TEST(echelon, lowerMat2x2) {
    vector<xor_func>arr {
            {false, {1,1}},
            {false, {0,1}},
        };
    vector<xor_func>arr2 = arr;
    to_lower_echelon(2, 2, arr, arr2);

    EXPECT_EQ(arr, arr2);
    EXPECT_EQ(true, arr[0][0]);
    EXPECT_EQ(false, arr[0][1]);
    EXPECT_EQ(false, arr[1][0]);
    EXPECT_EQ(true, arr[1][1]);
}
TEST(echelon, lowerMat2x2EmptyCol1) {
    vector<xor_func>arr {
            {false, {0,1,1}},
            {false, {0,0,1}},
        };
    vector<xor_func>arr2 = arr;
    to_lower_echelon(3, 2, arr, arr2);

    EXPECT_EQ(arr, arr2);
    EXPECT_EQ(false, arr[0][0]);
    EXPECT_EQ(true,  arr[0][1]);
    EXPECT_EQ(false, arr[0][2]);

    EXPECT_EQ(false, arr[1][0]);
    EXPECT_EQ(false, arr[1][1]);
    EXPECT_EQ(true,  arr[1][2]);
}
TEST(echelon, lowerMat2x2EmptyCol2) {
    vector<xor_func>arr {
            {false, {1,0,1}},
            {false, {0,0,1}},
        };
    vector<xor_func>arr2 = arr;
    to_lower_echelon(3, 2, arr, arr2);

    EXPECT_EQ(arr, arr2);
    EXPECT_EQ(true,  arr[0][0]);
    EXPECT_EQ(false, arr[0][1]);
    EXPECT_EQ(false, arr[0][2]);

    EXPECT_EQ(false, arr[1][0]);
    EXPECT_EQ(false, arr[1][1]);
    EXPECT_EQ(true,  arr[1][2]);
}
TEST(echelon, lowerMat2x2EmptyCol3) {
    vector<xor_func>arr {
            {false, {1,1,0}},
            {false, {0,1,0}},
        };
    vector<xor_func>arr2 = arr;
    to_lower_echelon(3, 2, arr, arr2);

    EXPECT_EQ(arr, arr2);
    EXPECT_EQ(true,  arr[0][0]);
    EXPECT_EQ(false, arr[0][1]);
    EXPECT_EQ(false, arr[0][2]);

    EXPECT_EQ(false, arr[1][0]);
    EXPECT_EQ(true,  arr[1][1]);
    EXPECT_EQ(false, arr[1][2]);
}
TEST(echelon, lowerMat2x3) {
    vector<xor_func>arr {
            {false, {1,1,1}},
            {false, {0,1,1}},
        };
    vector<xor_func>arr2 = arr;
    to_lower_echelon(3, 2, arr, arr2);

    EXPECT_EQ(arr, arr2);
    EXPECT_EQ(true, arr[0][0]);
    EXPECT_EQ(false, arr[0][1]);
    EXPECT_EQ(false, arr[0][2]);
    EXPECT_EQ(false, arr[1][0]);
    EXPECT_EQ(true, arr[1][1]);
    EXPECT_EQ(true, arr[1][2]);
}
TEST(echelon, lowerMat3x3) {
    vector<xor_func>arr {
            {false, {1,1,1}},
            {false, {0,1,1}},
            {false, {0,0,1}},
        };
    vector<xor_func>arr2 = arr;
    to_lower_echelon(3, 3, arr, arr2);

    EXPECT_EQ(arr, arr2);
    EXPECT_EQ(true,  arr[0][0]);
    EXPECT_EQ(false, arr[0][1]);
    EXPECT_EQ(false, arr[0][2]);

    EXPECT_EQ(false, arr[1][0]);
    EXPECT_EQ(true,  arr[1][1]);
    EXPECT_EQ(false, arr[1][2]);

    EXPECT_EQ(false, arr[2][0]);
    EXPECT_EQ(false, arr[2][1]);
    EXPECT_EQ(true,  arr[2][2]);
}
TEST(echelon, lowerMat3x5) {
    vector<xor_func>arr {
            {false, {1,1,1,1,0}},
            {false, {0,1,1,0,1}},
            {false, {0,0,0,1,0}},
        };
    vector<xor_func>arr2 = arr;
    to_lower_echelon(5, 3, arr, arr2);

    EXPECT_EQ(arr, arr2);
    EXPECT_EQ(true,  arr[0][0]);
    EXPECT_EQ(false, arr[0][1]);
    EXPECT_EQ(false, arr[0][2]);
    EXPECT_EQ(false, arr[0][3]);
    EXPECT_EQ(true,  arr[0][4]);

    EXPECT_EQ(false, arr[1][0]);
    EXPECT_EQ(true,  arr[1][1]);
    EXPECT_EQ(true,  arr[1][2]);
    EXPECT_EQ(false, arr[1][3]);
    EXPECT_EQ(true,  arr[1][4]);

    EXPECT_EQ(false, arr[2][0]);
    EXPECT_EQ(false, arr[2][1]);
    EXPECT_EQ(false, arr[2][2]);
    EXPECT_EQ(true,  arr[2][3]);
    EXPECT_EQ(false, arr[2][4]);

    // And the tag along matrix too
    EXPECT_EQ(true,  arr2[0][0]);
    EXPECT_EQ(false, arr2[0][1]);
    EXPECT_EQ(false, arr2[0][2]);
    EXPECT_EQ(false, arr2[0][3]);
    EXPECT_EQ(true,  arr2[0][4]);

    EXPECT_EQ(false, arr2[1][0]);
    EXPECT_EQ(true,  arr2[1][1]);
    EXPECT_EQ(true,  arr2[1][2]);
    EXPECT_EQ(false, arr2[1][3]);
    EXPECT_EQ(true,  arr2[1][4]);

    EXPECT_EQ(false, arr2[2][0]);
    EXPECT_EQ(false, arr2[2][1]);
    EXPECT_EQ(false, arr2[2][2]);
    EXPECT_EQ(true,  arr2[2][3]);
    EXPECT_EQ(false, arr2[2][4]);
}

// TODO fix_basis
TEST(DISABLED_fixBasis, basic) {
    EXPECT_EQ(1,0);
}
void compose(int num, vector<xor_func>& A, const vector<xor_func>& B);
TEST(compose, basic) {

    vector<xor_func>A {
            {false, {1,0,0}},
            {false, {0,1,0}},
            {false, {0,0,1}},
        };
    vector<xor_func>B {
            {false, {1,0,0}},
            {false, {0,0,1}},
            {false, {0,1,0}},
        };
    compose(3, A, B);
    EXPECT_EQ(true,  A[0][0]);
    EXPECT_EQ(false, A[0][1]);
    EXPECT_EQ(false, A[0][2]);

    EXPECT_EQ(false, A[1][0]);
    EXPECT_EQ(false, A[1][1]);
    EXPECT_EQ(true,  A[1][2]);

    EXPECT_EQ(false, A[2][0]);
    EXPECT_EQ(true,  A[2][1]);
    EXPECT_EQ(false, A[2][2]);
}

gatelist Lwr_CNOT_synth(int n, int m, vector<xor_func>& bits, const vector<string>& names, bool rev);
TEST(LwrCNotSynth, id3x3) {
    auto arr = vector<xor_func>{
            {false, {1,0,0}},
            {false, {0,1,0}},
            {false, {0,0,1}},
        };
    auto gates = Lwr_CNOT_synth(3, 1, arr, {"A", "B", "C"}, false);
    EXPECT_EQ(0, gates.size());
}
TEST(LwrCNotSynth, revId3x3) {
    auto arr = vector<xor_func>{
            {false, {0,0,1}},
            {false, {0,1,0}},
            {false, {1,0,0}},
        };
    auto gates = Lwr_CNOT_synth(3, 1, arr, {"A", "B", "C"}, false);
    // Should be just a swap
    ASSERT_EQ(2, gates.size());
    auto g = gates.begin();
    // Parial swap
    EXPECT_EQ("tof: A, C", stringify_gate(*g));
    g++;
    EXPECT_EQ("tof: C, A", stringify_gate(*g));
    g++;
}
TEST(LwrCNotSynth, 4x4) {
    auto arr = vector<xor_func>{
            {false, {1,1,0,0}},
            {false, {1,0,0,1}},
            {false, {0,1,0,0}},
            {false, {1,1,1,0}},
        };
    auto gates = Lwr_CNOT_synth(4, 2, arr, {"A", "B", "C", "D"}, false);
    for(const auto& g : gates){
        cout << stringify_gate(g) << endl;
    }
    // Fix last 2
    // 4 -> 1 + 4
    auto g = gates.begin();
    // Fix bottom row
    EXPECT_EQ("tof: D, A", stringify_gate(*g));
    g++;
    // Fix box 0,0
    EXPECT_EQ("tof: B, A", stringify_gate(*g));
    g++;
    EXPECT_EQ("tof: C, B", stringify_gate(*g));
    g++;

}
TEST(LwrCNotSynth, 6x6) {
    auto arr = vector<xor_func>{
            {false, {1,1,0,0,0,0}},
            {false, {1,0,0,1,1,0}},
            {false, {0,1,0,0,1,0}},
            {false, {1,1,1,1,1,1}},
            {false, {1,1,0,1,1,1}},
            {false, {0,0,1,1,1,0}},
        };
    auto gates = Lwr_CNOT_synth(6, 2, arr,
            {"A", "B", "C", "D", "E", "F"}, false);
    for(const auto& g : gates){
        cout << stringify_gate(g) << endl;
    }
    // Fix last 2
    // 4 -> 1 + 4
    auto g = gates.begin();
    // A -> D
    EXPECT_EQ("tof: D, A", stringify_gate(*g));
    g++;
    // A -> E
    EXPECT_EQ("tof: E, A", stringify_gate(*g));
    g++;
    // A -> B
    EXPECT_EQ("tof: B, A", stringify_gate(*g));
    g++;
    // B -> C
    EXPECT_EQ("tof: C, B", stringify_gate(*g));
    g++;
    // C -> E
    EXPECT_EQ("tof: E, C", stringify_gate(*g));
    g++;
    // D -> F
    EXPECT_EQ("tof: F, D", stringify_gate(*g));
    g++;
    // D -> C
    EXPECT_EQ("tof: C, D", stringify_gate(*g));
    g++;
    // C -> D
    EXPECT_EQ("tof: D, C", stringify_gate(*g));
    g++;
    EXPECT_EQ(g, gates.end());
}


gatelist CNOT_synth(int n,
        vector<xor_func>& bits,
        const vector<string> names);
TEST(DISABLED_CNotSynth, I3x3) {
    vector<xor_func> arr{
        {false, {1,0,0}},
        {false, {0,1,0}},
        {false, {0,0,1}},
    };
    gatelist gates = CNOT_synth(3, arr, {"A", "B", "C"});
    EXPECT_EQ(0, gates.size());
}
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

TEST(stringifyGate, basic) {
    gatelist g = {
        {"tof", {"A", "B", "C"}},
        {"X", {"A", "D", "boom"}}};
    auto it = g.begin();
    EXPECT_EQ("tof: A, B, C", stringify_gate(*it));
    it++;
    EXPECT_EQ("X: A, D, boom", stringify_gate(*it));
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
