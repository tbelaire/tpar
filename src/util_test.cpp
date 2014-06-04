
#include <boost/dynamic_bitset.hpp>
#include "gtest/gtest.h"


#include "util.h"

using namespace std;


TEST(utilTest, init_maxtrix) {
    xor_func f = init_xor_func({1,0,0,1});
    EXPECT_EQ(true, f.any());
    EXPECT_EQ(false, f.test(1));

    vector<xor_func> bits = init_matrix_transpose({
                {1,0,1,0},
                {0,1,1,0}
            });
    EXPECT_EQ(2, bits.size());
    EXPECT_EQ(true, bits.at(0).test(0));
    EXPECT_EQ(false, bits.at(0).test(1));
    EXPECT_EQ(false, bits.at(1).test(0));
}

TEST(utilTest, compute_rank) {

    EXPECT_EQ(1, compute_rank(2, 4, init_matrix_transpose({
                {1,0,1,0},
                {1,0,1,0},
            })));
    EXPECT_EQ(2, compute_rank(3, 4, init_matrix_transpose({
                {1,0,1,0},
                {1,1,1,0},
                {1,0,1,0},
            })));
    EXPECT_EQ(0, compute_rank(3, 4, init_matrix_transpose({
                {0,0,0,0},
                {0,0,0,0},
                {0,0,0,0},
            })));
    EXPECT_EQ(3, compute_rank(3, 4, init_matrix_transpose({
                {1,1,1,1},
                {1,0,0,0},
                {1,0,1,1},
            })));
    EXPECT_EQ(2, compute_rank(3, 4, init_matrix_transpose({
                {1,1,1,1},
                {1,0,1,0},
                {0,1,0,1},
            })));
}

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

// num = N = |A'|
bool construct_and_test(int dim, int num, initializer_list<initializer_list<int>> lst) {

    const auto arr = init_matrix_transpose(lst);
    set<xor_func> set;
    for(const auto& f_lst : lst){
        set.insert(init_xor_func(f_lst));
    }
    const int length = arr[0].size();
    ind_oracle oracle{num, dim, length};
    return oracle(set);
}

TEST(oracle, simple) {
    /* The first param is the target dimension, dim or n.
     * The second is the maximum number of vectors, num or N.
     * Next is a set of vectors that must be included in every possible
     * solution.  If it's possible to solve, then it's true.
     *
     * Another way of phrasing it is there are b dependent vectors in set,
     * then return num - b >= dim.
     */
    EXPECT_EQ(true,
        construct_and_test(3, 4, {
            {1,0,0,0},
            {0,1,0,0},
            {1,1,0,0}
            }));
    EXPECT_EQ(false,
        construct_and_test(5, 7, {
            {1,0,0,0,0},
            {0,1,0,0,0},
            {1,1,0,0,0},
            {1,1,1,0,0},
            {0,1,1,0,0},
            {0,0,1,0,0},
            }));
    EXPECT_EQ(false,
        construct_and_test(5, 5, {
            {1,0,0,0,0},
            {0,1,0,0,0},
            {1,1,0,0,0},
            }));
}

/*  Needs to be manually inspected.
void print_wires(const vector<xor_func> wires);
TEST(utilTest, printing) {
    cout << "Should see:" <<endl;
    cout << "1111" << endl;
    cout << "1000" << endl;
    cout << "1011" << endl << endl;

    print_wires( init_matrix_transpose({
                {1,1,1,1},
                {1,0,0,0},
                {1,0,1,1},
            }));
}
*/

