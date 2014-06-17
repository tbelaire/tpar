#include "gtest/gtest.h"

#include "util.h"
#include "oracle.h"
using namespace std;
// num = N = |A'|
bool construct_and_test(int dim, int num, initializer_list<initializer_list<int>> lst) {

    // Causing a segfault on a operator= in xor_func
    const auto arr = init_matrix(lst);
    set<xor_func> set;
    for(const auto& f_lst : lst){
        set.insert(xor_func{false, f_lst});
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

// optional<xor_func> ind_oracle::retrieve_lin_dep(const set<xor_func> & lst) const;
TEST(oracle, retrieveLinDep) {
    ind_oracle oracle(3, 3, 4);

    boost::optional<xor_func> res;

    res = oracle.retrieve_lin_dep({
            {false, {0,0,1,1}},
            {false, {0,0,1,0}},
            {false, {0,0,0,1}},
        });
    EXPECT_EQ(true, (bool)res);

    res = oracle.retrieve_lin_dep({
            {false, {0,1,1,1}},
            {false, {0,0,1,0}},
            {false, {0,0,0,1}},
        });
    EXPECT_EQ(false, (bool)res);
}
