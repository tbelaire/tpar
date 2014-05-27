
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

