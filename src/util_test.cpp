
#include <boost/dynamic_bitset.hpp>
#include "gtest/gtest.h"


#include "util.h"

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
