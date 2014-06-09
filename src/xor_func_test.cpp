#include "gtest/gtest.h"

#include "xor_func.h"


TEST(xorFuncSlice, one) {
    xor_func bits{false, {0, 0, 0, 1}};
    EXPECT_EQ(8, bits.slice(0, 4));
    /* EXPECT_EQ(8, slow_slice(bits, 0, 4)); */
}

TEST(xorFuncSlice, x101000) {
    xor_func bits{false, {1, 0, 1, 0, 0, 0}};
    EXPECT_EQ(5, bits.slice(0,3));
    /* EXPECT_EQ(5, slow_slice(bits,0,3)); */
}

TEST(xorFuncSlice, mid) {
    xor_func bits{false, {0, 1, 1, 0, 0}};
    EXPECT_EQ(3, bits.slice(1,3));
    /* EXPECT_EQ(3, slow_slice(bits,1,3)); */
}

