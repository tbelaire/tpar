#include <gtest/gtest.h>

#include <vector>

#include "xor_func.h"

using namespace std;

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

TEST(constructors, copyConstructor) {
    const xor_func a{false, {0, 1, 1, 0, 0}};
    xor_func b{a};

    EXPECT_EQ(a,b);
    EXPECT_EQ(5, b.size());
}

TEST(vectorInit, size2) {
    const xor_func temp{false, {0,1}};
    vector<xor_func> arr{2, temp};
    for(xor_func& f : arr) {
        EXPECT_EQ(2, f.size()) << f;
    }
    for(int i = 0; i < 2; i++) {
        EXPECT_EQ(2, arr[i].size()) << "for i=" << i;
        EXPECT_EQ(true, arr[i].test(1));
    }
}

TEST(vectorInit, size3) {
    const xor_func temp{false, {0,0,0}};
    vector<xor_func> arr( 3, temp );
    for(xor_func& f : arr) {
        EXPECT_EQ(3, f.size()) << f;
    }
    for(int i = 0; i < 3; i++) {
        ASSERT_EQ(3, arr[i].size()) << "for i=" << i;
        arr[i].set(i);
    }
}
