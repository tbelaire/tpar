
#include "gtest/gtest.h"

#include "util.h"
#include "partition.h"
#include "types.h"
#include "matroid.h"
#include "oracle.h"

#include <boost/dynamic_bitset.hpp>

using namespace std;

TEST(partitions, creation) {
    xor_func a{false, {0,0,1,0}};
    xor_func b{false, {0,0,0,1}};

    partitioning p = create(set<xor_func>{a,b});
    ASSERT_EQ(1, p.size());
}

TEST(partitions, n2d2) {
    xor_func a{false, {0,0,1,0}};
    xor_func b{false, {0,0,0,1}};

    ind_oracle oracle(2,2,4);
    partitioning p = create(set<xor_func>{a,b});

    add_to_partition(p, {false, {0,0,1,1}}, oracle);
    EXPECT_EQ(2, p.size());
    add_to_partition(p, {false, {0,1,0,1}}, oracle);
    EXPECT_EQ(2, p.size());
    add_to_partition(p, {false, {0,1,1,0}}, oracle);
    EXPECT_EQ(3, p.size());
    add_to_partition(p, {false, {0,1,1,1}}, oracle);
    EXPECT_EQ(3, p.size());
    /* cout << p << endl; */
}

TEST(partitions, repartition) {
    xor_func a{false, {0,0,1,0}};
    xor_func b{false, {0,0,0,1}};

    ind_oracle oracle(3,2,4);
    partitioning p = create(set<xor_func>{a,b});

    add_to_partition(p, {false, {0,1,0,0}}, oracle);
    add_to_partition(p, {false, {0,1,1,0}}, oracle);
    EXPECT_EQ(2, p.size());
    /* cout << p << endl; */
    oracle.set_dim(2);
    repartition(p, oracle);
    EXPECT_EQ(2, p.size());
    /* cout << p << endl; */
}
