#include <sstream>

#include <gtest/gtest.h>

#include "util.h"
#include "dotqc.h"

using namespace std;
// Dotqc testing
TEST(dotqc, superSimpleInput) {
    stringstream input;
    input << ".v 1" << endl;
    input << ".i 1" << endl;
    input << ".o 1" << endl
        << endl;
    input << "BEGIN" << endl
        << endl;
    input << "Z 1" << endl
        << endl;
    input << "END" << endl;

    dotqc input_dotqc;
    input_dotqc.input(input);

    EXPECT_EQ(input_dotqc.n, 1);
    EXPECT_EQ(input_dotqc.m, 0);
    list<string> names{"1"};
    EXPECT_EQ(input_dotqc.names, names);
    map<string, bool> zero{{"1", false}};
    EXPECT_EQ(input_dotqc.zero, zero);
    gatelist circ{{"Z", {"1"} } };
    EXPECT_EQ(input_dotqc.circ, circ);

    dotqc initalized_dotqc {.n = 1, .m = 0,
        .names = {"1"},
        .zero = {{"1", false}},
        .input_wires = {"1"},
        .output_wires = {"1"},
        .circ = {{"Z", {"1"} }},
    };
    EXPECT_EQ(input_dotqc, initalized_dotqc);
}

TEST(dotqc, somewhatSimpleInput) {
    stringstream input;
    input << ".v 1 2 3 4" << endl;
    input << ".i 1 2 3 4" << endl;
    input << ".o 1 2 3 4" << endl
        << endl;
    input << "BEGIN" << endl
        << endl;
    input << "Z 1" << endl;
    input << "Z 2 3" << endl
        << endl;
    input << "END" << endl;

    dotqc input_dotqc;
    input_dotqc.input(input);

    dotqc initalized_dotqc {.n = 4, .m = 0,
        .names = {"1", "2", "3", "4"},
        .zero = {{"1", false}, {"2", false}, {"3", false}, {"4", false}},
        .input_wires = {"1", "2", "3", "4"},
        .output_wires = {"1", "2", "3", "4"},
        .circ = {{"Z", {"1"}}, {"Z", {"2", "3"}}},
    };
    EXPECT_EQ(input_dotqc, initalized_dotqc);
}

TEST(removeIds, zz) {
    dotqc zz {.n = 1, .m = 0,
        .names = {"1"},
        .zero = {{"1", false}},
        .input_wires = {},
        .output_wires = {},
        .circ = {
            {"Z", {"1"} },
            {"Z", {"1"} },
        }
    };
    zz.remove_ids();
    EXPECT_EQ(0, zz.circ.size());
}


TEST(removeIds, czcz) {
    dotqc czcz {.n = 2, .m = 0,
        .names = {"A", "B"},
        .zero = {{"A", false}, {"B", false}},
        .input_wires = {"A", "B"},
        .output_wires = {"A", "B"},
        .circ = {
            {"Z", {"A", "B"} },
            {"Z", {"A", "B"} },
        },
    };
    czcz.remove_ids();
    EXPECT_EQ(0, czcz.circ.size());
}

TEST(removeIds, zazbza) {
    dotqc zazbza {.n = 2, .m = 0,
        .names = {"A", "B"},
        .zero = {{"A", false}, {"B", false}},
        .input_wires = {"A", "B"},
        .output_wires = {"A", "B"},
        .circ = {
            {"Z", {"A"} },
            {"Z", {"B"} },
            {"Z", {"A"} },
        },
    };
    zazbza.remove_ids();
    EXPECT_EQ(1, zazbza.circ.size());
    auto gate = *zazbza.circ.begin();
    EXPECT_EQ("Z", gate.first);
    EXPECT_EQ(1, gate.second.size());
    EXPECT_EQ("B", *(gate.second.begin()));
}

TEST(removeIds, ztz) {
    dotqc ztz {.n = 2, .m = 0,
        .names = {"A", "B"},
        .zero = {{"A", false}, {"B", false}},
        .input_wires = {"A", "B"},
        .output_wires = {"A", "B"},
        .circ = {
            {"Z", {"A"} },
            {"T", {"A"} },
            {"Z", {"A"} },
        },
    };
    dotqc ztz_copy = ztz;
    ztz.remove_ids();
    EXPECT_EQ(3, ztz.circ.size());
    EXPECT_EQ(ztz, ztz_copy);
}

TEST(removeIds, ztz2) {
    dotqc ztz2 {.n = 2, .m = 0,
        .names = {"A", "B"},
        .zero = {{"A", false}, {"B", false}},
        .input_wires = {"A", "B"},
        .output_wires = {"A", "B"},
        .circ = {
            {"Z", {"A"} },
            {"T", {"A", "B"} },
            {"Z", {"A"} },
        },
    };
    dotqc ztz2_copy = ztz2;
    ztz2.remove_ids();
    EXPECT_EQ(3, ztz2.circ.size());
    EXPECT_EQ(ztz2, ztz2_copy);
}

TEST(removeIDs, tt) {
    dotqc tt {.n = 1, .m = 0,
        .names = {"1"},
        .zero = {{"1", false}},
        .input_wires = {"1"},
        .output_wires = {"1"},
        .circ = {
            {"T", {"1"} },
            {"T", {"1"} },
        },
    };
    tt.remove_ids();
    EXPECT_EQ(2, tt.circ.size());
}

TEST(removeIDs, tts) {
    dotqc tt {.n = 1, .m = 0,
        .names = {"1"},
        .zero = {{"1", false}},
        .input_wires = {"1"},
        .output_wires = {"1"},
        .circ = {
            {"T", {"1"} },
            {"T*", {"1"} },
        },
    };
    tt.remove_ids();
    EXPECT_EQ(0, tt.circ.size());
}

TEST(removeIDs, superset) {
    dotqc zz {.n = 3, .m = 0,
        .names = {"A", "B", "C"},
        .zero = {{"A", false}, {"B", false}, {"C", false}},
        .input_wires = {"A", "B", "C"},
        .output_wires = {"A", "B", "C"},
        .circ = {
            {"Z", {"A", "B"} },
            {"Z", {"A", "B", "C"} },
        },
    };
    dotqc zz_copy = zz;
    zz.remove_ids();
    EXPECT_EQ(2, zz.circ.size());
    const bool is_equal = zz == zz_copy;
    EXPECT_EQ(true, is_equal);
    EXPECT_EQ(zz, zz_copy);
    // There was a crash cause .zero wasn't initalized properly.
}

TEST(removeIDs, yny) {
    dotqc yny {.n = 2, .m = 0,
        .names = {"A", "B"},
        .zero = {{"A", false}, {"B", false}},
        .input_wires = {"A", "B"},
        .output_wires = {"A", "B"},
        .circ = {
            {"Y", {"A", "B"} },
            {"Y", {"B", "A"} },
        },
    };
    yny.remove_ids();
    EXPECT_EQ(2, yny.circ.size());
}

TEST(removeIDs, yabc) {
    dotqc yny {.n = 3, .m = 0,
        .names = {"A", "B", "C"},
        .zero = {{"A", false}, {"B", false}, {"C", false}},
        .input_wires = {"A", "B"},
        .output_wires = {"A", "B"},
        .circ = {
            {"Y", {"A", "B", "C"} },
            {"Y", {"A", "C", "B"} },
        },
    };
    yny.remove_ids();
    EXPECT_EQ(2, yny.circ.size());
}

int max_depth(const map<string, int> & depths, const list<string> & names);
TEST(depth, maxDepth) {
    map<string, int> depths { {"Apple", 4}, {"Pear", 5}, {"Pinecone", 3} };
    EXPECT_EQ(5, max_depth(depths, {"Apple", "Pear"}));
    EXPECT_EQ(4, max_depth(depths, {"Apple"}));
    EXPECT_EQ(0, max_depth(depths, {}));
    // Throws an error, that might be a good thing.
    // Could mitigate by checking with map::count()
    /* EXPECT_EQ(0, max_depth(depths, {"BEAR"})); */
}

TEST(depth, countT) {
    dotqc t {.n = 3, .m = 0,
        .names = {"A", "B", "C"},
        .zero = {{"A", false}, {"B", false}, {"C", false}},
        .input_wires = {"A", "B"},
        .output_wires = {"A", "B"},
        .circ = {
            {"T", {"A"} },
            {"Y", {"A", "C", "B"} },
        },
    };
    EXPECT_EQ(1, t.count_t_depth());
}

TEST(depth, countNone) {
    dotqc no_t {.n = 3, .m = 0,
        .names = {"A", "B", "C"},
        .zero = {{"A", false}, {"B", false}, {"C", false}},
        .input_wires = {"A", "B"},
        .output_wires = {"A", "B"},
        .circ = {
            {"Z", {"A"} },
            {"Y", {"A", "C", "B"} },
        },
    };
    EXPECT_EQ(0, no_t.count_t_depth());
}

TEST(depth, countTT) {
    dotqc tt {.n = 3, .m = 0,
        .names = {"A", "B", "C"},
        .zero = {{"A", false}, {"B", false}, {"C", false}},
        .input_wires = {"A", "B"},
        .output_wires = {"A", "B"},
        .circ = {
            {"T", {"A"} },
            {"T", {"A"} },
        },
    };
    EXPECT_EQ(2, tt.count_t_depth());
}

TEST(depth, transitive) {
    dotqc cir {.n = 3, .m = 0,
        .names = {"A", "B", "C"},
        .zero = {{"A", false}, {"B", false}, {"C", false}},
        .input_wires = {"A", "B"},
        .output_wires = {"A", "B"},
        .circ = {
            {"T", {"A"} },
            {"Y", {"A", "C", "B"} },
            {"T", {"C"} },
        },
    };
    EXPECT_EQ(2, cir.count_t_depth());
}

TEST(depth, transitive2) {
    dotqc cir {.n = 3, .m = 0,
        .names = {"A", "B", "C"},
        .zero = {{"A", false}, {"B", false}, {"C", false}},
        .input_wires = {"A", "B"},
        .output_wires = {"A", "B"},
        .circ = {
            {"T", {"A"} },
            {"T", {"C"} },
            {"Y", {"A", "C", "B"} },
            {"T", {"C"} },
        },
    };
    EXPECT_EQ(2, cir.count_t_depth());
}

TEST(depth, countH) {
    dotqc cir {.n = 3, .m = 0,
        .names = {"A", "B", "C"},
        .zero = {{"A", false}, {"B", false}, {"C", false}},
        .input_wires = {"A", "B"},
        .output_wires = {"A", "B"},
        .circ = {
            {"T", {"A"} },
            {"H", {"C"} },
            {"Y", {"A", "C", "B"} },
            {"H", {"C"} },
        },
    };
    EXPECT_EQ(2, cir.count_h());
}

TEST(dotqc, append) {
    dotqc cir {.n = 3, .m = 0,
        .names = {"A", "B", "C"},
        .zero = {{"A", false}, {"B", false}, {"C", false}},
        .input_wires = {"A", "B"},
        .output_wires = {"A", "B"},
        .circ = {
            {"T", {"A"} },
            {"H", {"C"} },
            {"Y", {"A", "C", "B"} },
            {"H", {"C"} },
        },
    };
    cir.append({{"T"}, {"A", "B"}});
    EXPECT_EQ(5, cir.circ.size());
    EXPECT_EQ(3, cir.names.size());
    EXPECT_EQ(2, cir.count_t_depth());
}

