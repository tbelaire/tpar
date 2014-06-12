#ifndef XOR_FUNC_H
#define XOR_FUNC_H

#include <boost/dynamic_bitset.hpp>
#include "types.h"

class xor_func {
    private:
        boost::dynamic_bitset<> bitset;
        bool negated;
    public:
        xor_func(const bool neg, const std::initializer_list<int> lst);
        explicit xor_func(const size_t size) : bitset(size), negated(false) {}
        /* xor_func() : bitset(), negated(false) {} */
        xor_func(const bool negated, const boost::dynamic_bitset<> bitset)
            : bitset(bitset), negated(negated) {}

        bool is_negated() const { return negated; }
        void negate() { negated = !negated; }

        unsigned long slice(size_t start, size_t offset);

        xor_func operator^(const xor_func& b) const {
            return xor_func(negated ^ b.negated,
                    bitset ^ b.bitset);
        }

        xor_func& operator^=(const xor_func& b) {
            this->bitset ^= b.bitset;
            this->negated ^= b.negated; // double check this
            return *this;
        }
        bool operator==(const xor_func& b) const {
            return this->bitset == b.bitset
                && this->negated == b.negated;
        }
        bool operator<(const xor_func& b) const {
            if(this->negated == b.negated) {
                return this->bitset < b.bitset;
            } else {
                return this->negated < b.negated;
            }
        }
        friend std::ostream& operator<<(std::ostream& out, const xor_func& f);
        bool contains(const xor_func& b) const {
            // Explicitly don't care about `negated`.
            return ((~ this->bitset) & b.bitset).none();
        }
        // Pass through functions
        size_t size() const { return bitset.size(); }
        bool test(const size_t i) const { return bitset.test(i); }
        void set(const size_t i, const bool val = true) {
            bitset.set(i, val);
        }
        void reset(const size_t i) { bitset.reset(i); }
        void reset() { bitset.reset(); }
        void flip(const size_t i) { bitset.flip(i); }
        bool none() const { return bitset.none(); }
        bool any() const { return bitset.any(); }

        bool operator[](const size_t& i) const {
            return this->bitset[i];
        }
};

std::ostream& operator<<(std::ostream& out, const std::vector<xor_func>& arr);
#endif // XOR_FUNC_H
