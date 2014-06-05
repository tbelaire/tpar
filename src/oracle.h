#ifndef ORACLE_H
#define ORACLE_H
#include <set>
#include <boost/optional.hpp>

#include "xor_func.h"

class ind_oracle {
  private:
    int num;
    int dim;
    int length;
  public:
    ind_oracle() { num = 0; dim = 0; length = 0; }
    ind_oracle(int numin, int dimin, int lengthin) { num = numin; dim = dimin; length = lengthin; }

    void set_dim(int newdim) { dim = newdim; }
    boost::optional<xor_func> retrieve_lin_dep(const std::set<xor_func> & lst) const;

    bool operator()(const std::set<xor_func> & lst) const;
};
#endif // ORACLE_H
