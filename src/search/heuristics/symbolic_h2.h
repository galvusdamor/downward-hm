#ifndef HEURISTICS_SYMBOLIC_H2
#define HEURISTICS_SYMBOLIC_H2

#include "../heuristic.h"
#include "symb_h2_bdds.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace plugins {
class Options;
}


namespace symbolic_h2 {

class SymbolicH2Heuristic : public Heuristic {
    std::shared_ptr<symbolic_2::SymbolicH2BDDs> bdds; // The symbolic BDDs are declared

protected:
    virtual int compute_heuristic(const State &ancestor_state) override;

public:
    explicit SymbolicH2Heuristic(const plugins::Options &opts);

};
}

#endif
