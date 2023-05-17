#ifndef HEURISTICS_SYMBOLIC_H1
#define HEURISTICS_SYMBOLIC_H1

#include "../heuristic.h"
#include "symb_h1_bdds.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace plugins {
class Options;
}


namespace symbolic_h1 {

class SymbolicH1Heuristic : public Heuristic {
    std::shared_ptr<symbolic_1::SymbolicH1BDDs> bdds; // The symbolic BDDs are declared

protected:
    virtual int compute_heuristic(const State &ancestor_state) override;

public:
    explicit SymbolicH1Heuristic(const plugins::Options &opts);

};
}

#endif
