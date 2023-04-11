#ifndef HEURISTICS_SYMBOLIC_HM_H
#define HEURISTICS_SYMBOLIC_HM_H

#include "../heuristic.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace plugins {
class Options;
}

namespace symbolic_hm {

class SymbolicHMHeuristic : public Heuristic {
	int m;

protected:
    virtual int compute_heuristic(const State &ancestor_state) override;

public:
    explicit SymbolicHMHeuristic(const plugins::Options &opts);

};
}

#endif
