#ifndef LANDMARKS_LANDMARK_H
#define LANDMARKS_LANDMARK_H

#include "../task_proxy.h"

#include <set>

namespace landmarks {
class Landmark {
public:
    Landmark(std::vector<FactPair> &facts, bool disjunctive, bool conjunctive,
             bool is_true_in_goal = false, bool is_derived = false)
        : facts(facts), disjunctive(disjunctive), conjunctive(conjunctive),
          is_true_in_goal(is_true_in_goal), is_derived(is_derived), cost(1) {
    }

    std::vector<FactPair> facts;
    bool disjunctive;
    bool conjunctive;
    bool is_true_in_goal;
    bool is_derived;

    // Cost of achieving the landmark (as determined by the landmark factory)
    int cost;

    std::set<int> first_achievers;
    std::set<int> possible_achievers;

    bool is_true_in_state(const State &state) const;
};
}
#endif
