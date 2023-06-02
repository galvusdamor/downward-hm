#include "symbolic_h2.h"
#include "symb_h2_bdds.h"

#include "../plugins/plugin.h"

#include "../task_utils/task_properties.h"
#include "../utils/logging.h"

#include <cuddObj.hh>

#include <cassert>

using namespace std;

namespace symbolic_h2 {
SymbolicH2Heuristic::SymbolicH2Heuristic(const plugins::Options &opts)
    : Heuristic(opts) {

    bdds = make_shared<symbolic_2::SymbolicH2BDDs>(task_proxy);

    bdds->init();

	if (log.is_at_least_normal()) {
        log << "Using symbolic h^2." << endl;
    }
}



int SymbolicH2Heuristic::compute_heuristic(const State &ancestor_state) {
    State state = convert_ancestor_state(ancestor_state);
    if (task_properties::is_goal_state(task_proxy, state)) {
        return 0;
    } else {
        int count = bdds->calculate_heuristic(state);
        if (log.is_at_least_verbose()) {
            log << "h^2 value: " << count << endl;
        }
        return count;
    }
}

class SymbolicH2HeuristicFeature : public plugins::TypedFeature<Evaluator, SymbolicH2Heuristic> {
public:
    SymbolicH2HeuristicFeature() : TypedFeature("symb_h2") {
        document_title("symbolic h^2 heuristic");
        add_option<int>("m", "subset size", "1", plugins::Bounds("1", "infinity"));
        Heuristic::add_options_to_feature(*this);

        document_language_support("action costs", "supported");
        document_language_support("conditional effects", "ignored");
        document_language_support("axioms", "ignored");

        document_property(
            "admissible",
            "yes for tasks without conditional effects or axioms");
        document_property(
            "consistent",
            "yes for tasks without conditional effects or axioms");
        document_property(
            "safe",
            "yes for tasks without conditional effects or axioms");
        document_property("preferred operators", "no");
    }
};


static plugins::FeaturePlugin<SymbolicH2HeuristicFeature> _plugin;
}
