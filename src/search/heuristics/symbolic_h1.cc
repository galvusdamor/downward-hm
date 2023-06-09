#include "symbolic_h1.h"
#include "symb_h1_bdds.h"

#include "../plugins/plugin.h"

#include "../task_utils/task_properties.h"
#include "../utils/logging.h"

#include <cuddObj.hh>

#include <cassert>

using namespace std;

namespace symbolic_h1 {
SymbolicH1Heuristic::SymbolicH1Heuristic(const plugins::Options &opts)
    : Heuristic(opts), var_order(opts.get<int>("var_order")) {

    bdds = make_shared<symbolic_1::SymbolicH1BDDs>(task_proxy, var_order);

    bdds->init();

	if (log.is_at_least_normal()) {
        log << "Using symbolic h^1." << endl;
    }
}



int SymbolicH1Heuristic::compute_heuristic(const State &ancestor_state) {
    State state = convert_ancestor_state(ancestor_state);
    if (task_properties::is_goal_state(task_proxy, state)) {
        return 0;
    } else {
        int count = bdds->calculate_heuristic(state);
        if (log.is_at_least_verbose()) {
            log << "h^1 value: " << count << endl;
        }
        return count;
    }
}

class SymbolicH1HeuristicFeature : public plugins::TypedFeature<Evaluator, SymbolicH1Heuristic> {
public:
    SymbolicH1HeuristicFeature() : TypedFeature("symb_h1") {
        document_title("symbolic h^1 heuristic");
        add_option<int>("var_order", "subset size", "1", plugins::Bounds("1", "2"));
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


static plugins::FeaturePlugin<SymbolicH1HeuristicFeature> _plugin;
}
