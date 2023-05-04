#include "symbolic_hm.h"
#include "symb_hm_bdds.h"

#include "../plugins/plugin.h"

#include "../task_utils/task_properties.h"
#include "../utils/logging.h"

#include <cuddObj.hh>

#include <cassert>

using namespace std;

namespace symbolic_hm {
SymbolicHMHeuristic::SymbolicHMHeuristic(const plugins::Options &opts)
    : Heuristic(opts),
      m(opts.get<int>("m")) {

    // get the max number of preconditions
    int max_num_preconditions = 0;
    for (OperatorProxy op : task_proxy.get_operators()) {
        int num_preconditions = 0;
        for (FactProxy pre : op.get_preconditions()) {
            num_preconditions++;
        }
        if (num_preconditions > max_num_preconditions) {
            max_num_preconditions = num_preconditions;
        }
    }
    // get the amount of operators
    int num_operators = task_proxy.get_operators().size();


    bdds = make_shared<symbolic::SymbolicHMBDDs>(task_proxy, m);

    bdds->init();


	if (log.is_at_least_normal()) {
        log << "Using symbolic h^" << m << "." << endl;
    }
}



int SymbolicHMHeuristic::compute_heuristic(const State &ancestor_state) {
    State state = convert_ancestor_state(ancestor_state);
    if (task_properties::is_goal_state(task_proxy, state)) {
        return 0;
    } else {
        int count = bdds->calculate_heuristic(state);
        if (log.is_at_least_verbose()) {
            log << "h^" << m << " value: " << count << endl;
        }
        return count;
	}
}

class SymbolicHMHeuristicFeature : public plugins::TypedFeature<Evaluator, SymbolicHMHeuristic> {
public:
    SymbolicHMHeuristicFeature() : TypedFeature("symb_hm") {
        document_title("symbolic h^m heuristic");

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


static plugins::FeaturePlugin<SymbolicHMHeuristicFeature> _plugin;
}
