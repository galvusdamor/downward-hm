#include "symbolic_hm.h"
#include "transition_relation.h"

#include "../plugins/plugin.h"

#include "../task_utils/task_properties.h"
#include "../utils/logging.h"

#include <cassert>

using namespace std;

namespace symbolic_hm {
SymbolicHMHeuristic::SymbolicHMHeuristic(const plugins::Options &opts)
    : Heuristic(opts),
      m(opts.get<int>("m")) {
	assert(m == 1); // We can for now only handle m=1

	vars = make_shared<symbolic::SymVariables>(opts, task_proxy);
	vars->init();

	symbolic::TransitionRelation tr_test (vars,OperatorID(1),task_proxy);


	if (log.is_at_least_normal()) {
        log << "Using symbolic h^" << m << "." << endl;
    }
}



int SymbolicHMHeuristic::compute_heuristic(const State &ancestor_state) {
    State state = convert_ancestor_state(ancestor_state);
    if (task_properties::is_goal_state(task_proxy, state)) {
        return 0;
    } else {
		cout << "I am computing something!" << endl;
		cout << "  There are " << task_proxy.get_operators().size() << " operators." << endl;
		return 1;
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
