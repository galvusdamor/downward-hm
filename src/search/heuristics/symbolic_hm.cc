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
	assert(m == 1); // We can for now only handle m=1

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


    bdds = make_shared<symbolic::SymbolicHMBDDs>(task_proxy);

    bdds->init();

    //vars to dot file
    

    cout << "  There are " << task_proxy.get_operators().size() << " operators." << endl;
    cout << "  There are " << task_proxy.get_variables().size() << " variables." << endl;
    cout << "  There are " << task_proxy.get_axioms().size() << " axioms." << endl;
    cout << "  There are " << task_proxy.get_goals().size() << " goals." << endl;
    // print the amount of facts
    int num_facts = 0;
    for (VariableProxy var : task_proxy.get_variables()) {
        num_facts += var.get_domain_size();
    }
    cout << "  There are " << num_facts << " facts." << endl;

    State state = task_proxy.get_initial_state();
    // print the amount of facts that are true
    int num_true_facts = 0;
    for (FactProxy fact : state) {
        num_true_facts++;
    }
    cout << "  There are " << num_true_facts << " true facts." << endl;

    // print the symbols of the true facts
    for (FactProxy fact : state) {
        cout << "  " << fact.get_variable().get_name() << "=" << fact.get_value() << endl;
    }

    // print the amount of operators that are applicable
    int num_applicable_ops = 0;
    for (OperatorProxy op : task_proxy.get_operators()) {
        // get preconditions
        vector<FactPair> preconditions;
        for (FactProxy pre : op.get_preconditions()) {
            preconditions.push_back(pre.get_pair());
        }
        // check if preconditions are true
        bool applicable = true;
        for (FactPair pre : preconditions) {
            if ((!state[pre.var].get_value()) == pre.value) {
                applicable = false;
                break;
            }
        }
        if (applicable) {
            num_applicable_ops++;
        }
        cout << "  " << op.get_name() << " is applicable: " << applicable << endl;
        cout << "    amount effects: " << op.get_effects().size() << endl;



    }
    cout << "  There are " << num_applicable_ops << " applicable operators." << endl;




	if (log.is_at_least_normal()) {
        log << "Using symbolic h^" << m << "." << endl;
    }
}



int SymbolicHMHeuristic::compute_heuristic(const State &ancestor_state) {
    State state = convert_ancestor_state(ancestor_state);
    if (task_properties::is_goal_state(task_proxy, state)) {
        return 0;
    } else {
		// cout << "I am computing something!" << endl;
		// cout << "  There are " << task_proxy.get_operators().size() << " operators." << endl;
        // cout << "  There are " << task_proxy.get_variables().size() << " variables." << endl;
        // cout << "  There are " << task_proxy.get_axioms().size() << " axioms." << endl;
        // cout << "  There are " << task_proxy.get_goals().size() << " goals." << endl;
        // // print the amount of facts
        // int num_facts = 0;
        // for (VariableProxy var : task_proxy.get_variables()) {
        //     num_facts += var.get_domain_size();
        // }
        // cout << "  There are " << num_facts << " facts." << endl;
        // // print the amount of facts that are true
        // int num_true_facts = 0;
        // for (FactProxy fact : state) {
        //     num_true_facts++;
        // }
        // cout << "  There are " << num_true_facts << " true facts." << endl;
        // // print the amount of operators that are applicable
        // int num_applicable_ops = 0;
        // for (OperatorProxy op : task_proxy.get_operators()) {
        //     // get preconditions
        //     vector<FactPair> preconditions;
        //     for (FactProxy pre : op.get_preconditions()) {
        //         preconditions.push_back(pre.get_pair());
        //     }
        //     // check if preconditions are true
        //     bool applicable = true;
        //     for (FactPair pre : preconditions) {
        //         if (!state[pre.var].get_value() == pre.value) {
        //             applicable = false;
        //             break;
        //         }
        //     }
        //     if (applicable) {
        //         num_applicable_ops++;
        //     }
        // }
        // cout << "  There are " << num_applicable_ops << " applicable operators." << endl;
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
