#include "symb_hm_bdds.h"

#include "../plugins/options.h"
#include <cmath>
#include <cstdio>
#include "../utils/logging.h"


using namespace std;
using plugins::Options;

namespace symbolic {

void exceptionError(string /*message*/) {
    // utils::g_log << message << endl;
    throw BDDError();
}

SymbolicHMBDDs::SymbolicHMBDDs(const TaskProxy &task)
    : task_proxy(task),
      cudd_init_nodes(16000000L), cudd_init_cache_size(16000000L),
      cudd_init_available_memory(0L) {}

void SymbolicHMBDDs::init() {
    // get the number of facts
    num_facts = 0;
    for (size_t i = 0; i < task_proxy.get_variables().size(); ++i) {
        num_facts += task_proxy.get_variables()[i].get_domain_size();
    }

    // get the number of fact bits
    num_fact_bits = ceil(log2(num_facts));

    // get the max amount of preconditions
    max_preconditions = 0;
    for (size_t i = 0; i < task_proxy.get_operators().size(); ++i) {
        int num_preconditions = task_proxy.get_operators()[i].get_preconditions().size();
        if (num_preconditions > max_preconditions) {
            max_preconditions = num_preconditions;
        }
    }

    // create cudd manager with variable size (max_preconditions+1) * ceil(log2(facts))
    manager = new Cudd(num_facts * (max_preconditions + 1) * num_fact_bits, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);

    // populate fact_bdd_vars
    for (int i = 0; i < max_preconditions + 1; ++i) {
        std::vector<BDD> fact_bdd_vars_copy;
        for (int j = 0; j < num_fact_bits; ++j) {
            fact_bdd_vars_copy.push_back(manager->bddVar(get_var_num(i, j)));
        }
        fact_bdd_vars.push_back(fact_bdd_vars_copy);
    }

    create_goal_bdd();

    // create pre_true_cube
    pre_true_cube = manager->bddVar(max_preconditions * num_fact_bits - 1);
    for (int i = 0; i < max_preconditions * num_fact_bits - 1; ++i) {
        pre_true_cube *= manager->bddVar(i);
    }

    return;
}


int SymbolicHMBDDs::calculate_heuristic(State state) {
    // create current state BDD
    create_current_state_bdd(state);
    BDD previous_state = current_state;

    int count = 0;
    while (true) {
        count++;
        // create state copy BDD
        create_state_copy_bdd();

        // create BDDs for each operator
        create_operators_bdd();

        // get effects and append to current state
        current_state += operators.AndAbstract(state_copy, pre_true_cube).SwapVariables(fact_bdd_vars[max_preconditions], fact_bdd_vars[0]);

        // check if goal is reachable
        if (goal <= current_state) {
            return count;
        }

        // check if current state is equal to previous state
        if (current_state == previous_state) {
            return -1;
        }

        // set previous state to current state
        previous_state = current_state;
    }


}

/**
 * Creates a dot file for the given BDD.
*/
void SymbolicHMBDDs::to_dot(const std::string &filename, BDD bdd) {
    ADD add = bdd.Add();
    FILE *fp = fopen(filename.c_str(), "w");
    DdNode **ddnodearray = (DdNode **)malloc(sizeof(add.getNode()));
    ddnodearray[0] = add.getNode();

    vector<string> var_names(num_fact_bits * (max_preconditions + 1));
    for (int i = 0; i < num_fact_bits; ++i) {
        for (int j = 0; j < max_preconditions; ++j) {
            std::string name = "pre" + to_string(j) + "_bit" + to_string(i);
            var_names[get_var_num(i, j)] = name;
        }
        std::string name = "eff_bit" + to_string(i);
        var_names[get_var_num(i, max_preconditions)] = name;
    }

    vector<char *> inames(num_fact_bits * (max_preconditions + 1));
    for (int i = 0; i < num_fact_bits * (max_preconditions + 1); ++i) {
        inames[i] = &var_names[i].front();
    }

    Cudd_DumpDot(manager->getManager(), 1, ddnodearray, inames.data(), NULL, fp);
    fclose(fp);
    free(ddnodearray);
}

/**
 * Converts a fact to a BDD (by encoding the bits of the fact)
*/
BDD SymbolicHMBDDs::fact_to_bdd(int fact, int copy) {
    BDD bdd = manager->bddOne();
    for (int i = 0; i < num_fact_bits; ++i) {
        if (fact & (1 << i)) {
            bdd *= manager->bddVar(get_var_num(num_fact_bits-i-1, copy));
        } else {
            bdd *= !manager->bddVar(get_var_num(num_fact_bits-i-1, copy));
        }
    }
    return bdd;
}

int SymbolicHMBDDs::get_var_num(int bit, int copy) {
    return bit + (copy * num_fact_bits);
}

/**
 * Creates the current_state BDD
 * size: num_fact_bits
 * use: It can be used to check if a fact is in the current state
 * made: for each fact, create a BDD with the fact bit set to 1 and the rest to 0
*/
void SymbolicHMBDDs::create_current_state_bdd(State state) {
    current_state = manager->bddZero();
    int fact = 0;
    for (int i = 0; i < task_proxy.get_variables().size(); ++i) {
        int varValue = state[i].get_value();
        for (int j = 0; j < task_proxy.get_variables()[i].get_domain_size(); ++j) {
            if (j == varValue) {
                current_state += fact_to_bdd(fact, 0);
            }
            fact++;
        }
    }

    // to_dot("current_state.dot", current_state);
}

/**
 * Creates the state_copy BDD
 * size: num_fact_bits * max_preconditions
 * use: It can be conjuncted with the operators to get the operators that are applicable in the current state
 * made: copy the current state to all copies
*/
void SymbolicHMBDDs::create_state_copy_bdd() {

    state_copy = manager->bddOne();
    for (int i = 0; i < max_preconditions; ++i) {
        state_copy *= current_state.SwapVariables(fact_bdd_vars[0], fact_bdd_vars[i]);
    }

    // to_dot("state_copy.dot", state_copy);
}

void SymbolicHMBDDs::create_goal_bdd() {
    goal = manager->bddZero();
    int fact = 0;
    for (FactProxy fact : task_proxy.get_goals()) {
        int var = fact.get_variable().get_id();
        int val = fact.get_value();
        int fact_num = 0;
        for (int i = 0; i < var; ++i) {
            fact_num += task_proxy.get_variables()[i].get_domain_size();
        }
        fact_num += val;
        goal += fact_to_bdd(fact_num, 0);
    }

    // to_dot("goal.dot", goal);
}


/**
 * Creates the operators BDD
 * size: num_fact_bits * (max_preconditions + 1)
 * use: It show all the operators in the problem, and the effects that they have
 * made: for each precondition within an operator, create a BDD using the fact bits.
 *    The effects are all encoded together in one part at the end.
 *    Repeat this for each operator 
*/
void SymbolicHMBDDs::create_operators_bdd() {
    operators = manager->bddZero();
    for (size_t i = 0; i < task_proxy.get_operators().size(); ++i) {
        OperatorProxy op = task_proxy.get_operators()[i];
        BDD preconditionBDD = manager->bddOne();
        for (size_t i = 0; i < op.get_preconditions().size(); ++i) {
            FactProxy fact = op.get_preconditions()[i];
            int var = fact.get_variable().get_id();
            int val = fact.get_value();
            int fact_num = 0;
            for (size_t k = 0; k < var; ++k) {
                fact_num += task_proxy.get_variables()[k].get_domain_size();
            }
            fact_num += val;
            preconditionBDD *= fact_to_bdd(fact_num, i);
        }
        BDD effectBDD = manager->bddZero();
        int fact_num = 0;
        int effect = 0;
        int amount_effects = op.get_effects().size();
        for (size_t i = 0; i < task_proxy.get_variables().size(); ++i) {
            for (size_t j = 0; j < task_proxy.get_variables()[i].get_domain_size(); ++j) {
                if (op.get_effects()[effect].get_fact().get_value() == j && op.get_effects()[effect].get_fact().get_variable().get_id() == i) {
                    effectBDD += fact_to_bdd(fact_num, max_preconditions);
                    effect++;
                }
                fact_num++;
            }
        }

        // for (size_t i = 0; i < op.get_effects().size(); ++i) {
        //     FactProxy fact = op.get_effects()[i].get_fact();
        //     int var = fact.get_variable().get_id();
        //     int val = fact.get_value();
        //     int fact_num = 0;
        //     for (size_t k = 0; k < var; ++k) {
        //         fact_num += task_proxy.get_variables()[k].get_domain_size();
        //     }
        //     fact_num += val;
        //     cout << "fact_num: " << fact_num << endl;
        //     effectBDD += fact_to_bdd(fact_num, max_preconditions);
        // }


        operators += preconditionBDD * effectBDD;
    }

    // to_dot("operators.dot", operators);

}

}
