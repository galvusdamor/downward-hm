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
      cudd_init_available_memory(0L) {
        m=1;
}

SymbolicHMBDDs::SymbolicHMBDDs(const TaskProxy &task, int mm)
    : task_proxy(task),
      cudd_init_nodes(16000000L), cudd_init_cache_size(16000000L),
      cudd_init_available_memory(0L) {
        m = mm;
}

void SymbolicHMBDDs::init() {
    // get the number of facts and create fact_map
    num_facts = 0;
    implicit_removes = std::vector<std::vector<int>>(task_proxy.get_variables().size(), std::vector<int>());
    fact_map = std::map<std::pair<int, int>, int>();
    int fact_num = 0;
    for (size_t i = 0; i < task_proxy.get_variables().size(); ++i) {
        num_facts += task_proxy.get_variables()[i].get_domain_size();
        for (int j = 0; j < task_proxy.get_variables()[i].get_domain_size(); ++j) {
            fact_map[make_pair(i, j)] = fact_num;
            implicit_removes[i].push_back(fact_num);
            fact_num++;
        }
    }

    // get the number of fact bits
    num_fact_bits = ceil(log2(num_facts));
    num_set_bits = num_fact_bits * m;

    // get the max amount of preconditions
    max_preconditions = 0;
    for (size_t i = 0; i < task_proxy.get_operators().size(); ++i) {
        int num_preconditions = task_proxy.get_operators()[i].get_preconditions().size();
        if (num_preconditions > max_preconditions) {
            max_preconditions = num_preconditions;
        }
    }


    // create cudd manager with variable size (max_preconditions+1) * ceil(log2(facts))
    manager = new Cudd((max_preconditions + 1) * num_set_bits, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);

    // populate fact_bdd_vars
    for (int i = 0; i < max_preconditions + 1; ++i) {
        std::vector<BDD> fact_bdd_vars_copy;
        for (int k = 0; k < m; ++k) {
            for (int j = 0; j < num_fact_bits; ++j) {
                fact_bdd_vars_copy.push_back(manager->bddVar(get_var_num(j, k, i)));
            }
        }
        fact_bdd_vars.push_back(fact_bdd_vars_copy);
    }

    create_goal_bdd();

    // create pre_true_cube
    pre_true_cube = manager->bddVar(max_preconditions * num_set_bits - 1);
    for (int i = 0; i < max_preconditions * num_set_bits - 1; ++i) {
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

    vector<string> var_names(num_set_bits * (max_preconditions + 1));
    for (int i = 0; i < num_fact_bits; ++i) {
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < max_preconditions; ++k) {
                std::string name = "pre" + to_string(k) + "_fact" + to_string(j) + "_bit" + to_string(i);
                var_names[get_var_num(i, j, k)] = name;
            }
            std::string name = "eff_fact" + to_string(j) + "_bit" + to_string(i);
            var_names[get_var_num(i, j, max_preconditions)] = name;
        }
    }

    vector<char *> inames(num_set_bits * (max_preconditions + 1));
    for (int i = 0; i < num_set_bits * (max_preconditions + 1); ++i) {
        inames[i] = &var_names[i].front();
    }

    Cudd_DumpDot(manager->getManager(), 1, ddnodearray, inames.data(), NULL, fp);
    fclose(fp);
    free(ddnodearray);
}

/**
 * Converts a fact to a BDD (by encoding the bits of the fact)
*/
BDD SymbolicHMBDDs::fact_to_bdd(int fact, int fact_place, int copy) {
    BDD bdd = manager->bddOne();
    for (int i = 0; i < num_fact_bits; ++i) {
        if (fact & (1 << i)) {
            bdd *= manager->bddVar(get_var_num(num_fact_bits-i-1, fact_place, copy));
        } else {
            bdd *= !manager->bddVar(get_var_num(num_fact_bits-i-1, fact_place, copy));
        }
    }
    return bdd;
}

int SymbolicHMBDDs::get_var_num(int bit, int fact_place, int copy) {
    return bit + (fact_place * num_fact_bits) + (copy * num_set_bits);
}

/**
 * Creates the current_state BDD
 * size: num_fact_bits
 * use: It can be used to check if a fact is in the current state
 * made: for each fact, create a BDD with the fact bit set to 1 and the rest to 0
*/
void SymbolicHMBDDs::create_current_state_bdd(State state) {
    current_state = manager->bddZero();
    std::vector<int> facts;
    for (FactProxy fact : state) {
        facts.push_back(fact_map[make_pair(fact.get_variable().get_id(), fact.get_value())]);
    }

    // create all possible sets of m facts
    std::vector<std::vector<int>> all_sets = get_all_sets(facts);

    // create BDD for each set
    for (std::vector<int> set : all_sets) {
        BDD set_bdd = manager->bddOne();
        for (int i = 0; i < m; ++i) {
            set_bdd *= fact_to_bdd(set[i], i, 0);
        }
        current_state += set_bdd;
    }

    to_dot("current_state.dot", current_state);
}

std::vector<std::vector<int>> SymbolicHMBDDs::get_all_sets(std::vector<int> set) {
    std::vector<std::vector<int>> all_sets;
    std::vector<int> current_set;
    get_all_sets_rec(set, 0, current_set, all_sets);
    return all_sets;
}

void SymbolicHMBDDs::get_all_sets_rec(std::vector<int> set, int index, std::vector<int> current_set, std::vector<std::vector<int>> &all_sets) {
    if (current_set.size() == m) {
        all_sets.push_back(current_set);
        return;
    }
    if (index == set.size()) {
        return;
    }
    current_set.push_back(set[index]);
    get_all_sets_rec(set, index, current_set, all_sets);
    current_set.pop_back();
    get_all_sets_rec(set, index+1, current_set, all_sets);
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
    std::vector<int> facts;

    for (FactProxy fact : task_proxy.get_goals()) {
        facts.push_back(fact_map[make_pair(fact.get_variable().get_id(), fact.get_value())]);
    }

    // create all possible sets of m facts
    std::vector<std::vector<int>> all_sets = get_all_sets(facts);

    // create BDD for each set
    for (std::vector<int> set : all_sets) {
        BDD set_bdd = manager->bddOne();
        for (int i = 0; i < m; ++i) {
            set_bdd *= fact_to_bdd(set[i], i, 0);
        }
        goal += set_bdd;
    }

    to_dot("goal.dot", goal);
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
            preconditionBDD *= fact_to_bdd(fact_map[make_pair(fact.get_variable().get_id(), fact.get_value())], 0, i);
        }
        BDD effectBDD = manager->bddZero();
        for (size_t i = 0; i < op.get_effects().size(); ++i) {
            FactProxy fact = op.get_effects()[i].get_fact();
            effectBDD += fact_to_bdd(fact_map[make_pair(fact.get_variable().get_id(), fact.get_value())], 0, max_preconditions);
        }


        operators += preconditionBDD * effectBDD;
    }

    // to_dot("operators.dot", operators);

}

}
