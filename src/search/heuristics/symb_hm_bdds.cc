#include "symb_hm_bdds.h"

#include "../plugins/options.h"
#include <cmath>
#include <cstdio>


using namespace std;
using plugins::Options;

namespace symbolic {

void exceptionError(string /*message*/) {
    // utils::g_log << message << endl;
    throw BDDError();
}

SymbolicHMBDDs::SymbolicHMBDDs(const Options &opts,
                           const TaskProxy &task)
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
    int max_preconditions = 0;
    for (size_t i = 0; i < task_proxy.get_operators().size(); ++i) {
        int num_preconditions = task_proxy.get_operators()[i].get_preconditions().size();
        if (num_preconditions > max_preconditions) {
            max_preconditions = num_preconditions;
        }
    }

    // create cudd manager with variable size (max_preconditions+1) * ceil(log2(facts))
    manager = new Cudd(num_facts * (max_preconditions + 1) * ceil(log2(num_facts)), 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);


    // // create a BDD to represent the current state
    // BDD state = manager->bddOne();
	// int varNum = 0;
	// for (size_t i = 0; i < task_proxy.get_variables().size(); ++i) {
    //     int varValue = task_proxy.get_initial_state()[i].get_value(); //   get_variables()[i].get_id();
    //     int bddVariableNumber = varNum + varValue;
	// 	state = state * manager->bddVar(bddVariableNumber);
	// 	for (size_t j = 0; j < task_proxy.get_variables()[i].get_domain_size(); ++j) {
	// 		if (j != varValue)
	// 			state = state * ! manager->bddVar(varNum + j);

	// 	}

	// 	varNum += task_proxy.get_variables()[i].get_domain_size();

    // }

    // create current state BDD
    current_state = manager->bddZero();
    int fact = 0;
    State state = task_proxy.get_initial_state();
    for (size_t i = 0; i < task_proxy.get_variables().size(); ++i) {
        int varValue = state[i].get_value();
        for (size_t j = 0; j < task_proxy.get_variables()[i].get_domain_size(); ++j) {
            if (j == varValue) {
                current_state += fact_to_bdd(fact);
                to_dot("asdfadsf.dot", fact_to_bdd(fact));
                cout << "fact: " << fact << endl;
            }
            fact++;
        }
    }
    current_state_to_dot("current_state.dot");



    return;
}

void SymbolicHMBDDs::current_state_to_dot(const std::string &filename) {
    ADD add = current_state.Add();
    FILE *fp = fopen(filename.c_str(), "w");
    DdNode **ddnodearray = (DdNode **)malloc(sizeof(add.getNode()));
	ddnodearray[0] = add.getNode();
	Cudd_DumpDot(manager->getManager(), 1, ddnodearray, NULL, NULL, fp);
    free(ddnodearray);
    fclose(fp);
}

void SymbolicHMBDDs::to_dot(const std::string &filename, BDD bdd) {
    ADD add = bdd.Add();
    FILE *fp = fopen(filename.c_str(), "w");
    DdNode **ddnodearray = (DdNode **)malloc(sizeof(add.getNode()));
    ddnodearray[0] = add.getNode();
    Cudd_DumpDot(manager->getManager(), 1, ddnodearray, NULL, NULL, fp);
    free(ddnodearray);
    fclose(fp);
}

BDD SymbolicHMBDDs::fact_to_bdd(int fact) {
    BDD bdd = manager->bddOne();
    for (int i = 0; i < num_fact_bits; ++i) {
        if (fact & (1 << i)) {
            bdd *= manager->bddVar(i);
        } else {
            bdd *= !manager->bddVar(i);
        }
    }
    return bdd;
}

}
