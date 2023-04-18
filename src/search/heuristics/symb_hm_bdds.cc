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
    // get the number of variables
    int num_fd_vars = task_proxy.get_variables().size();
    // get the number of facts
    int num_facts = 0;
    for (size_t i = 0; i < task_proxy.get_variables().size(); ++i) {
        num_facts += task_proxy.get_variables()[i].get_domain_size();
    }
    // get the max amount of preconditions
    int max_preconditions = 0;
    for (size_t i = 0; i < task_proxy.get_operators().size(); ++i) {
        int num_preconditions = task_proxy.get_operators()[i].get_preconditions().size();
        if (num_preconditions > max_preconditions) {
            max_preconditions = num_preconditions;
        }
    }

    // create cudd manager with variable size (max_preconditions+1) * ceil(log2(facts))
    manager = Cudd_Init(num_fd_vars * (max_preconditions + 1) * ceil(log2(num_facts)), 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);


    // create a BDD to represent the current state
    DdNode *var = NULL;
    DdNode *first_var = NULL;
    for (size_t i = 0; i < task_proxy.get_variables().size(); ++i) {
        int varValue = task_proxy.get_variables()[i].get_id();
        for (size_t j = 0; j < task_proxy.get_variables()[i].get_domain_size(); ++j) {
            int factValue = task_proxy.get_variables()[i].get_fact(j).get_value();
            cout << "  " << task_proxy.get_variables()[i].get_name() << "=" << factValue << endl;
            if (j == 0 && i == 0) {
                first_var = Cudd_bddNewVar(manager);
                continue;
            }
            if (var == NULL) {
                var = Cudd_bddNewVar(manager);
                if (factValue == varValue) {
                    current_state = Cudd_bddAnd(manager, first_var, var);
                } else {
                    current_state = Cudd_bddAnd(manager, first_var, Cudd_Not(var));
                }
            } else {
                var = Cudd_bddNewVar(manager);
                if (factValue == varValue) {
                    current_state = Cudd_bddAnd(manager, current_state, var);
                } else {
                    current_state = Cudd_bddAnd(manager, current_state, Cudd_Not(var));
                }
            }
        }
    }

    // to dot file
    FILE *fp = fopen("current_state.dot", "w");
    Cudd_DumpDot(manager, 1, &current_state, NULL, NULL, fp);
    fclose(fp);




    // // set the BDD values to the current true facts
    // for (size_t i = 0; i < task_proxy.get_variables().size(); ++i) {
    //     int var = task_proxy.get_variables()[i].get_id();
    //     int val = task_proxy.get_variables()[i].get_fact(task_proxy.get_variables()[i].get_domain_size() - 1).get_value();
    //     int index = var * num_facts + val;
    //     current_state = current_state * manager->bddVar(index);
    // }

    return;
    

}
}