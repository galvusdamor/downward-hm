#include "symb_hm_bdds.h"




using namespace std;
using plugins::Options;

namespace symbolic {


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

    // create cudd manager
    manager = new Cudd(cudd_init_nodes, cudd_init_cache_size, cudd_init_available_memory);
    
    // create a BDD to represent the current state
    current_state = manager->bddOne();

    // set the BDD values to the current true facts
    for (size_t i = 0; i < task_proxy.get_variables().size(); ++i) {
        int var = task_proxy.get_variables()[i].get_id();
        int val = task_proxy.get_variables()[i].get_fact(task_proxy.get_variables()[i].get_domain_size() - 1).get_value();
        int index = var * num_facts + val;
        current_state = current_state * manager->bddVar(index);
    }

    return;
    

}
}