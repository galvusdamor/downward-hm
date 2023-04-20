#ifndef SYMBOLIC_SYM_VARIABLES_H
#define SYMBOLIC_SYM_VARIABLES_H


#include "../task_proxy.h"
#include "../tasks/root_task.h"
#include "../plugins/options.h"
#include <cuddObj.hh>


class GlobalState;


namespace options {
class Options;
class OptionParser;
} // namespace options


namespace symbolic {

struct BDDError {};
extern void exceptionError(std::string message);


class SymbolicHMBDDs {
    // Use task_proxy to access task information.
    TaskProxy task_proxy;

    const long cudd_init_nodes;          // Number of initial nodes
    const long cudd_init_cache_size;     // Initial cache size
    const long cudd_init_available_memory; // Maximum available memory (bytes)

    Cudd *manager; // manager associated with this symbolic search

    BDD current_state;   // BDDs associated with the precondition of a predicate

    int max_preconditions; // maximum number of preconditions in the task
    int max_effects; // maximum number of effects in the task

    int num_facts;
    int num_fact_bits;



public:
    SymbolicHMBDDs(const TaskProxy &task);

    void init();

    void current_state_to_dot(const std::string &filename);
    void to_dot(const std::string &filename, BDD bdd);
    BDD fact_to_bdd(int fact, int copy);

    int get_var_num(int bit, int copy);

};

}


#endif
