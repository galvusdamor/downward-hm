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

    BDD current_state;   // BDD for the current state
    BDD state_copy;      // BDD for the state copy (current state copied into each copy)
    BDD goal;            // BDD for the goal
    BDD operators;   // BDD for the preconditions
    BDD pre_true_cube;   // BDD for a true cube the size of the preconditions

    int max_preconditions; // maximum number of preconditions in the task

    int num_facts;
    int num_fact_bits;

    std::vector<std::vector<BDD>> fact_bdd_vars; // BDDs for each fact



public:
    SymbolicHMBDDs(const TaskProxy &task);

    void init();

    int calculate_heuristic(State state);

    void to_dot(const std::string &filename, BDD bdd);
    BDD fact_to_bdd(int fact, int copy);


    void create_current_state_bdd(State state);
    void create_state_copy_bdd();
    void create_goal_bdd();
    void create_operators_bdd();


    int get_var_num(int bit, int copy);

};

}


#endif
