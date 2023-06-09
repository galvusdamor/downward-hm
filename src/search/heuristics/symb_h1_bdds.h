#ifndef SYMBOLIC_SYM_VARIABLES_H1
#define SYMBOLIC_SYM_VARIABLES_H1


#include "../task_proxy.h"
#include "../tasks/root_task.h"
#include <cuddObj.hh>
#include <map>


class GlobalState;



namespace symbolic_1 {

struct BDDError {};
extern void exceptionError(std::string message);


class SymbolicH1BDDs {
    // Use task_proxy to access task information.
    TaskProxy task_proxy;
    int var_order;

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
    int num_facts; // number of facts in the task
    int num_fact_bits; // number of bits needed to represent the facts

    std::vector<std::vector<BDD>> fact_bdd_vars; // BDDs for each fact

    std::map<std::pair<int, int>, int> fact_map; // map from var and if to fact number

public:
    SymbolicH1BDDs(const TaskProxy &task, int var_order);

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
