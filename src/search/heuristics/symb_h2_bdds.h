#ifndef SYMBOLIC_SYM_VARIABLES_H2
#define SYMBOLIC_SYM_VARIABLES_H2


#include "../task_proxy.h"
#include "../tasks/root_task.h"
#include "../plugins/options.h"
#include <cuddObj.hh>
#include <map>


class GlobalState;




namespace symbolic_2 {

struct BDDError {};
extern void exceptionError(std::string message);


class SymbolicH2BDDs {
    // Use task_proxy to access task information.
    TaskProxy task_proxy;
    int m = 2;

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
    int num_precondition_sets; // number of sets of preconditions
    int num_implicit_precondition_sets; // number of sets of implicit preconditions
    int num_facts; // number of facts in the task
    int num_fact_bits; // number of bits needed to represent the facts
    int num_set_bits; // number of set bits in the current state

    std::vector<std::vector<BDD>> sets_bdd_vars; // BDDs for each set
    std::vector<std::vector<std::vector<BDD>>> facts_bdd_vars ;// BDDs for each fact

    std::map<std::pair<int, int>, int> fact_map; // map from var and if to fact number

    std::vector<std::vector<int>> implicit_removes; // implicit removes for each fact


    

    // void generate_sets_util(std::vector<int>& nums, int m, std::vector<int>& current_set, std::vector<std::vector<int>>& result, int start_index);
    // std::vector<std::vector<int>> generate_sets(std::vector<int>& nums, int m);

    // void generate_sets_without_permutations_util(std::vector<int>& nums, int m, std::vector<int>& current_set, std::vector<std::vector<int>>& result, int start_index);
    // std::vector<std::vector<int>> generate_sets_without_permutations(std::vector<int>& nums, int m);


public:
    SymbolicH2BDDs(const TaskProxy &task);
    SymbolicH2BDDs(const TaskProxy &task, int m);

    void init();

    int calculate_heuristic(State state);

    void to_dot(const std::string &filename, BDD bdd);
    BDD fact_to_bdd(int fact, int fact_place, int copy);
    BDD set_to_bdd(std::vector<int> facts, int copy);
    BDD createBiimplicationBDD(std::vector<std::vector<BDD>> vars);

    void create_current_state_bdd(State state);
    void create_state_copy_bdd();
    void create_goal_bdd();
    void create_operators_bdd();


    int get_var_num(int bit, int fact_place, int copy);

};

}


#endif
