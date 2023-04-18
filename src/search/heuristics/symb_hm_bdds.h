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



public:
    SymbolicHMBDDs(const plugins::Options &opts,
                 const TaskProxy &task);

    void init();

};

}


#endif
