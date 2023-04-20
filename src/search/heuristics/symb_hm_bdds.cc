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

    // get the max amount of effects
    max_effects = 0;
    for (size_t i = 0; i < task_proxy.get_operators().size(); ++i) {
        int num_effects = task_proxy.get_operators()[i].get_effects().size();
        if (num_effects > max_effects) {
            max_effects = num_effects;
        }
    }
    cout << "max preconditions: " << max_preconditions << endl;
    cout << "max effects: " << max_effects << endl;

    // create cudd manager with variable size (max_preconditions+1) * ceil(log2(facts))
    manager = new Cudd(num_facts * (max_preconditions + max_effects) * num_fact_bits, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);

    // create current state BDD
    current_state = manager->bddZero();
    int fact = 0;
    State state = task_proxy.get_initial_state();
    for (int i = 0; i < task_proxy.get_variables().size(); ++i) {
        int varValue = state[i].get_value();
        for (int j = 0; j < task_proxy.get_variables()[i].get_domain_size(); ++j) {
            if (j == varValue) {
                current_state += fact_to_bdd(fact, 0);
            }
            fact++;
        }
    }
    current_state_to_dot("current_state.dot");

    // create BDDs for each operator
    BDD preconditions = manager->bddZero();
    for (OperatorProxy op : task_proxy.get_operators()) {
        BDD precondition = manager->bddOne();
        for (int i = 0; i < op.get_preconditions().size(); ++i) {
            FactProxy fact = op.get_preconditions()[i];
            int var = fact.get_variable().get_id();
            int val = fact.get_value();
            int fact_num = 0;
            for (size_t k = 0; k < var; ++k) {
                fact_num += task_proxy.get_variables()[k].get_domain_size();
            }
            fact_num += val;
            precondition *= fact_to_bdd(fact_num, i);
            if (i == 0) {
                cout << "fact: " << fact_num << endl;
            }
        }
        // combine all preconditions
        preconditions += precondition;
    }
    to_dot("preconditions.dot", preconditions);


    return;
}

void SymbolicHMBDDs::current_state_to_dot(const std::string &filename) {
    to_dot(filename, current_state);
}

void SymbolicHMBDDs::to_dot(const std::string &filename, BDD bdd) {
    ADD add = bdd.Add();
    FILE *fp = fopen(filename.c_str(), "w");
    DdNode **ddnodearray = (DdNode **)malloc(sizeof(add.getNode()));
    ddnodearray[0] = add.getNode();

    vector<string> var_names(num_fact_bits * (max_preconditions + max_effects));
    for (int i = 0; i < num_fact_bits; ++i) {
        for (int j = 0; j < max_preconditions; ++j) {
            std::string name = "pre" + to_string(j) + "_bit" + to_string(i);
            var_names[get_var_num(i, j)] = name;
        }
        for (int j = max_preconditions; j < (max_preconditions + max_effects); ++j) {
            std::string name = "pre" + to_string(j) + "_bit" + to_string(i);
            var_names[get_var_num(i, j)] = name;
        }
    }

    vector<char *> inames(num_fact_bits * (max_preconditions + max_effects));
    for (int i = 0; i < num_fact_bits * (max_preconditions + max_effects); ++i) {
        inames[i] = &var_names[i].front();
    }

    Cudd_DumpDot(manager->getManager(), 1, ddnodearray, inames.data(), NULL, fp);
    free(ddnodearray);
    fclose(fp);
}

BDD SymbolicHMBDDs::fact_to_bdd(int fact, int copy) {
    BDD bdd = manager->bddOne();
    for (int i = 0; i < num_fact_bits; ++i) {
        if (fact & (1 << i)) {
            bdd *= manager->bddVar(get_var_num(i, copy));
        } else {
            bdd *= !manager->bddVar(get_var_num(i, copy));
        }
    }
    return bdd;
}

int SymbolicHMBDDs::get_var_num(int bit, int copy) {
    return bit + (copy * num_fact_bits);
}


}
