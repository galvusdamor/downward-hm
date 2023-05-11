#include "symb_hm_bdds.h"

#include "../plugins/options.h"
#include <cmath>
#include <cstdio>
#include "../utils/logging.h"
#include <stdlib.h>

using namespace std;
using plugins::Options;

// Returns factorial of n
int fact(int n)
{
      if(n==0)
      return 1;
    int res = 1;
    for (int i = 2; i <= n; i++)
        res = res * i;
    return res;
}


int nCr(int n, int r)
{
    return fact(n) / (fact(r) * fact(n - r));
}



vector<vector<int>> generate_permutations(vector<int> nums) {
    vector<vector<int>> result;
    if (nums.size() == 0) {
        return result;
    }
    if (nums.size() == 1) {
        result.push_back(nums);
        return result;
    }

    vector<int> nums_copy = nums;
    nums_copy.erase(nums_copy.begin());
    vector<vector<int>> permutations = generate_permutations(nums_copy);
    for (int j = 0; j < permutations.size(); j++) {
        for (int k = 0; k < permutations[j].size(); k++) {
            vector<int> permutation_copy = permutations[j];
            permutation_copy.insert(permutation_copy.begin() + k, nums[0]);
            result.push_back(permutation_copy);
        }
        permutations[j].push_back(nums[0]);
        result.push_back(permutations[j]);
    }

    return result;
}


void generate_sets_util(vector<int>& nums, int m, vector<int>& current_set, vector<vector<int>>& result, int start_index) {
    if (current_set.size() == m) {
        result.push_back(current_set);
        return;
    }

    for (int i = start_index; i < nums.size(); i++) {
        current_set.push_back(nums[i]);
        generate_sets_util(nums, m, current_set, result, i + 1);
        current_set.pop_back();
    }
}

vector<vector<int>> generate_sets(vector<int>& nums, int m) {
    vector<vector<int>> result;
    if (m >= nums.size()) {
        // return single answer
        result.push_back(nums);
        if (m > nums.size()) {
            // add padding
            for (int i = 0; i < m - nums.size(); i++) {
                result[0].push_back(nums[0]);
            }
        }
    } else {
        vector<int> current_set;
        generate_sets_util(nums, m, current_set, result, 0);
    }


    vector<vector<int>> permutations_result;
    for (int i = 0; i < result.size(); i++) {
        vector<vector<int>> permutations = generate_permutations(result[i]);
        for (int j = 0; j < permutations.size(); j++) {
            permutations_result.push_back(permutations[j]);
        }
    }
    return permutations_result;
}

void generate_sets_without_permutations_util(vector<int>& nums, int m, vector<int>& current_set, vector<vector<int>>& result, int start_index) {
    if (current_set.size() == m) {
        result.push_back(current_set);
        return;
    }

    for (int i = start_index; i < nums.size(); i++) {
        current_set.push_back(nums[i]);
        generate_sets_without_permutations_util(nums, m, current_set, result, i + 1);
        current_set.pop_back();
    }

}

vector<vector<int>> generate_sets_without_permutations(vector<int>& nums, int m) {
    if (m >= nums.size()) {
        // return single answer
        vector<vector<int>> result;
        result.push_back(nums);
        if (m > nums.size()) {
            // add padding
            for (int i = 0; i < m - nums.size(); i++) {
                result[0].push_back(nums[0]);
            }
        }
        return result;
    }
    vector<vector<int>> result;
    vector<int> current_set;
    generate_sets_without_permutations_util(nums, m, current_set, result, 0);
    return result;
}




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
#include <stdlib.h>
    // This has not yet been implemented
    assert(max_preconditions < m);

    if (m == 1) {
        num_precondition_sets = max_preconditions;
        num_implicit_precondition_sets = 0;
    } else {
        // numer of sets is ncr of max_preconditions and m
        num_precondition_sets = nCr(max_preconditions, m);
        num_implicit_precondition_sets = nCr(max_preconditions, m - 1);
        for (int i = 0; i < m - 1; ++i) {
            num_implicit_precondition_sets += nCr(max_preconditions, i);
        }
    }


    // create cudd manager with variable size (max_preconditions+num_implicit_precondition_sets+1) * ceil(log2(facts))
    manager = new Cudd((num_precondition_sets + num_implicit_precondition_sets + 1) * num_set_bits, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);

    cout << "num_precondition_sets: " << num_precondition_sets << endl;

    // populate fact_bdd_vars and set_bdd_vars
    for (int i = 0; i < max_preconditions + 1; ++i) {
        std::vector<BDD> set_bdd_vars_copy;
        std::vector<std::vector<BDD>> fact_bdd_vars_copy1;
        for (int k = 0; k < m; ++k) {
            std::vector<BDD> fact_bdd_vars_copy;
            for (int j = 0; j < num_fact_bits; ++j) {
                fact_bdd_vars_copy.push_back(manager->bddVar(get_var_num(j, k, i)));
                set_bdd_vars_copy.push_back(manager->bddVar(get_var_num(j, k, i)));
            }
            // append to fact_bdd_vars
            fact_bdd_vars_copy1.push_back(fact_bdd_vars_copy);
        }
        facts_bdd_vars.push_back(fact_bdd_vars_copy1);
        sets_bdd_vars.push_back(set_bdd_vars_copy);
    }



    vector<int> nums = {1, 2, 3};

    cout << "Generating sets with permutations:" << endl;
    vector<vector<int>> sets = generate_sets(nums, m);

    for (const auto& set : sets) {
        for (const auto& num : set) {
            cout << num << " ";
        }
        cout << endl;
    }

    cout << "Generating sets without permutations:" << endl;
    vector<vector<int>> sets_without_permutations = generate_sets_without_permutations(nums, m);

    for (const auto& set : sets_without_permutations) {
        for (const auto& num : set) {
            cout << num << " ";
        }
        cout << endl;
    }
    exit(0);






    create_goal_bdd();

    // create pre_true_cube
    pre_true_cube = manager->bddVar(max_preconditions * num_set_bits - 1);
    for (int i = 0; i < max_preconditions * num_set_bits - 1; ++i) {
        pre_true_cube *= manager->bddVar(i);
    }
    
    // create BDDs for each operator
    create_operators_bdd();

    return;
}


int SymbolicHMBDDs::calculate_heuristic(State state) {
    // create current state BDD
    create_current_state_bdd(state);
    BDD previous_state = current_state;

    // TESTING

    // OperatorProxy op = task_proxy.get_operators()[0];
    // // get list of precondition numbers
    // std::vector<int> precondition_nums;
    // for (size_t i = 0; i < op.get_preconditions().size(); ++i) {
    //     precondition_nums.push_back(fact_map[make_pair(op.get_preconditions()[i].get_variable().get_id(), op.get_preconditions()[i].get_value())]);
    // }
    // // make list of precondition sets
    // std::vector<std::vector<int>> precondition_sets = get_all_ncr_sets(precondition_nums);

    // // print precondition sets
    // for (size_t i = 0; i < precondition_sets.size(); ++i) {
    //     for (size_t j = 0; j < precondition_sets[i].size(); ++j) {
    //         std::cout << precondition_sets[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // // create bdd for the precondition
    // BDD precondition_bdd = manager->bddOne();
    // for (size_t i = 0; i < precondition_sets.size(); ++i) {
    //     BDD precondition_set_bdd = manager->bddOne();
    //     for (size_t j = 0; j < precondition_sets[i].size(); ++j) {
    //         precondition_set_bdd *= fact_to_bdd(precondition_sets[i][j], j, i);
    //     }
    //     precondition_bdd *= precondition_set_bdd;
    // }

    // // get effects
    // std::vector<int> effect_nums;
    // for (size_t i = 0; i < op.get_effects().size(); ++i) {
    //     effect_nums.push_back(fact_map[make_pair(op.get_effects()[i].get_fact().get_variable().get_id(), op.get_effects()[i].get_fact().get_value())]);
    // }
    // // make list of effect sets
    // std::vector<std::vector<int>> effect_sets = get_all_npr_sets(effect_nums);






    int count = 0;
    while (true) {
        count++;
        // create state copy BDD
        create_state_copy_bdd();


        // get effects and append to current state
        current_state += operators.AndAbstract(state_copy, pre_true_cube).SwapVariables(sets_bdd_vars[max_preconditions], sets_bdd_vars[0]);
        
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

/**
 * Converts a set of facts to a BDD (by encoding the bits of the facts)
*/
BDD SymbolicHMBDDs::set_to_bdd(std::vector<int> facts, int copy) {
    BDD bdd = manager->bddOne();
    for (int i = 0; i < facts.size(); ++i) {
        bdd *= fact_to_bdd(facts[i], i, copy);
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
    current_state = manager->bddOne();
    std::vector<int> facts;
    for (FactProxy fact : state) {
        facts.push_back(fact_map[make_pair(fact.get_variable().get_id(), fact.get_value())]);
    }

    // create a bdd with each fact
    BDD fact_bdd = manager->bddZero();
    for (int fact : facts) {
        fact_bdd += fact_to_bdd(fact, 0, 0);
    }

    // copy over to each fact in the set
    for (int i = 0; i < m; ++i) {
    }

    to_dot("current_state.dot", current_state);
    exit(0);
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
        state_copy *= current_state.SwapVariables(sets_bdd_vars[0], sets_bdd_vars[i]);
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
    std::vector<std::vector<int>> all_sets = generate_sets_without_permutations(facts, m);

    // create BDD for each set
    for (std::vector<int> set : all_sets) {
        goal += set_to_bdd(set, 0);
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
