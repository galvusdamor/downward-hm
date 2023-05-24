#include "symb_h2_bdds.h"

#include <cmath>
#include <cstdio>
#include "../utils/logging.h"
#include <stdlib.h>

using namespace std;

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
        generate_sets_util(nums, m, current_set, result, i+1);
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


namespace symbolic_2 {

void exceptionError(string /*message*/) {
    // utils::g_log << message << endl;
    throw BDDError();
}

SymbolicH2BDDs::SymbolicH2BDDs(const TaskProxy &task)
    : task_proxy(task),
      cudd_init_nodes(16000000L), cudd_init_cache_size(16000000L),
      cudd_init_available_memory(0L) {
}

void SymbolicH2BDDs::init() {
    // get the number of facts and create fact_map
    num_facts = 0;
    fact_map = std::map<std::pair<int, int>, int>();
    int fact_num = 0;
    for (size_t i = 0; i < task_proxy.get_variables().size(); ++i) {
        num_facts += task_proxy.get_variables()[i].get_domain_size();
        int first_fact = fact_num;
        for (int j = 0; j < task_proxy.get_variables()[i].get_domain_size(); ++j) {
            fact_map[make_pair(i, j)] = fact_num;
            fact_num++;
        }
    }

    implicit_removes = vector<vector<int>>(num_facts, vector<int>());
    for (size_t i = 0; i < task_proxy.get_variables().size(); ++i) {
        for (int j = 0; j < task_proxy.get_variables()[i].get_domain_size(); ++j) {
            int fact = fact_map[make_pair(i, j)];
            for (int k = 0; k < task_proxy.get_variables()[i].get_domain_size(); ++k) {
                if (k != j) {
                    implicit_removes[fact].push_back(fact_map[make_pair(i, k)]);
                }
            }
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

    // This has not yet been implemented
    assert(max_preconditions < m);


    // numer of sets is ncr of max_preconditions and m
    num_precondition_sets = nCr(max_preconditions, m);
    num_implicit_precondition_sets = max_preconditions;

    // create cudd manager with variable size (max_preconditions+num_implicit_precondition_sets+1) * ceil(log2(facts))
    manager = new Cudd((num_precondition_sets + num_implicit_precondition_sets + 1) * num_set_bits, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);

    // populate fact_bdd_vars and set_bdd_vars
    for (int i = 0; i < num_precondition_sets + num_implicit_precondition_sets + 1; ++i) {
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

    create_goal_bdd();

    // create pre_true_cube
    pre_true_cube = manager->bddVar((num_precondition_sets + num_implicit_precondition_sets) * num_set_bits - 1);
    for (int i = 0; i < (num_precondition_sets + num_implicit_precondition_sets) * num_set_bits - 1; ++i) {
        pre_true_cube *= manager->bddVar(i);
    }
    
    // create BDDs for each operator
    create_operators_bdd();

    return;
}


int SymbolicH2BDDs::calculate_heuristic(State state) {
    // create current state BDD
    create_current_state_bdd(state);
    BDD previous_state = current_state;

    int count = 0;
    while (true) {
        // check if goal is reachable
        if (goal <= current_state) {
            return count;
        }
        
        count++;
        // create state copy BDD
        create_state_copy_bdd();

        // get effects and append to current state
        current_state += operators.AndAbstract(state_copy, pre_true_cube).SwapVariables(sets_bdd_vars[num_precondition_sets + num_implicit_precondition_sets], sets_bdd_vars[0]);
        
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
void SymbolicH2BDDs::to_dot(const std::string &filename, BDD bdd) {
    ADD add = bdd.Add();
    FILE *fp = fopen(filename.c_str(), "w");
    DdNode **ddnodearray = (DdNode **)malloc(sizeof(add.getNode()));
    ddnodearray[0] = add.getNode();

    vector<string> var_names(num_set_bits * (num_precondition_sets + num_implicit_precondition_sets + 1));
    for (int i = 0; i < num_fact_bits; ++i) {
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < num_precondition_sets; ++k) {
                std::string name = "pre" + to_string(k) + "_fact" + to_string(j) + "_bit" + to_string(i);
                var_names[get_var_num(i, j, k)] = name;
            }
            for (int k = num_precondition_sets; k < num_precondition_sets + num_implicit_precondition_sets; ++k) {
                std::string name = "impl_pre" + to_string(k - num_precondition_sets) + "_fact" + to_string(j) + "_bit" + to_string(i);
                var_names[get_var_num(i, j, k)] = name;
            }
            std::string name = "eff_fact" + to_string(j) + "_bit" + to_string(i);
            var_names[get_var_num(i, j, num_precondition_sets + num_implicit_precondition_sets)] = name;
        }
    }

    vector<char *> inames(num_set_bits * (num_precondition_sets + num_implicit_precondition_sets + 1));
    for (int i = 0; i < num_set_bits * (num_precondition_sets + num_implicit_precondition_sets + 1); ++i) {
        inames[i] = &var_names[i].front();
    }

    Cudd_DumpDot(manager->getManager(), 1, ddnodearray, inames.data(), NULL, fp);
    fclose(fp);
    free(ddnodearray);
}

/**
 * Converts a fact to a BDD (by encoding the bits of the fact)
*/
BDD SymbolicH2BDDs::fact_to_bdd(int fact, int fact_place, int copy) {
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
BDD SymbolicH2BDDs::set_to_bdd(std::vector<int> facts, int copy) {
    BDD bdd = manager->bddOne();
    for (int i = 0; i < facts.size(); ++i) {
        bdd *= fact_to_bdd(facts[i], i, copy);
    }
    return bdd;
}

int SymbolicH2BDDs::get_var_num(int bit, int fact_place, int copy) {
    return bit + (fact_place * num_fact_bits) + (copy * num_set_bits);
}

/**
 * Creates the current_state BDD
 * size: num_fact_bits
 * use: It can be used to check if a fact is in the current state
 * made: for each fact, create a BDD with the fact bit set to 1 and the rest to 0
*/
void SymbolicH2BDDs::create_current_state_bdd(State state) {
    std::vector<int> facts;
    for (FactProxy fact : state) {
        facts.push_back(fact_map[make_pair(fact.get_variable().get_id(), fact.get_value())]);
    }

    // create a bdd with each fact
    BDD fact_bdd = manager->bddZero();
    for (int fact : facts) {
        fact_bdd += fact_to_bdd(fact, 0, 0);
    }

    current_state = manager->bddOne();
    // copy over to each fact in the set
    for (int i = 0; i < m; ++i) {
        current_state *= fact_bdd.SwapVariables(facts_bdd_vars[0][0], facts_bdd_vars[0][i]);
    }

    // to_dot("current_state.dot", current_state);
}


BDD SymbolicH2BDDs::createBiimplicationBDD(std::vector<std::vector<BDD>> vars) {
    BDD bdd = manager->bddOne();
    for (int i = 0; i < vars[0].size(); ++i) {
        // biimp all with eachother
        BDD xnor = manager->bddOne();
        for (int j = 0; j < vars.size() - 1; ++j) {
            xnor *= vars[j][i].Xnor(vars[j+1][i]);
        }
        bdd *= xnor;

    }
    return bdd;
}

/**
 * Creates the state_copy BDD
 * size: num_fact_bits * (num_precondition_sets + num_implicit_precondition_sets)
 * use: It can be conjuncted with the operators to get the operators that are applicable in the current state
 * made: copy the current state to all copies
*/
void SymbolicH2BDDs::create_state_copy_bdd() {

    state_copy = manager->bddOne();
    for (int i = 0; i < num_precondition_sets + num_implicit_precondition_sets; ++i) {
        state_copy *= current_state.SwapVariables(sets_bdd_vars[0], sets_bdd_vars[i]);
    }

    // to_dot("state_copy.dot", state_copy);
}

void SymbolicH2BDDs::create_goal_bdd() {
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

    // to_dot("goal.dot", goal);
}


/**
 * Creates the operators BDD
 * size: num_fact_bits * (max_preconditions + 1)
 * use: It show all the operators in the problem, and the effects that they have
 * made: for each precondition within an operator, create a BDD using the fact bits.
 *    The effects are all encoded together in one part at the end.
 *    Repeat this for each operator 
*/
void SymbolicH2BDDs::create_operators_bdd() {
    operators = manager->bddZero();
    for (OperatorProxy op : task_proxy.get_operators()) {
        // get preconditions
        vector<int> preconditions;
        for (size_t i = 0; i < op.get_preconditions().size(); ++i) {
            FactProxy fact = op.get_preconditions()[i];
            preconditions.push_back(fact_map[make_pair(fact.get_variable().get_id(), fact.get_value())]);
        }
        // get precondition sets
        vector<vector<int>> precondition_sets = generate_sets_without_permutations(preconditions, m);
        // create precondition BDD
        BDD preconditionBDD = manager->bddOne();
        for (size_t i = 0; i < precondition_sets.size(); ++i) {
            preconditionBDD *= set_to_bdd(precondition_sets[i], i);
        }

        // create biimps
        vector<vector<BDD>> biimp_vars_1;
        vector<vector<BDD>> biimp_vars_2;
        for (size_t i = num_precondition_sets; i < num_precondition_sets + op.get_preconditions().size(); ++i) {
            biimp_vars_1.push_back(facts_bdd_vars[i][1]);
            biimp_vars_2.push_back(facts_bdd_vars[i][1]);
        }
        biimp_vars_1.push_back(facts_bdd_vars[num_precondition_sets + num_implicit_precondition_sets][1]);
        biimp_vars_2.push_back(facts_bdd_vars[num_precondition_sets + num_implicit_precondition_sets][0]);
        BDD implicit_preconditionBDD_1 = createBiimplicationBDD(biimp_vars_1) * preconditionBDD;
        BDD implicit_preconditionBDD_2 = createBiimplicationBDD(biimp_vars_2) * preconditionBDD;

        // implicit removes
        BDD implicit_removeBDD = manager->bddZero();
        for (size_t i = 0; i < preconditions.size(); ++i) {
            for (size_t j = 0; j < implicit_removes[preconditions[i]].size(); ++j) {
                implicit_removeBDD += fact_to_bdd(implicit_removes[preconditions[i]][j], 0, 0);
            }
        }

        // create implicit preconditions
        for (size_t i = 0; i < preconditions.size(); ++i) {
            implicit_preconditionBDD_1 *= fact_to_bdd(preconditions[i], 0, num_precondition_sets + i);
            implicit_preconditionBDD_2 *= fact_to_bdd(preconditions[i], 0, num_precondition_sets + i);
        }

        // create fact effect BDD
        BDD effectBDD = manager->bddZero();
        for (size_t i = 0; i < op.get_effects().size(); ++i) {
            FactProxy fact = op.get_effects()[i].get_fact();
            effectBDD += fact_to_bdd(fact_map[make_pair(fact.get_variable().get_id(), fact.get_value())], 0, num_precondition_sets + num_implicit_precondition_sets);
        }

        // add the effects to the biimplication
        implicit_preconditionBDD_1 *= effectBDD;
        implicit_preconditionBDD_2 *= effectBDD.SwapVariables(facts_bdd_vars[num_precondition_sets + num_implicit_precondition_sets][0], facts_bdd_vars[num_precondition_sets + num_implicit_precondition_sets][1]);
        // remove the implicited deletes
        implicit_preconditionBDD_1 -= implicit_removeBDD.SwapVariables(facts_bdd_vars[0][0], facts_bdd_vars[num_precondition_sets + num_implicit_precondition_sets][1]);
        implicit_preconditionBDD_2 -= implicit_removeBDD.SwapVariables(facts_bdd_vars[0][0], facts_bdd_vars[num_precondition_sets + num_implicit_precondition_sets][0]);
        // create set effect BDD
        effectBDD *= effectBDD.SwapVariables(facts_bdd_vars[num_precondition_sets + num_implicit_precondition_sets][0], facts_bdd_vars[num_precondition_sets + num_implicit_precondition_sets][1]);
        // add to operators
        operators += (preconditionBDD * effectBDD) + (implicit_preconditionBDD_1 + implicit_preconditionBDD_2);
    }
}

}
