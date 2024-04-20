#include <string>
#include <sstream>
#include <fstream>
#include "varga.hpp"

const size_t g_n_generations = 1000000;
const size_t g_population_size = 10;
const size_t g_individual_n_genes = 100;
const size_t g_n_parents = 3;
const double g_p_mutation = 0.05;


struct MyIndividual : varga::Individual<std::vector<double>>
{
    MyIndividual() {
        genes.resize(g_individual_n_genes);
    }

    std::string str(size_t n_tabs = 0)
    {
        std::stringstream ss;
        for (size_t n = 0; n < n_tabs; n++) {
            ss << "\t";
        }
        ss << "genes:" << std::endl;
        for (size_t i = 0; i < genes.size(); i++) {
            for (size_t n = 0; n < n_tabs; n++) {
                ss << "\t";
            }
            ss << "\t[" << i << "]:" << genes[i];
            if ((i + 1) != genes.size()) {
                ss << std::endl;
            }
        }
        return ss.str();
    }

    void randomize(const std::function<double(void)> &rnd01)
    {
        for (auto& g : genes) {
            g = rnd01();
        }
    }

    double get_fitness(void)
    {
        double fitness = 0;
        for (auto& g : genes) {
            fitness += g;
        }
        return fitness;
    }

    void crossover(const std::function<double(void)> &rnd01,
                   Individual<std::vector<double>> &parent_a,
                   Individual<std::vector<double>> &parent_b)
    {
        assert(genes.size() != 0);
        assert(parent_a.genes.size() != 0);
        assert(parent_b.genes.size() != 0);
        assert(genes.size() == parent_a.genes.size());
        assert(parent_a.genes.size() == parent_b.genes.size());
        size_t crossover_i = rnd01() * genes.size();
        for (size_t i = 0; i < genes.size(); i++) {
            if (i < crossover_i) {
                genes[i] = parent_a.genes[i];
            } else {
                genes[i] = parent_b.genes[i];
            }
        }
    }

    void random_mutation(const std::function<double(void)> &rnd01)
    {
        assert(genes.size() != 0);
        size_t mutation_i = rnd01() * genes.size();
        genes[mutation_i] = rnd01();
    }
};


int main()
{
    varga::Context<MyIndividual> c{g_population_size, g_n_generations};
    c.n_parents = g_n_parents;
    c.p_mutation = g_p_mutation;

    varga::StateMachine<MyIndividual> sm{c};
    sm.init_functions = {varga::randomize_next_generation<MyIndividual>};
    sm.state_functions = {varga::evaluate_next_generation<MyIndividual>,
                          varga::sort_next_generation_by_fitness<MyIndividual>,
//                          varga::print_fitness<MyIndividual>,
//                          varga::print_context<MyIndividual>,
                          varga::print_progress<MyIndividual>,
                          varga::change_generations<MyIndividual>,
                          varga::select_next_generation_parents_as_prev_generation_best<MyIndividual>,
                          varga::add_next_generation_individuals_from_parents<MyIndividual>,
                          varga::add_next_generation_individuals_from_crossover<MyIndividual>,
                          varga::next_generation_random_mutation<MyIndividual>};
    sm.closure_functions = {varga::print_result<MyIndividual>,
                            varga::create_best_fitness_log_csv<MyIndividual>};
    sm.run();
    return 0;
}
