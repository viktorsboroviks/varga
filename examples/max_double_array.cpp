#include <string>
#include <sstream>
#include "varga.hpp"

const size_t g_n_generations = 100;
const size_t g_population_size = 10;
const size_t g_individual_n_genes = 2;
const size_t g_n_parents_mating = 3;


struct MyIndividual : varga::Individual<std::vector<double>>
{
    MyIndividual() {
        genes.resize(g_individual_n_genes);
    }

    std::string to_string(size_t n_tabs = 0)
    {
        std::stringstream ss;
        for (size_t i = 0; i < genes.size(); i++) {
            for (size_t n = 0; n < n_tabs; n++) {
                ss << "\t";
            }
            ss
                << "genes[" << i << "]:"
                << genes[i] << std::endl;
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
};


// TODO
// - Parent selection = steady state selection
//     - num_parents_mating = 3
//     - population = 10
// - Crossover = single point
//     - for every child
//     - select random crossover point
//     - select random 2 parents from set
//     - crossover
// - Mutation = random
int main()
{
    varga::Context<MyIndividual> c{g_population_size};
    c.n_generations = g_n_generations;
    c.n_parents_mating = g_n_parents_mating;

    varga::StateMachine<MyIndividual> sm{c};
    sm.init_functions = {varga::randomize_prev_generation<MyIndividual>};
    sm.state_functions = {varga::evaluate<MyIndividual>,
                          varga::select_best_parents<MyIndividual>,
                          varga::move_parents_to_next_generation<MyIndividual>,
//                          varga::single_point_crossover<MyIndividual>,
//                          varga::random_mutation<MyIndividual>,
//                          varga::print_context<MyIndividual>,
                          varga::change_generation<MyIndividual>};
    sm.run();
    return 0;
}
