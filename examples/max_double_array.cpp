#include "varga.hpp"

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

struct MyIndividual : varga::Individual<std::vector<double>>
{
    void randomize(const std::function<double(void)> &rnd01)
    {
        for (auto& g : genes) {
            g = rnd01();
        }
    }
    double get_fitness(void) {
        double fitness = 0;
        for (auto& g : genes) {
            fitness += g;
        }
        return fitness;
    }
};

int main()
{
    varga::Context<MyIndividual> context;
    varga::Runner<MyIndividual> runner{context};
    return 0;
}
