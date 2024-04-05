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

typedef varga::Individual<std::vector<double>> my_individual_t;

template <typename TGenes>
class MyIndividual : varga::Individual<TGenes>
{
    void randomize(void) {}
        virtual void single_point_crossover(void) {}
        virtual void evaluate_fitness(void) {}
        virtual void random_mutation(void) {}
};

int main()
{
    varga::Context<my_individual_t> context;
    varga::Runner<my_individual_t> runner{context};
    return 0;
}
