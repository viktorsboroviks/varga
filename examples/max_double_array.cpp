#include "varga.hpp"

// TODO
// - Individual = vector<double>, size=parameter
// - Parent selection = steady state selection
//     - num_parents_mating = 3
//     - population = 10
// - Crossover = single point
//     - for every child
//     - select random crossover point
//     - select random 2 parents from set
//     - crossover
// - Mutation = random

typedef varga::IndividualExt<std::vector<double>> individual_t;

int main()
{
    varga::ContextExt<individual_t> context;
    varga::Evaluator<individual_t> evaluator;
    varga::Generator<individual_t> generator;
    varga::Runner<individual_t> runner{context, evaluator, generator};
    return 0;
}
