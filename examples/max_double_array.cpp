#include "varga.hpp"

size_t g_individual_ngenes = 100;

// TODO
// + MyIndividual(context?)
//   + not needed, defined via template
// + MyEvaluator(context)
//   + evaluate()
// - Generator(context)
//   - init_first_generation() <- virtual?
//   - from_rnd01() <- in generator?
//   - select_parents() <- virtual?
//   - crossover() <- virtual?
//   - mutate() <- virtual?
// - MyGenerator()
//
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

template <typename TIndividual>
struct MyEvaluator : varga::Evaluator<TIndividual>
{
    void evaluate(varga::Context<TIndividual> &context)
    {
        (void) context;
        for (auto &i : context.this_generation.individuals) {
            i.fitness = 0;
            for (auto &ig : i.genes) {
                i.fitness += ig;
            }
        }
    }
};


int main()
{
    varga::Context<my_individual_t> context;
    MyEvaluator<my_individual_t> evaluator;
    varga::Generator<my_individual_t> generator;
    varga::Runner<my_individual_t> runner{context, evaluator, generator};
    return 0;
}
