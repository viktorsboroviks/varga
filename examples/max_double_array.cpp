#include "varga.hpp"

size_t g_individual_ngenes = 100;

//struct MyIndividual : varga::Individual
//{
//    std::vector<double> genes;
//
//    void from_rnd01(const std::function<double(void)> &rnd01)
//    {
//        for (auto &g : genes) {
//            g = rnd01();
//        }
//    }
//
//};
//
//struct MyPopulation : varga::Population
//{
//    std::vector<MyIndividual> individuals;
//
//    void from_rnd01(const std::function<double(void)> &rnd01)
//    {
//        for (auto &i : individuals) {
//            i.from_rnd01(rnd01);
//        }
//    }
//
//    MyPopulation(size_t ngenes)
//    {
//        individuals.resize(ngenes);
//    }
//};
//
//struct MyEvaluator : varga::Evaluator
//{
//    void evaluate(varga::ContextExt<MyPopulation> &context)
//    {
//        for (MyIndividual &i : context.this_generation.individuals) {
//            i.fitness = 0;
//            for (auto &iv : i.value) {
//                i.fitness += i.value;
//            }
//        }
//    }
//};

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

//typedef varga::IndividualExt<std::vector<double>> individual_t;
//typedef std::vector<individual_t> population_t;

int main()
{
//    MyPopulation initial_population;
//    varga::ContextExt context(initial_population);
//    varga::Evaluator evaluator;
//    varga::Generator generator;
//    varga::Runner runner{context, evaluator, generator};
    return 0;
}
