#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <mutex>
#include <functional>


namespace varga
{
    // base classes:
    // - Individual
    //   - stores values and fitness
    // - Polulation(Individual)
    // - Context(Populations)
    //   - holds all calculation data, including Populations
    //   - initializes first Population
    // - Runner(Context, Evaluator, Generator)
    //   - init Context
    //   - execute in loop
    //     - Evaluator.evaluate(Context)
    //     - Generator.generate(Context)

    template <typename TIndividualGenes>
    struct Individual
    {
        // define your own Individual gene(s)
        // TODO: consider using c++20 'requires' keyword
        TIndividualGenes genes;
        double fitness;
    };


    template <typename TIndividual>
    struct Population
    {
        // TODO: consider using c++20 'requires' keyword instead
        std::vector<TIndividual> individuals;
    };


    template <typename TPopulation, typename TIndividual>
    struct Context
    {
        // properties
        size_t ngeneration;
        // TODO: make this a private array and swap 2 public pointers
        TPopulation this_generation;
        TPopulation prev_generation;

        Context(TPopulation initial_population) :
            this_generation(initial_population),
            prev_generation(initial_population)
            {}
    };


    template <typename TPopulation, typename TIndividual>
    struct Evaluator
    {
        // evaluate the last population
        // modify the content of `context`
        virtual void evaluate(Context<TPopulation, TIndividual> &context)
        {
            (void) context;
        }
    };


    template <typename TPopulation, typename TIndividual>
    struct Generator
    {
        // generate new population
        // modify the content of `context`
        // return `true` if more generate() calls required
        // else return `false`
        virtual bool generate(Context<TPopulation, TIndividual> &context)
        {
            (void) context;
            std::cout << "error: abstract method used" << std::endl;
            return false;
        }
    };


    template <typename TPopulation, typename TIndividual>
    struct Runner {
        Context<TPopulation, TIndividual> context;
        Evaluator<TPopulation, TIndividual> evaluator;
        Generator<TPopulation, TIndividual> generator;

        // TODO: do we want to store pointers instead of copies?
        Runner(Context<TPopulation, TIndividual> &in_context,
               Evaluator<TPopulation, TIndividual> &in_evaluator,
               Generator<TPopulation, TIndividual> &in_generator) :
            context(in_context),
            evaluator(in_evaluator),
            generator(in_generator)
            {}

        void run()
        {
            bool continue_generating = true;
            while (continue_generating) {
                evaluator.evaluate(context);
                continue_generating = generator.generate(context);
            }
        }
    };


    // extended classes:
    // - added some utility methods

    template <typename TIndividualGenes>
    struct IndividualExt : Individual<TIndividualGenes>
    {
        // initialize the genes using rnd01() function returning random double in range [0, 1]
        virtual void from_rnd01(const std::function<double(void)> &rnd01) {
            (void) rnd01();
        }
    };


    struct Random
    {
        // return random double in range [0, 1]
        // extremely useful when producing random value over range [0, val_max):
        // val_rnd = rnd01() * val_max
        virtual double rnd01(void)
        {
            std::cout << "error: abstract method used" << std::endl;
            return 0;
        }
    };


    // TODO: try to remove public
    class RandomOpenGA : public Random
    {
        // thread-safe implementation of rnd01() borrowed from
        // https://github.com/Arash-codedev/openGA/blob/master/README.md
        // assuming those people knew what they were doing
        public:
            RandomOpenGA()
            {
                // initialize the random number generator with time-dependent seed
                uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
                std::seed_seq ss {uint32_t(timeSeed & 0xffffffff),
                                  uint32_t(timeSeed>>32)};
                rng.seed(ss);
                std::uniform_real_distribution<double> unif(0, 1);
            }

            double rnd01(void)
            {
                // prevent data race between threads
                std::lock_guard<std::mutex> lock(mtx_rand);
                return unif_dist(rng);
            }

        private:
            std::mutex mtx_rand;
            std::mt19937_64 rng;
            std::uniform_real_distribution<double> unif_dist;
    };


    template <typename TPopulation, typename TIndividual>
    struct ContextExt : Context<TPopulation, TIndividual>
    {
        // tools
        Random random = RandomOpenGA();
    };
}
