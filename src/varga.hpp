#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <mutex>
#include <functional>


namespace varga
{
    // tools

    class Random
    {
        // thread-safe implementation of rnd01() borrowed from
        // https://github.com/Arash-codedev/openGA/blob/master/README.md
        // assuming those people knew what they were doing
        private:
            std::mutex mtx_rand;
            std::mt19937_64 rng;
            std::uniform_real_distribution<double> unif_dist;

        public:
            Random()
            {
                std::cout << "init Random()" << std::endl;
                // initialize the random number generator with time-dependent seed
                uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
                std::seed_seq ss {uint32_t(timeSeed & 0xffffffff),
                                  uint32_t(timeSeed>>32)};
                rng.seed(ss);
                std::uniform_real_distribution<double> unif(0, 1);
            }

            Random(const Random &other)
            {
                // construct a new object on copy
                (void) other;
                Random();
            }

            double rnd01(void)
            {
                // prevent data race between threads
                std::lock_guard<std::mutex> lock(mtx_rand);
                return unif_dist(rng);
            }
    };


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


    template <typename TIndividual>
    struct Context
    {
        // tools
        Random random{};

        // properties
        size_t ngeneration;
        // TODO: make this a private array and swap 2 public pointers
        Population<TIndividual> this_generation;
        Population<TIndividual> prev_generation;
    };


    template <typename TIndividual>
    struct Evaluator
    {
        // evaluate the last population
        // modify the content of `context`
        virtual void evaluate(Context<TIndividual> &context)
        {
            (void) context;
        }
    };


    template <typename TIndividual>
    struct Generator
    {
        // generate new population
        // modify the content of `context`
        // return `true` if more generate() calls required
        // else return `false`
        virtual bool generate(Context<TIndividual> &context)
        {
            (void) context;
            std::cout << "error: abstract method used" << std::endl;
            return false;
        }
    };


    template <typename TIndividual>
    struct Runner {
        Context<TIndividual> context;
        Evaluator<TIndividual> evaluator;
        Generator<TIndividual> generator;

        // TODO: do we want to store pointers instead of copies?
        Runner(Context<TIndividual> &in_context,
               Evaluator<TIndividual> &in_evaluator,
               Generator<TIndividual> &in_generator) :
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

}
