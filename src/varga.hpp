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

    template <typename TIndividualValue>
    struct Individual
    {
        TIndividualValue value;
        double fitness;
    };


    template <typename TIndividualType>
    struct Population
    {
        std::vector<TIndividualType> individuals;
    };


    template <typename TIndividualType>
    struct Context
    {
        size_t ngeneration;

        Population<TIndividualType> prev_generation;
        Population<TIndividualType> this_generation;

        // init of this_generation for the first generation
        // is also intended to happen in this class.
    };


    template <typename TIndividualType>
    struct Evaluator
    {
        // evaluate the last population
        // modify the content of `context`
        virtual void evaluate(Context<TIndividualType> &context)
        {
            (void) context;
        }
    };


    template <typename TIndividualType>
    struct Generator
    {
        // generate new population
        // modify the content of `context`
        // return `true` if more generate() calls required
        // else return `false`
        virtual bool generate(Context<TIndividualType> &context)
        {
            (void) context;
            std::cout << "error: abstract method used" << std::endl;
            return false;
        }
    };


    template <typename TIndividualType>
    struct Runner {
        Context<TIndividualType> context;
        Evaluator<TIndividualType> evaluator;
        Generator<TIndividualType> generator;

        // TODO: do we want to store pointers instead of copies?
        Runner(Context<TIndividualType> &in_context,
               Evaluator<TIndividualType> &in_evaluator,
               Generator<TIndividualType> &in_generator) :
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

    template <typename TIndividualValue>
    struct IndividualExt : Individual<TIndividualValue>
    {
        TIndividualValue value;
        double fitness;

        // initialize the `value` using rnd01() function returning random double in range [0, 1]
        void from_rnd01(const std::function<double(void)> &rnd01)
        {
            value = rnd01();
        };
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


    template <typename TIndividualType>
    struct ContextExt : Context<TIndividualType>
    {
        // tools
        Random random;

        // properties
        size_t ngeneration;
        Population<TIndividualType> prev_generation;
        Population<TIndividualType> this_generation;

        // TODO: ask if this can be made something like
        // ContextExt(Random in_random = RandomOpenGA()) : random(in_random) {}
        ContextExt() : random(RandomOpenGA()) {}
        ContextExt(Random in_random) : random(in_random) {}
    };
}
