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

    template <typename TGenes>
    struct Individual
    {
        // virtual destructor is required if virtual methods are used
        virtual ~Individual() {}

        virtual void randomize(const std::function<double(void)> &rnd01)
        {
            (void) rnd01();
            std::cout << "error: method not implemented" << std::endl;
        }
        virtual double get_fitness(void) {
            std::cout << "error: method not implemented" << std::endl;
            return -1.0;
        }
        virtual void single_point_crossover(void) {}
        virtual void random_mutation(void) {}

        TGenes genes;
    };


    template <typename TIndividual>
    struct Population
    {
        std::vector<TIndividual> individuals;
        std::vector<double> fitness;
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
    void init_first_generation(Context<TIndividual>& context)
    {
        for (auto &i : context.this_generation.individuals) {
            i.randomize([&context](){return context.random.rnd01();});
        }
    }

    template <typename TIndividual>
    void evaluate(Context<TIndividual>& context)
    {
        // calculate fitness
        for (size_t i = 0; i < context.this_generation.individuals.size(); i++) {
            context.this_generation.fitness[i] = \
                context.this_generation.individuals[i].get_fitness();
        }
    }

    template <typename TIndividual>
    class Runner {
        private:
            Context<TIndividual> *p_context;

        public:
            typedef std::function<void(Context<TIndividual>&)> function_t;

            Runner(Context<TIndividual>& in_context) : p_context(&in_context) {}

            function_t f_init_first_generation = init_first_generation<TIndividual>;
            function_t f_evaluate = evaluate<TIndividual>;
            function_t f_select_parents = nullptr;
            function_t f_crossover = nullptr;
            function_t f_mutate = nullptr;
            function_t f_generate = nullptr;

            void run()
            {
                f_init_first_generation(*p_context);
                while (p_context->continue_generating) {
                    f_evaluate(*p_context);
                    f_select_parents(*p_context);
                    f_crossover(*p_context);
                    f_mutate(*p_context);
                    f_generate(*p_context);
                }
            }
    };
}
