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
    // - StateMachine(Context, Evaluator, Generator)
    //   - init Context
    //   - execute in loop
    //     - Evaluator.evaluate(Context)
    //     - Generator.generate(Context)

    template <typename TGenes>
    struct Individual
    {
        // virtual destructor is required if virtual methods are used
        virtual ~Individual() {}

        virtual std::string to_string(size_t n_tabs = 0) {
            (void) n_tabs;
            std::cout << "error: method not implemented" << std::endl;
            return "error: method not implemented";
        }

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

        // genetic data
        TGenes genes;
    };


    template <typename TIndividual>
    struct Population
    {
        Population(size_t in_size) : individuals(in_size), fitness(in_size) {}

        std::vector<TIndividual> individuals;
        std::vector<double> fitness;
    };


    template <typename TIndividual>
    struct Context
    {
        // TODO: properly init generations, make prev_ const
        Context(size_t in_n_generations, size_t in_population_size) :
            n_generations(in_n_generations),
            this_generation(in_population_size),
            prev_generation(in_population_size) {}

        // tools
        Random random{};

        // properties
        const size_t n_generations;
        size_t i_generation = 0;
        bool stop_state_machine = false;

        // TODO: make this a private array and swap 2 public pointers
        Population<TIndividual> this_generation;
        Population<TIndividual> prev_generation;
    };


    template <typename TIndividual>
    class StateMachine {
        private:
            typedef std::function<void(Context<TIndividual>&)> state_function_t;
            Context<TIndividual> *p_context;

        public:
            StateMachine(Context<TIndividual>& in_context) : p_context(&in_context) {}

            void run()
            {
                while (!p_context->stop_state_machine) {
                    for (state_function_t &f : state_functions) {
                        f(*p_context);
                    }
                }
            }

            std::vector<state_function_t> state_functions;
    };


    template <typename TIndividual>
    void init_first_generation(Context<TIndividual>& context)
    {
        // run only for the 0th generation
        if (context.i_generation != 0) {
            return;
        }

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
    void print_context(Context<TIndividual>& context)
    {
        std::cout << "i_geneation: " << context.i_generation << std::endl;
        std::cout << "this_generation: " << std::endl;
        for (size_t i = 0; i < context.this_generation.individuals.size(); i++) {
            std::cout
                << "\tindividuals[" << i << "]:" << std::endl
                << context.this_generation.individuals[i].to_string(2);
        }
    }

    template <typename TIndividual>
    void inc_ngeneration(Context<TIndividual>& context)
    {
        context.i_generation++;

        if (context.i_generation >= context.n_generations) {
            context.stop_state_machine = true;
            return;
        }
    }

}
