#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <mutex>
#include <functional>
#include <cassert>
#include <ostream>
#include <iostream>
#include <iomanip>


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

    class Progress
    {
        private:
            std::chrono::time_point<std::chrono::steady_clock> last_time = std::chrono::steady_clock::now();

        public:
            std::ostream& os = std::cerr;
            char c_opening_bracket = '[';
            char c_closing_bracket = ']';
            char c_fill = '.';
            char c_no_fill = ' ';
            size_t bar_len = 20;
            size_t n_max;
            size_t n = 0;
            size_t update_period = 0;

            Progress(size_t in_n_max) : n_max(in_n_max) {}

            ~Progress()
            {
                os_clean();
            }

            void update()
            {
                update("");
            }

            void update(std::string text)
            {
                assert(n <= n_max);
                // to not overload the console - update only at update_period
                // or once per every progress bar step
                size_t n_per_c;
                if (update_period) {
                    n_per_c = update_period;
                } else {
                    n_per_c = n_max/bar_len;
                }
                if (n % n_per_c != 0) {
                    n++;
                    return;
                }

                // generate a string first and then writhe the whole string to `os`
                // to prevent blinking cursor from jumping all over the place
                std::stringstream ss;

                ss << c_opening_bracket;
                size_t n_fill = (double)n/n_max * bar_len;
                for (size_t i = 0; i < bar_len; i++) {
                    if (i < n_fill) {
                        ss << c_fill;
                    } else {
                        ss << c_no_fill;
                    }
                }
                ss << c_closing_bracket;
                ss << " " << std::fixed << std::setprecision(1) << (double)n/n_max * 100 << "%";
                ss << " (" << std::scientific << std::setprecision(1) << iter_s(n_per_c) << " iter/s)";
                if (text != "") {
                    ss << ", " << text;
                }
                ss << "\r";
                os << ss.str();
                n++;
            }

            double iter_s(size_t n_iter)
            {
                assert(n_iter > 0);
                std::chrono::time_point<std::chrono::steady_clock> now = std::chrono::steady_clock::now();
                double duration = std::chrono::duration_cast<std::chrono::milliseconds>(now-last_time).count();
                last_time = now;
                return n_iter / (duration/1000);
            }

            void os_clean(size_t n_chars = 100)
            {
                // overwrite command line with n_chars ' '
                for (size_t i = 0; i < n_chars; i++) {
                    os << " ";
                }
                os << "\r" << std::flush;
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
        TGenes genes;

        // virtual destructor is required if virtual methods are used
        virtual ~Individual() {}

        virtual std::string str(size_t n_tabs = 0)
        {
            (void) n_tabs;
            std::cout << "error: method not implemented" << std::endl;
            return "error: method not implemented";
        }

        virtual void randomize(const std::function<double(void)> &rnd01)
        {
            (void) rnd01();
            std::cout << "error: method not implemented" << std::endl;
        }

        virtual double get_fitness(void)
        {
            std::cout << "error: method not implemented" << std::endl;
            return -1.0;
        }

        virtual void crossover(const std::function<double(void)> &rnd01,
                               Individual<TGenes>& parent_a,
                               Individual<TGenes>& parent_b)
        {
            (void) rnd01();
            (void) parent_a;
            (void) parent_b;
            std::cout << "error: method not implemented" << std::endl;
        }

        void random_mutation(const std::function<double(void)> &rnd01)
        {
            (void) rnd01();
            std::cout << "error: method not implemented" << std::endl;
        }
    };


    template <typename TIndividual>
    struct Population
    {
        std::vector<TIndividual> individuals;
        std::vector<double> fitness;
        double best_fitness;
        std::vector<size_t> sorted_idx;
        std::vector<size_t> parents_idx;

        Population(size_t in_size) :
            individuals(in_size),
            fitness(in_size),
            best_fitness(0),
            sorted_idx(in_size),
            parents_idx(0) {}
        
        void update_best_fitness()
        {
            best_fitness = *std::max_element(fitness.begin(),
                                             fitness.end());
        }
    };


    template <typename TIndividual>
    class Context
    {
        private:
            Population<TIndividual> population_storage_a;
            Population<TIndividual> population_storage_b;

        public:
            Random random{};
            Progress progress;
            Population<TIndividual>& prev_generation;
            Population<TIndividual>& next_generation;

            size_t n_generations;
            size_t n_parents;
            size_t keep_n_parents;
            double p_mutation;

            size_t i_generation;
            bool stop_state_machine;

            Context(size_t in_population_size, size_t in_n_generations) :
                population_storage_a(in_population_size),
                population_storage_b(in_population_size),
                progress(in_n_generations),
                prev_generation(population_storage_a),
                next_generation(population_storage_b),
                n_generations(in_n_generations),
                n_parents(0),
                keep_n_parents(0),
                p_mutation(0.0),
                i_generation(0),
                stop_state_machine(false)
            {
                // clean the space for the next generation
                next_generation.parents_idx.resize(0);
                next_generation.individuals.resize(0);
            }

            void swap_generations()
            {
                if (&population_storage_a == &prev_generation) {
                    prev_generation = population_storage_b;
                    next_generation = population_storage_a;
                } else {
                    prev_generation = population_storage_a;
                    next_generation = population_storage_b;
                }
            }
    };


    template <typename TIndividual>
    class StateMachine {
        private:
            typedef std::function<void(Context<TIndividual>&)> state_function_t;
            Context<TIndividual> *p_context;

        public:
            std::vector<state_function_t> init_functions;
            std::vector<state_function_t> state_functions;
            std::vector<state_function_t> closure_functions;

            StateMachine(Context<TIndividual>& in_context) :
                p_context(&in_context),
                init_functions(0),
                state_functions(0),
                closure_functions(0) {}

            void run()
            {
                for (state_function_t &f : init_functions) {
                    f(*p_context);
                }
                if (state_functions.size() > 0) {
                    while (!p_context->stop_state_machine) {
                        for (state_function_t &f : state_functions) {
                            f(*p_context);
                        }
                    }
                }
                p_context->progress.os_clean();
                for (state_function_t &f : closure_functions) {
                    f(*p_context);
                }
            }
    };

    // state machine states
    // functions ot type StateMachinte::state_function_t

    template <typename TIndividual>
    void change_generations(Context<TIndividual>& c)
    {
        c.i_generation++;
        if (c.i_generation >= c.n_generations) {
            c.stop_state_machine = true;
            return;
        }

        c.swap_generations();

        // clean the space for the next generation
        c.next_generation.parents_idx.resize(0);
        c.next_generation.individuals.resize(0);
    }


    template <typename TIndividual>
    void print_context(Context<TIndividual>& c)
    {
        std::cout << "i_geneation: " << c.i_generation << std::endl;
        std::cout << "prev_generation: " << std::endl;
        for (size_t i = 0; i < c.prev_generation.individuals.size(); i++) {
            std::cout
                << "\t[" << i << "]:" << std::endl
                << "\t\tindividuals:" << std::endl
                << c.prev_generation.individuals[i].str(3)
                << "\t\tfitness: " << c.prev_generation.fitness[i] << std::endl;
        }
        std::cout << "\tsorted_idx:" << std::endl;
        for (size_t i = 0; i < c.prev_generation.sorted_idx.size(); i++) {
            std::cout
                << "\t\t[" << i << "]:" << c.prev_generation.sorted_idx[i] << std::endl;
        }
        std::cout << "next_generation: " << std::endl;
        for (size_t i = 0; i < c.next_generation.individuals.size(); i++) {
            std::cout
                << "\t[" << i << "]:" << std::endl
                << "\t\tindividuals:" << std::endl
                << c.next_generation.individuals[i].str(3);
        }
        std::cout << "\tparents_idx:" << std::endl;
        for (size_t i = 0; i < c.next_generation.parents_idx.size(); i++) {
            std::cout
                << "\t\t[" << i << "]:" << c.next_generation.parents_idx[i] << std::endl;
        }
    }

    template <typename TIndividual>
    void print_fitness(Context<TIndividual>& c)
    {
        std::cout << "i_geneation: " << c.i_generation << std::endl;
        std::cout << "\tprev_generation:" << std::endl;
        std::cout << "\t\tfitness:" << std::endl;
        for (size_t i = 0; i < c.prev_generation.fitness.size(); i++) {
            std::cout
                << "\t\t\t[" << i << "]:" << c.prev_generation.fitness[i] << std::endl;
        }
        std::cout << "\t\tbest_fitness: " << c.prev_generation.best_fitness << std::endl;
    }

    template <typename TIndividual>
    void print_progress(Context<TIndividual>& c)
    {
        std::stringstream ss;
        ss << "best fitness: " << c.prev_generation.best_fitness;
        c.progress.update(std::string(ss.str()));
    }


    template <typename TIndividual>
    void print_result(Context<TIndividual>& c)
    {
        std::stringstream ss;
        size_t best_idx = c.prev_generation.sorted_idx[0];
        std::cout << "best result:" << std::endl;
        std::cout << c.prev_generation.individuals[best_idx].str(1) << std::endl;
        std::cout << "best fitness: " << c.prev_generation.best_fitness << std::endl;
    }


    template <typename TIndividual>
    void randomize_prev_generation(Context<TIndividual>& c)
    {
        assert (c.i_generation == 0);

        for (auto &i : c.prev_generation.individuals) {
            i.randomize([&c](){return c.random.rnd01();});
        }
    }


    template <typename TIndividual>
    void evaluate_prev_generation(Context<TIndividual>& c)
    {
        assert(c.prev_generation.individuals.size() != 0);

        // calculate fitness
        for (size_t i = 0; i < c.prev_generation.individuals.size(); i++) {
            c.prev_generation.fitness[i] = \
                c.prev_generation.individuals[i].get_fitness();
        }

        c.prev_generation.update_best_fitness();
    }


    template <typename TIndividual>
    void select_parents_as_most_fit(Context<TIndividual>& c)
    {
        // init .sorted_idx
        size_t population_size = c.prev_generation.individuals.size();
        assert(population_size != 0);
        assert(c.prev_generation.sorted_idx.size() == population_size);
        assert(c.prev_generation.fitness.size() == population_size);
        for (size_t i = 0; i < population_size; i++) {
            c.prev_generation.sorted_idx[i] = i;
        }

        // fill .sorted_idx
        std::sort(
            c.prev_generation.sorted_idx.begin(),
            c.prev_generation.sorted_idx.end(),
            [&c](size_t a, size_t b)->bool
            {
                return c.prev_generation.fitness[a] > c.prev_generation.fitness[b];
            });

        // select parents
        assert(c.next_generation.parents_idx.size() == 0);
        for (size_t i = 0; i < c.n_parents; i++) {
            size_t parent_i = c.prev_generation.sorted_idx[i];
            c.next_generation.parents_idx.push_back(parent_i);
        }
    }


    template <typename TIndividual>
    void move_parents_to_next_generation(Context<TIndividual>& c)
    {
        assert(c.next_generation.parents_idx.size() >= c.keep_n_parents);
        assert(c.next_generation.individuals.size() == 0);
        for (size_t i = 0; i < c.keep_n_parents; i++) {
            size_t parent_i = c.next_generation.parents_idx[i];
            TIndividual parent = c.prev_generation.individuals[parent_i];
            c.next_generation.individuals.push_back(parent);
        }
        assert(c.next_generation.individuals.size() == c.keep_n_parents);
    }


    template <typename TIndividual>
    void create_children_from_crossover(Context<TIndividual>& c)
    {
        assert(c.n_parents >= 2);
        assert(c.next_generation.parents_idx.size() == c.n_parents);
        assert(c.next_generation.individuals.size() == c.keep_n_parents);
        const size_t population_size = c.prev_generation.individuals.size();
        for (size_t i = c.keep_n_parents; i < population_size; i++) {
            size_t parent_a_i = c.n_parents * c.random.rnd01();
            size_t parent_b_i;
            do {
                parent_b_i = c.n_parents * c.random.rnd01();
            } while (parent_b_i == parent_a_i);
            TIndividual& parent_a = c.prev_generation.individuals[c.next_generation.parents_idx[parent_a_i]];
            TIndividual& parent_b = c.prev_generation.individuals[c.next_generation.parents_idx[parent_b_i]];
            TIndividual child{};
            child.crossover([&c](){return c.random.rnd01();},
                            parent_a,
                            parent_b);
            c.next_generation.individuals.push_back(child);
        }
        assert(c.next_generation.individuals.size() == c.prev_generation.individuals.size());
    }


    template <typename TIndividual>
    void random_mutation(Context<TIndividual>& c)
    {
        for (auto &ind : c.next_generation.individuals) {
            if (c.random.rnd01() < c.p_mutation) {
                ind.random_mutation([&c](){return c.random.rnd01();});
            }
        }
    }
}
