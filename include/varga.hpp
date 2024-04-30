#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <mutex>
#include <functional>
#include <cassert>
#include <ostream>
#include <fstream>
#include <iostream>
#include <iomanip>


namespace varga
{
    // settings for Context and StateMachine
    struct Settings
    {
        size_t n_generations;
        size_t population_size;
        size_t n_parents = 0;
        size_t n_keep_parents = 0;
        double p_mutation_gene = 0.0;
        double p_mutation_gene_swap = 0.0;
        double p_mutation_bad_gene = 0.0;

        size_t progress_update_period = 0;

        std::string best_fitness_log_filename{""};

        size_t best_individual_csv_creation_period = n_generations;
        std::string best_individual_filename_prefix{"best_individual_gen"};

        std::string stats_filename{"stats.txt"};

        Settings(size_t in_population_size, size_t in_n_generations) :
                n_generations(in_n_generations),
                population_size(in_population_size) {}
    };

    // tools
    std::string seconds_to_hhmmss_string(const double seconds)
    {
        std::stringstream ss{};
        ss << std::setw(2) << std::setfill('0') << ((int)seconds / 60 / 60) % 60
            << ":"
            << std::setw(2) << std::setfill('0') << ((int)seconds / 60) % 60
            << ":"
            << std::setw(2) << std::setfill('0') << (int)seconds % 60;
        return ss.str();
    }

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
            size_t n = 0;
            size_t n_max;
            size_t update_period;

            Progress(size_t in_n_max, size_t in_update_period = 0) :
                n_max(in_n_max),
                update_period(in_update_period) {}

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

                // generate a string first and then write the whole string to `os`
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
                ss << " " << std::scientific << std::setprecision(1) << get_iter_s(n_per_c) << "/s";
                double eta_s = get_eta_s(n_per_c);
                ss << " ETA " << seconds_to_hhmmss_string(eta_s);
                ss << text;
                // overwrite remalining command line with ' '
                const size_t n_chars = 10;
                for (size_t i = 0; i < n_chars; i++) {
                    ss << " ";
                }
                ss << "\r";
                os << ss.str();
                n++;

                update_last_time();
            }

            void update_last_time()
            {
                last_time = std::chrono::steady_clock::now();
            }

            double get_iter_s(const size_t n_iter)
            {
                const std::chrono::time_point<std::chrono::steady_clock> now = std::chrono::steady_clock::now();
                const double one_iter_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now-last_time).count()/n_iter;
                return 1 / (one_iter_ms/1000);
            }

            double get_eta_s(const size_t n_iter)
            {
                const std::chrono::time_point<std::chrono::steady_clock> now = std::chrono::steady_clock::now();
                const auto one_iter_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now-last_time)/n_iter;
                const size_t remaining_n_iter = n_max - n;
                const double eta_s = std::chrono::duration_cast<std::chrono::seconds>(one_iter_ms * remaining_n_iter).count();
                return eta_s;
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

        virtual void create_csv(const std::string filename)
        {
            (void) filename;
            std::cout << "error: method not implemented" << std::endl;
        }

        virtual void randomize(Settings& s, const std::function<double(void)> &rnd01)
        {
            (void) s;
            (void) rnd01();
            std::cout << "error: method not implemented" << std::endl;
        }

        virtual double get_fitness(Settings& s)
        {
            (void) s;
            std::cout << "error: method not implemented" << std::endl;
            return -1.0;
        }

        virtual void crossover(Settings& s,
                               const std::function<double(void)> &rnd01,
                               Individual<TGenes>& parent_a,
                               Individual<TGenes>& parent_b)
        {
            (void) s;
            (void) rnd01();
            (void) parent_a;
            (void) parent_b;
            std::cout << "error: method not implemented" << std::endl;
        }

        void random_mutation(Settings& s, const std::function<double(void)> &rnd01)
        {
            (void) s;
            (void) rnd01;
            std::cout << "error: method not implemented" << std::endl;
        }
    };


    template <typename TIndividual>
    struct Population
    {
        std::vector<TIndividual> individuals{};
        std::vector<double> fitness{};
        double best_fitness = 0;
        std::vector<size_t> sorted_idx{};
        std::vector<size_t> parents_idx{};

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
            Population<TIndividual> population_storage_a{};
            Population<TIndividual> population_storage_b{};

        public:
            Settings settings;
            Random random{};
            Progress progress;
            Population<TIndividual>& prev_generation;
            Population<TIndividual>& next_generation;

            // important! generations begin with 1st, not 0th
            size_t generation = 1;
            bool stop_state_machine = false;
            std::vector<double> best_fitness_log;
            std::chrono::time_point<std::chrono::steady_clock> start_time;
            std::chrono::time_point<std::chrono::steady_clock> stop_time;

            Context(Settings& s) :
                settings(s),
                progress(s.n_generations, s.progress_update_period),
                prev_generation(population_storage_a),
                next_generation(population_storage_b) {}

            void swap_generations()
            {
                if (&population_storage_a == &prev_generation) {
                    prev_generation = population_storage_b;
                    next_generation = population_storage_a;
                } else {
                    prev_generation = population_storage_a;
                    next_generation = population_storage_b;
                }

                // clean the space for the next generation
                next_generation.individuals.resize(0);
                next_generation.fitness.resize(0);
                next_generation.sorted_idx.resize(0);
                next_generation.parents_idx.resize(0);
            }

            const std::string get_best_individual_csv_filename()
            {
                std::stringstream ss;
                ss << settings.best_individual_filename_prefix << generation << ".csv";
                return ss.str();
            }

            std::string get_stats()
            {
                const double runtime_s = std::chrono::duration_cast<std::chrono::seconds>(stop_time-start_time).count();
                const double individual_s = (settings.population_size * settings.n_generations) / runtime_s;
                const double best_fitness = next_generation.best_fitness;
                std::stringstream ss{};
                ss << "generations           \t" << settings.n_generations << std::endl;
                ss << "population            \t" << settings.population_size << std::endl;
                ss << "parents               \t" << settings.n_parents << std::endl;
                ss << "keep parents          \t" << settings.n_keep_parents << std::endl;
                ss << "p(mutation gene)      \t" << settings.p_mutation_gene << std::endl;
                ss << "p(mutation bad gene)  \t" << settings.p_mutation_bad_gene << std::endl;
                ss << "p(mutation gene swap) \t" << settings.p_mutation_gene_swap << std::endl;
                ss << "---------------------"    << std::endl;
                ss << "runtime               \t" << seconds_to_hhmmss_string(runtime_s) << std::endl;
                ss << "individuals/s         \t" << individual_s << std::endl;
                ss << "best fitness          \t" << best_fitness << std::endl;
                return ss.str();
            }
    };


    template <typename TIndividual>
    class StateMachine {
        private:
            typedef std::function<void(Context<TIndividual>&)> state_function_t;
            Context<TIndividual> context;

        public:
            std::vector<state_function_t> init_functions{};
            std::vector<state_function_t> state_functions{};
            std::vector<state_function_t> closure_functions{};

            StateMachine(Settings& s) :
                context(s) {}

            void run()
            {
                context.start_time = std::chrono::steady_clock::now();
                for (state_function_t &f : init_functions) {
                    f(context);
                }
                while (!context.stop_state_machine) {
                    for (state_function_t &f : state_functions) {
                        if (context.stop_state_machine) {
                            context.stop_time = std::chrono::steady_clock::now();
                            break;
                        }
                        f(context);
                    }
                }
                context.progress.os_clean();
                for (state_function_t &f : closure_functions) {
                    f(context);
                }
            }
    };

    // state machine states
    // functions ot type StateMachinte::state_function_t

    template <typename TIndividual>
    void change_generations(Context<TIndividual>& c)
    {
        if (c.generation >= c.settings.n_generations) {
            c.stop_state_machine = true;
            return;
        }
        c.generation++;
        c.swap_generations();
    }


    template <typename TIndividual>
    void print_context(Context<TIndividual>& c)
    {
        std::cout << "geneation: " << c.generation << std::endl;
        std::cout << "next_generation: " << std::endl;
        for (size_t i = 0; i < c.next_generation.individuals.size(); i++) {
            std::cout
                << "\t[" << i << "]:" << std::endl
                << "\t\tindividuals:" << std::endl
                << c.next_generation.individuals[i].str(3)
                << "\t\tfitness: " << c.next_generation.fitness[i] << std::endl;
        }
        std::cout << "\tsorted_idx:" << std::endl;
        for (size_t i = 0; i < c.next_generation.sorted_idx.size(); i++) {
            std::cout
                << "\t\t[" << i << "]:" << c.next_generation.sorted_idx[i] << std::endl;
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
        std::cout << "geneation: " << c.generation << std::endl;
        std::cout << "\tnext_generation:" << std::endl;
        std::cout << "\t\tfitness:" << std::endl;
        for (size_t i = 0; i < c.next_generation.fitness.size(); i++) {
            std::cout
                << "\t\t\t[" << i << "]:" << c.next_generation.fitness[i] << std::endl;
        }
        std::cout << "\t\tbest_fitness: " << c.next_generation.best_fitness << std::endl;
    }

    template <typename TIndividual>
    void print_progress(Context<TIndividual>& c)
    {
        std::stringstream ss;
        ss << " best fitness " << c.next_generation.best_fitness;
        c.progress.update(std::string(ss.str()));
    }


    template <typename TIndividual>
    void print_result(Context<TIndividual>& c)
    {
        std::stringstream ss;
        size_t best_idx = c.next_generation.sorted_idx[0];
        std::cout << "best result:" << std::endl;
        std::cout << c.next_generation.individuals[best_idx].str(1) << std::endl;
    }


    template <typename TIndividual>
    void print_stats(Context<TIndividual>& c)
    {
        std::cout << c.get_stats();
    }


    template <typename TIndividual>
    void create_stats_file(Context<TIndividual>& c)
    {
        if (c.settings.stats_filename.empty()) {
            return;
        }

        std::ofstream f(c.settings.stats_filename);
        f.is_open();
        f << c.get_stats();
    }


    template <typename TIndividual>
    void create_best_fitness_log_csv(Context<TIndividual>& c)
    {
        if (c.settings.best_fitness_log_filename.empty()) {
            return;
        }

        std::ofstream f(c.settings.best_fitness_log_filename);
        f.is_open();
        f << "generation,best_fitness" << std::endl;
        for (size_t g = 0; g < c.settings.n_generations; g++) {
            f << g << "," << c.best_fitness_log[g] << std::endl;
        }
    }


    template <typename TIndividual>
    void create_best_individual_csv(Context<TIndividual>& c)
    {
        if (c.settings.best_individual_csv_creation_period != 0
        && (c.generation % c.settings.best_individual_csv_creation_period) != 0) {
            return;
        }

        TIndividual& best_individual = c.next_generation.individuals[c.next_generation.sorted_idx[0]];
        best_individual.create_csv(c.get_best_individual_csv_filename());
    }


    template <typename TIndividual>
    void randomize_next_generation(Context<TIndividual>& c)
    {
        assert (c.generation == 1);
        assert (c.next_generation.individuals.size() == 0);

        for (size_t i = 0; i < c.settings.population_size; i++) {
            TIndividual individual;
            individual.randomize(c.settings, [&c](){return c.random.rnd01();});
            c.next_generation.individuals.push_back(individual);
        }
    }


    template <typename TIndividual>
    void evaluate_next_generation(Context<TIndividual>& c)
    {
        assert(c.next_generation.individuals.size() == c.settings.population_size);
        assert(c.next_generation.fitness.size() == 0);

        // calculate fitness
        for (size_t i = 0; i < c.next_generation.individuals.size(); i++) {
            double fitness = c.next_generation.individuals[i].get_fitness(c.settings);
            c.next_generation.fitness.push_back(fitness);
        }

        c.next_generation.update_best_fitness();
        c.best_fitness_log.push_back(c.next_generation.best_fitness);
    }


    template <typename TIndividual>
    void sort_next_generation_by_fitness(Context<TIndividual>& c)
    {
        // init .sorted_idx
        assert(c.next_generation.individuals.size() == c.settings.population_size);
        assert(c.next_generation.fitness.size() == c.settings.population_size);
        assert(c.next_generation.sorted_idx.size() == 0);
        for (size_t i = 0; i < c.settings.population_size; i++) {
            c.next_generation.sorted_idx.push_back(i);
        }
        assert(c.next_generation.sorted_idx.size() == c.settings.population_size);
        // fill .sorted_idx
        std::sort(
            c.next_generation.sorted_idx.begin(),
            c.next_generation.sorted_idx.end(),
            [&c](size_t a, size_t b)->bool
            {
                return c.next_generation.fitness[a] > c.next_generation.fitness[b];
            });
    }


    template <typename TIndividual>
    void select_next_generation_parents_as_prev_generation_best(Context<TIndividual>& c)
    {
        assert(c.prev_generation.sorted_idx.size() == c.settings.population_size);
        assert(c.next_generation.parents_idx.size() == 0);
        for (size_t i = 0; i < c.settings.n_parents; i++) {
            size_t parent_i = c.prev_generation.sorted_idx[i];
            c.next_generation.parents_idx.push_back(parent_i);
        }
    }


    template <typename TIndividual>
    void add_next_generation_individuals_from_parents(Context<TIndividual>& c)
    {
        assert(c.next_generation.parents_idx.size() >= c.settings.n_keep_parents);
        assert(c.next_generation.individuals.size() == 0);
        for (size_t i = 0; i < c.settings.n_keep_parents; i++) {
            size_t parent_i = c.next_generation.parents_idx[i];
            TIndividual parent = c.prev_generation.individuals[parent_i];
            c.next_generation.individuals.push_back(parent);
        }
        assert(c.next_generation.individuals.size() == c.settings.n_keep_parents);
    }


    template <typename TIndividual>
    void add_next_generation_individuals_from_crossover(Context<TIndividual>& c)
    {
        assert(c.settings.n_parents >= 2);
        assert(c.next_generation.parents_idx.size() == c.settings.n_parents);
        assert(c.next_generation.individuals.size() == c.settings.n_keep_parents);
        for (size_t i = c.settings.n_keep_parents; i < c.settings.population_size; i++) {
            size_t parent_a_i = c.settings.n_parents * c.random.rnd01();
            size_t parent_b_i;
            do {
                parent_b_i = c.settings.n_parents * c.random.rnd01();
            } while (parent_b_i == parent_a_i);
            TIndividual& parent_a = c.prev_generation.individuals[c.next_generation.parents_idx[parent_a_i]];
            TIndividual& parent_b = c.prev_generation.individuals[c.next_generation.parents_idx[parent_b_i]];
            TIndividual child{};
            child.crossover(c.settings,
                            [&c](){return c.random.rnd01();},
                            parent_a,
                            parent_b);
            c.next_generation.individuals.push_back(child);
        }
        assert(c.next_generation.individuals.size() == c.prev_generation.individuals.size());
    }


    template <typename TIndividual>
    void next_generation_random_mutation(Context<TIndividual>& c)
    {
        for (auto &ind : c.next_generation.individuals) {
            ind.random_mutation(c.settings, [&c](){return c.random.rnd01();});
        }
    }
}
