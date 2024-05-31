#include <cassert>
#include <chrono>
#include <deque>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <ostream>
#include <random>
#include <vector>

namespace varga {

// settings for Context and StateMachine
struct Settings {
    std::map<std::string, double> custom_parameter;

    size_t n_generations;
    size_t population_size;
    size_t n_elite_best = 0;
    size_t n_elite_worst = 0;
    size_t n_elite_random = 0;
    size_t n_parents_best = 0;
    size_t n_parents_worst = 0;
    size_t n_parents_random = 0;
    size_t n_parents_randomized = 0;

    double p_replace_solution = 0.0;
    double p_replace_data = 0.0;
    double p_mutate_data = 0.0;
    double p_mutate_bad_data = 0.0;
    double p_swap_data = 0.0;

    size_t progress_update_period = 0;

    std::string log_filename{""};

    size_t best_solution_csv_creation_period = n_generations;
    std::string best_solution_filename_prefix{"best_solution_gen"};

    std::string stats_filename{"stats.txt"};

    Settings(size_t in_population_size, size_t in_n_generations) :
        n_generations(in_n_generations),
        population_size(in_population_size)
    {
    }
};

// tools
std::string seconds_to_hhmmss_string(const double seconds)
{
    std::stringstream ss{};
    ss << std::setw(2) << std::setfill('0') << ((int)seconds / 60 / 60) % 60;
    ss << ":" << std::setw(2) << std::setfill('0') << ((int)seconds / 60) % 60;
    ss << ":" << std::setw(2) << std::setfill('0') << (int)seconds % 60;
    return ss.str();
}

class Random {
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
        uint64_t timeSeed = std::chrono::high_resolution_clock::now()
                                    .time_since_epoch()
                                    .count();
        std::seed_seq ss{uint32_t(timeSeed & 0xffffffff),
                         uint32_t(timeSeed >> 32)};
        rng.seed(ss);
        std::uniform_real_distribution<double> unif(0, 1);
    }

    Random(const Random& other)
    {
        // construct a new object on copy
        (void)other;
        Random();
    }

    double rnd01(void)
    {
        // prevent data race between threads
        std::lock_guard<std::mutex> lock(mtx_rand);
        return unif_dist(rng);
    }
};

class Progress {
private:
    std::chrono::time_point<std::chrono::steady_clock> last_update_time =
            std::chrono::steady_clock::now();

public:
    std::ostream& os = std::cerr;
    char c_opening_bracket = '[';
    char c_closing_bracket = ']';
    char c_fill = '.';
    char c_no_fill = ' ';
    size_t bar_len = 10;
    size_t n = 0;
    size_t n_max;
    size_t update_period;

    Progress(size_t in_n_max, size_t in_update_period = 0) :
        n_max(in_n_max),
        update_period(in_update_period)
    {
    }

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
        }
        else {
            n_per_c = n_max / bar_len;
        }
        if (n % n_per_c != 0) {
            n++;
            return;
        }

        // generate a string first and then write the whole string to `os`
        // to prevent blinking cursor from jumping all over the place
        std::stringstream ss;

        ss << c_opening_bracket;
        size_t n_fill = (double)n / n_max * bar_len;
        for (size_t i = 0; i < bar_len; i++) {
            if (i < n_fill) {
                ss << c_fill;
            }
            else {
                ss << c_no_fill;
            }
        }

        const size_t n_max_strlen = std::to_string(n_max).length();

        ss << c_closing_bracket;
        ss << " " << std::setfill('0') << std::setw(n_max_strlen) << n;
        ss << "/" << n_max;
        ss << " " << std::fixed << std::setprecision(1)
           << (double)n / n_max * 100 << "%";
        double eta_s = get_eta_s(n_per_c);
        ss << " ETA " << seconds_to_hhmmss_string(eta_s);
        ss << text;
        // overwrite remalining command line with ' '
        const size_t n_chars = 5;
        for (size_t i = 0; i < n_chars; i++) {
            ss << " ";
        }
        ss << "\r";
        os << ss.str();
        n++;

        last_update_time = std::chrono::steady_clock::now();
    }

    double get_iter_s(const size_t n_iter)
    {
        const std::chrono::time_point<std::chrono::steady_clock> now =
                std::chrono::steady_clock::now();
        const double us_per_n_iter =
                std::chrono::duration_cast<std::chrono::microseconds>(
                        now - last_update_time)
                        .count();
        return n_iter / us_per_n_iter * 1000000;
    }

    double get_eta_s(const size_t n_iter)
    {
        const std::chrono::time_point<std::chrono::steady_clock> now =
                std::chrono::steady_clock::now();
        const double us_per_n_iter =
                std::chrono::duration_cast<std::chrono::microseconds>(
                        now - last_update_time)
                        .count();
        const double remaining_n_iter = n_max - n;
        const double eta_s =
                us_per_n_iter / n_iter * remaining_n_iter / 1000000;
        return eta_s;
    }

    void os_clean(size_t n_chars = 200)
    {
        // overwrite command line with n_chars ' '
        for (size_t i = 0; i < n_chars; i++) {
            os << " ";
        }
        os << "\r" << std::flush;
    }
};

// base classes:
// - Solution
//   - stores values and value
// - Polulation(Solution)
// - Context(Populations)
//   - holds all calculation data, including Populations
//   - initializes first Population
// - StateMachine(Context, Evaluator, Generator)
//   - init Context
//   - execute in loop
//     - Evaluator.evaluate(Context)
//     - Generator.generate(Context)

template <typename TData>
class Solution {
protected:
    // set to true every time solution changes
    bool _changed = true;
    double _value = 0.0;

public:
    TData data;

    // virtual destructor is required if virtual methods are used
    virtual ~Solution() {}

    virtual std::string str(size_t n_tabs = 0)
    {
        (void)n_tabs;
        std::cout << "error: method not implemented" << std::endl;
        return "error: method not implemented";
    }

    virtual void create_csv(const std::string filename)
    {
        (void)filename;
        std::cout << "error: create_csv not implemented" << std::endl;
    }

    virtual void randomize(varga::Settings& s,
                           const std::function<double(void)>& rnd01)
    {
        (void)s;
        (void)rnd01();
        _changed = true;
        std::cout << "error: randomize method not implemented" << std::endl;
    }

    virtual double get_value(Settings& s)
    {
        (void)s;
        std::cout << "error: get_value method not implemented" << std::endl;
        return -1.0;
    }

    template <typename TSolution>
    void crossover(Settings& s,
                   const std::function<double(void)>& rnd01,
                   TSolution& parent_a,
                   TSolution& parent_b)
    {
        (void)s;
        (void)rnd01();
        (void)parent_a;
        (void)parent_b;
        std::cout << "error: crossover method not implemented" << std::endl;
        std::cout << "hint: you can call one of existing crossover functions"
                  << std::endl;
    }

    template <typename TSolution>
    void uniform_crossover(Settings& s,
                           const std::function<double(void)>& rnd01,
                           TSolution& parent_a,
                           TSolution& parent_b)
    {
        (void)s;
        assert(data.size() != 0);
        assert(parent_a.data.size() != 0);
        assert(parent_b.data.size() != 0);
        assert(data.size() == parent_a.data.size());
        assert(parent_a.data.size() == parent_b.data.size());
        for (size_t i = 0; i < data.size(); i++) {
            if (rnd01() < 0.5) {
                data[i] = parent_a.data[i];
            }
            else {
                data[i] = parent_b.data[i];
            }
        }
        _changed = true;
    }

    template <typename TSolution>
    void one_point_crossover(Settings& s,
                             const std::function<double(void)>& rnd01,
                             TSolution& parent_a,
                             TSolution& parent_b)
    {
        (void)s;

        assert(data.size() != 0);
        assert(parent_a.data.size() != 0);
        assert(parent_b.data.size() != 0);
        assert(data.size() == parent_a.data.size());
        assert(parent_a.data.size() == parent_b.data.size());

        const size_t point = rnd01() * data.size();
        for (size_t i = 0; i < data.size(); i++) {
            if (i < point)
                data[i] = parent_a.data[i];
            else
                data[i] = parent_b.data[i];
        }
        _changed = true;
    }

    template <typename TSolution>
    void two_point_crossover(Settings& s,
                             const std::function<double(void)>& rnd01,
                             TSolution& parent_a,
                             TSolution& parent_b)
    {
        (void)s;

        assert(data.size() != 0);
        assert(parent_a.data.size() != 0);
        assert(parent_b.data.size() != 0);
        assert(data.size() == parent_a.data.size());
        assert(parent_a.data.size() == parent_b.data.size());
        size_t point_a;
        do {
            point_a = rnd01() * data.size();
            if (point_a < (data.size() - 1))
                break;
        } while (true);
        size_t point_b;
        do {
            point_b = rnd01() * data.size();
            if (point_b > point_a)
                break;
        } while (true);

        for (size_t i = 0; i < data.size(); i++) {
            if (i < point_a || i > point_b)
                data[i] = parent_a.data[i];
            else
                data[i] = parent_b.data[i];
        }
        _changed = true;
    }

    template <typename TSolution>
    void replace(Settings& s,
                 const std::function<double(void)>& rnd01,
                 const std::vector<TSolution>& all_solutions)
    {
        // replace solution (all data)
        if (rnd01() < s.p_replace_solution) {
            const size_t src_idx = rnd01() * all_solutions.size();
            for (size_t i = 0; i < data.size(); i++) {
                data[i] = all_solutions[src_idx].data[i];
            }
            _changed = true;
            return;
        }

        // replace data
        for (auto& d : data) {
            if (rnd01() < s.p_replace_data) {
                const size_t src_idx = rnd01() * data.size();
                d = data[src_idx];
                _changed = true;
            }
        }
    }

    template <typename TSolution>
    void swap(Settings& s,
              const std::function<double(void)>& rnd01,
              const std::vector<TSolution>& all_solutions)
    {
        // follows the same approach as replace() for simplicity

        // swap solutions
        // there is currently no use case for this, but it might be helpful
        // in some configurations in the future.
        (void)all_solutions;

        // swap data
        if (rnd01() < s.p_swap_data) {
            const size_t i1 = rnd01() * data.size();
            const size_t i2 = rnd01() * data.size();
            const auto store_data = data[i1];
            data[i1] = data[i2];
            data[i2] = store_data;
            _changed = true;
        }
    }

    virtual void mutate(Settings& s, const std::function<double(void)>& rnd01)
    {
        // this method is very solution-specific, so to not overthink it
        // I leave it virtual
        (void)s;
        (void)rnd01;
        _changed = true;
        std::cout << "error: mutate method not implemented" << std::endl;
    }
};

template <typename TSolution>
struct Population {
    std::vector<TSolution> solutions{};
    double best_value = 0;
    std::vector<size_t> best_idx{};
    std::vector<TSolution> parents{};
};

struct LogEntry {
    size_t generation;
    double value;
};

template <typename TSolution>
class Context {
private:
    Population<TSolution> population_storage_a{};
    Population<TSolution> population_storage_b{};

public:
    Settings settings;
    Random random{};
    Progress progress;
    Population<TSolution>& prev_generation;
    Population<TSolution>& next_generation;

    // important! generations begin with 1st, not 0th
    size_t generation = 1;
    bool stop_state_machine = false;
    std::chrono::time_point<std::chrono::steady_clock> start_time;
    std::chrono::time_point<std::chrono::steady_clock> stop_time;
    double cycle_time_us = 0;

    std::deque<LogEntry> log;
    std::ofstream log_f;

    Context(Settings& s) :
        settings(s),
        progress(s.n_generations, s.progress_update_period),
        prev_generation(population_storage_a),
        next_generation(population_storage_b)
    {
    }

    void swap_generations()
    {
        if (&population_storage_a == &prev_generation) {
            prev_generation = population_storage_b;
            next_generation = population_storage_a;
        }
        else {
            prev_generation = population_storage_a;
            next_generation = population_storage_b;
        }

        // clean the space for the next generation
        next_generation.solutions.resize(0);
        next_generation.parents.resize(0);
        next_generation.best_idx.resize(0);
    }

    const std::string get_best_solution_csv_filename()
    {
        std::stringstream ss;
        ss << settings.best_solution_filename_prefix << generation << ".csv";
        return ss.str();
    }

    std::string get_stats()
    {
        const double runtime_s =
                std::chrono::duration_cast<std::chrono::seconds>(stop_time -
                                                                 start_time)
                        .count();
        const double solution_s =
                (settings.population_size * settings.n_generations) /
                runtime_s;
        const double best_value = next_generation.best_value;
        const size_t first_col_width = 22;
        std::stringstream ss{};
        // standard parameters
        ss << std::left << std::setw(first_col_width) << "generations"
           << settings.n_generations << std::endl;
        ss << std::left << std::setw(first_col_width) << "population"
           << settings.population_size << std::endl;
        ss << std::left << std::setw(first_col_width) << "elite best"
           << settings.n_elite_best << std::endl;
        ss << std::left << std::setw(first_col_width) << "elite worst"
           << settings.n_elite_worst << std::endl;
        ss << std::left << std::setw(first_col_width) << "elite random"
           << settings.n_elite_random << std::endl;
        ss << std::left << std::setw(first_col_width) << "parents best"
           << settings.n_parents_best << std::endl;
        ss << std::left << std::setw(first_col_width) << "parents worst"
           << settings.n_parents_worst << std::endl;
        ss << std::left << std::setw(first_col_width) << "parents random"
           << settings.n_parents_random << std::endl;
        ss << std::left << std::setw(first_col_width) << "parents randomized"
           << settings.n_parents_randomized << std::endl;
        ss << std::left << std::setw(first_col_width) << "p replace solution"
           << settings.p_replace_solution << std::endl;
        ss << std::left << std::setw(first_col_width) << "p replace data"
           << settings.p_replace_data << std::endl;
        ss << std::left << std::setw(first_col_width) << "p swap data"
           << settings.p_swap_data << std::endl;
        ss << std::left << std::setw(first_col_width) << "p mutate data"
           << settings.p_mutate_data << std::endl;
        ss << std::left << std::setw(first_col_width) << "p mutate bad data"
           << settings.p_mutate_bad_data << std::endl;
        ss << std::endl;
        // custom parameters
        for (const auto& pair : settings.custom_parameter) {
            ss << std::left << std::setw(first_col_width) << pair.first
               << pair.second << std::endl;
        }
        ss << std::endl;
        // runtime stats
        ss << std::left << std::setw(first_col_width) << "runtime"
           << seconds_to_hhmmss_string(runtime_s) << std::endl;
        ss << std::left << std::setw(first_col_width) << "solutions/s"
           << solution_s << std::endl;
        ss << std::left << std::setw(first_col_width) << "best value"
           << best_value << std::endl;
        return ss.str();
    }
};

template <typename TSolution>
class StateMachine {
private:
    typedef std::function<void(Context<TSolution>&)> state_function_t;
    Context<TSolution> context;

public:
    std::vector<state_function_t> init_functions{};
    std::vector<state_function_t> state_functions{};
    std::vector<state_function_t> closure_functions{};

    StateMachine(Settings& s) :
        context(s)
    {
    }

    void run()
    {
        context.start_time = std::chrono::steady_clock::now();
        for (state_function_t& f : init_functions) {
            f(context);
        }
        auto cycle_begin_time = std::chrono::steady_clock::now();
        while (!context.stop_state_machine) {
            for (state_function_t& f : state_functions) {
                if (context.stop_state_machine) {
                    context.stop_time = std::chrono::steady_clock::now();
                    break;
                }
                f(context);
            }
            auto cycle_end_time = std::chrono::steady_clock::now();
            context.cycle_time_us =
                    std::chrono::duration_cast<std::chrono::microseconds>(
                            cycle_end_time - cycle_begin_time)
                            .count();
            cycle_begin_time = cycle_end_time;
        }
        context.progress.os_clean();
        for (state_function_t& f : closure_functions) {
            f(context);
        }
    }
};

// state machine states
// functions ot type StateMachinte::state_function_t

template <typename TSolution>
void change_generations(Context<TSolution>& c)
{
    if (c.generation >= c.settings.n_generations) {
        c.stop_state_machine = true;
        return;
    }
    c.generation++;
    c.swap_generations();
}

template <typename TSolution>
void print_context(Context<TSolution>& c)
{
    std::cout << "geneation: " << c.generation << std::endl;
    std::cout << "next_generation: " << std::endl;
    for (size_t i = 0; i < c.next_generation.solutions.size(); i++) {
        std::cout << "\t[" << i << "]:" << std::endl
                  << "\t\tsolutions:" << std::endl
                  << c.next_generation.solutions[i].str(3) << "\t\tvalue: "
                  << c.next_generation.solutions[i].get_value(c.settings)
                  << std::endl;
    }
    std::cout << "\tbest_idx:" << std::endl;
    for (size_t i = 0; i < c.next_generation.best_idx.size(); i++) {
        std::cout << "\t\t[" << i << "]:" << c.next_generation.best_idx[i]
                  << std::endl;
    }
    std::cout << "next_generation: " << std::endl;
    for (size_t i = 0; i < c.next_generation.solutions.size(); i++) {
        std::cout << "\t[" << i << "]:" << std::endl
                  << "\t\tsolutions:" << std::endl
                  << c.next_generation.solutions[i].str(3);
    }
    std::cout << "\tparents:" << std::endl;
    for (size_t i = 0; i < c.next_generation.parents.size(); i++) {
        std::cout << c.next_generation.parents[i].str(2);
    }
}

template <typename TSolution>
void print_value(Context<TSolution>& c)
{
    std::cout << "geneation: " << c.generation << std::endl;
    std::cout << "\tnext_generation:" << std::endl;
    std::cout << "\t\tvalue:" << std::endl;
    for (size_t i = 0; i < c.next_generation.solutions.size(); i++) {
        std::cout << "\t\t\t[" << i << "]:"
                  << c.next_generation.solutions[i].get_value(c.settings)
                  << std::endl;
    }
    std::cout << "\t\tbest_value: " << c.next_generation.best_value
              << std::endl;
}

template <typename TSolution>
void print_progress(Context<TSolution>& c)
{
    const double sol_s =
            c.settings.population_size / c.cycle_time_us * 1000000;

    std::stringstream ss;
    ss << " best " << c.next_generation.best_value;
    ss << " sol/s " << sol_s;

    c.progress.update(std::string(ss.str()));
}

template <typename TSolution>
void print_result(Context<TSolution>& c)
{
    std::stringstream ss;
    size_t best_idx = c.next_generation.best_idx[0];
    std::cout << "best result:" << std::endl;
    std::cout << c.next_generation.solutions[best_idx].str(1) << std::endl;
}

template <typename TSolution>
void print_stats(Context<TSolution>& c)
{
    std::cout << c.get_stats();
}

template <typename TSolution>
void create_stats_file(Context<TSolution>& c)
{
    if (c.settings.stats_filename.empty()) {
        return;
    }

    std::ofstream f(c.settings.stats_filename);
    f.is_open();
    f << c.get_stats();
}

template <typename TSolution>
void init_log(Context<TSolution>& c)
{
    assert(!c.log_f.is_open());
    if (c.settings.log_filename.empty()) {
        return;
    }

    c.log_f.open(c.settings.log_filename);
    c.log_f.is_open();
    c.log_f << "generation,best_value" << std::endl;
}

template <typename TSolution>
void update_log(Context<TSolution>& c)
{
    // update the log in program memory
    const LogEntry new_log_entry{c.generation, c.next_generation.best_value};
    c.log.push_back(new_log_entry);

    // write the log to file
    assert(!c.settings.log_filename.empty());
    if (!c.log_f.is_open()) {
        return;
    }
    while (!c.log.empty()) {
        const LogEntry val = c.log.front();
        c.log.pop_front();

        c.log_f << val.generation << "," << val.value << std::endl;
    }
    c.log_f << std::flush;
}

template <typename TSolution>
void create_best_solution_csv(Context<TSolution>& c)
{
    if (c.settings.best_solution_csv_creation_period != 0 &&
        (c.generation % c.settings.best_solution_csv_creation_period) != 0) {
        return;
    }

    TSolution& best_solution =
            c.next_generation.solutions[c.next_generation.best_idx[0]];
    best_solution.create_csv(c.get_best_solution_csv_filename());
}

template <typename TSolution>
void randomize_next_gen(Context<TSolution>& c)
{
    assert(c.generation == 1);
    assert(c.next_generation.solutions.size() == 0);

    for (size_t i = 0; i < c.settings.population_size; i++) {
        TSolution solution;
        solution.randomize(c.settings, [&c]() { return c.random.rnd01(); });
        c.next_generation.solutions.push_back(solution);
    }
}

template <typename TSolution>
void sort_next_gen_by_value(Context<TSolution>& c)
{
    // init .best_idx
    assert(c.next_generation.solutions.size() == c.settings.population_size);
    assert(c.next_generation.best_idx.size() == 0);
    for (size_t i = 0; i < c.settings.population_size; i++) {
        c.next_generation.best_idx.push_back(i);
    }
    assert(c.next_generation.best_idx.size() == c.settings.population_size);
    // fill .best_idx
    std::sort(c.next_generation.best_idx.begin(),
              c.next_generation.best_idx.end(),
              [&c](size_t a, size_t b) -> bool {
                  return c.next_generation.solutions[a].get_value(c.settings) >
                         c.next_generation.solutions[b].get_value(c.settings);
              });

    // update best value
    const size_t best_i = c.next_generation.best_idx[0];
    c.next_generation.best_value =
            c.next_generation.solutions[best_i].get_value(c.settings);
}

template <typename TSolution>
void select_next_gen_parents(Context<TSolution>& c)
{
    assert(c.prev_generation.best_idx.size() == c.settings.population_size);
    assert(c.next_generation.parents.size() == 0);

    // best
    assert(c.prev_generation.best_idx.size() >= c.settings.n_parents_best);
    for (size_t i = 0; i < c.settings.n_parents_best; i++) {
        const size_t parent_i = c.prev_generation.best_idx[i];
        c.next_generation.parents.push_back(
                c.prev_generation.solutions[parent_i]);
    }

    // worst
    assert(c.prev_generation.best_idx.size() >= c.settings.n_parents_worst);
    for (size_t i = 0; i < c.settings.n_parents_worst; i++) {
        const size_t parent_i =
                c.prev_generation
                        .best_idx[c.prev_generation.best_idx.size() - i - 1];
        c.next_generation.parents.push_back(
                c.prev_generation.solutions[parent_i]);
    }

    // random
    assert(c.prev_generation.best_idx.size() >= c.settings.n_parents_random);
    for (size_t i = 0; i < c.settings.n_parents_random; i++) {
        const size_t parent_i =
                c.random.rnd01() * c.prev_generation.best_idx.size();
        c.next_generation.parents.push_back(
                c.prev_generation.solutions[parent_i]);
    }

    // randomized
    for (size_t i = 0; i < c.settings.n_parents_randomized; i++) {
        TSolution solution;
        solution.randomize(c.settings, [&c]() { return c.random.rnd01(); });
        c.next_generation.parents.push_back(solution);
    }

    assert(c.next_generation.parents.size() ==
           (c.settings.n_parents_best + c.settings.n_parents_worst +
            c.settings.n_parents_random + c.settings.n_parents_randomized));
}

template <typename TSolution>
void add_next_gen_solutions_from_elite(Context<TSolution>& c)
{
#ifndef NDEBUG
    const size_t n_elite_total = c.settings.n_elite_best +
                                 c.settings.n_elite_worst +
                                 c.settings.n_elite_random;
    assert(c.next_generation.solutions.size() + n_elite_total <=
           c.settings.population_size);
#endif

    // best
    assert(c.prev_generation.best_idx.size() >= c.settings.n_elite_best);
    for (size_t i = 0; i < c.settings.n_elite_best; i++) {
        const size_t elite_i = c.prev_generation.best_idx[i];
        c.next_generation.solutions.push_back(
                c.prev_generation.solutions[elite_i]);
    }

    // worst
    assert(c.prev_generation.best_idx.size() >= c.settings.n_elite_worst);
    for (size_t i = 0; i < c.settings.n_elite_worst; i++) {
        const size_t elite_i =
                c.prev_generation
                        .best_idx[c.prev_generation.best_idx.size() - i - 1];
        c.next_generation.solutions.push_back(
                c.prev_generation.solutions[elite_i]);
    }

    // random
    assert(c.prev_generation.best_idx.size() >= c.settings.n_elite_random);
    for (size_t i = 0; i < c.settings.n_elite_random; i++) {
        const size_t elite_i =
                c.random.rnd01() * c.prev_generation.best_idx.size();
        c.next_generation.solutions.push_back(
                c.prev_generation.solutions[elite_i]);
    }

    assert(c.next_generation.solutions.size() <= c.settings.population_size);
}

template <typename TSolution>
void add_next_gen_solutions_from_crossover(Context<TSolution>& c)
{
    assert(c.next_generation.parents.size() >= 2);
    assert(c.next_generation.parents.size() ==
           (c.settings.n_parents_best + c.settings.n_parents_worst +
            c.settings.n_parents_random + c.settings.n_parents_randomized));

    const size_t n_elite_total = c.settings.n_elite_best +
                                 c.settings.n_elite_worst +
                                 c.settings.n_elite_random;
    assert(n_elite_total <= c.settings.population_size);
    const size_t n_children_total = c.settings.population_size - n_elite_total;
    assert(c.next_generation.solutions.size() + n_children_total <=
           c.settings.population_size);

    for (size_t i = 0; i < n_children_total; i++) {
        size_t parent_a_i =
                c.random.rnd01() * c.next_generation.parents.size();
        size_t parent_b_i;
        do {
            parent_b_i = c.random.rnd01() * c.next_generation.parents.size();
        } while (parent_b_i == parent_a_i);
        TSolution& parent_a = c.next_generation.parents[parent_a_i];
        TSolution& parent_b = c.next_generation.parents[parent_b_i];
        TSolution child{};
        child.crossover(
                c.settings, [&c]() { return c.random.rnd01(); }, parent_a,
                parent_b);
        c.next_generation.solutions.push_back(child);
    }

    assert(c.next_generation.solutions.size() <= c.settings.population_size);
}

template <typename TSolution>
void next_gen_replacements(Context<TSolution>& c)
{
    for (auto& s : c.next_generation.solutions) {
        s.replace(
                c.settings, [&c]() { return c.random.rnd01(); },
                c.next_generation.solutions);
    }
}

template <typename TSolution>
void next_gen_mutations(Context<TSolution>& c)
{
    for (auto& s : c.next_generation.solutions) {
        s.mutate(c.settings, [&c]() { return c.random.rnd01(); });
    }
}

template <typename TSolution>
void next_gen_swaps(Context<TSolution>& c)
{
    for (auto& s : c.next_generation.solutions) {
        s.swap(
                c.settings, [&c]() { return c.random.rnd01(); },
                c.next_generation.solutions);
    }
}

}  // namespace varga
