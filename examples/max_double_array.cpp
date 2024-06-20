#include <sstream>
#include <string>

#include "varga.hpp"

const size_t g_n_generations = 1000000;
const size_t g_population_size = 10;
const size_t g_solution_n_genes = 100;
const size_t g_n_parents_best = 3;
const double g_p_mutate_gene = 0.05;
const std::string g_log_filename = "max_double_array_log.csv";

struct MySolution : varga::Solution<std::vector<double> > {
    MySolution()
    {
        genes.resize(g_solution_n_genes);
    }

    std::string str(size_t n_tabs = 0)
    {
        std::stringstream ss;
        for (size_t n = 0; n < n_tabs; n++) {
            ss << "\t";
        }
        ss << "genes:" << std::endl;
        for (size_t i = 0; i < genes.size(); i++) {
            for (size_t n = 0; n < n_tabs; n++) {
                ss << "\t";
            }
            ss << "\t[" << i << "]:" << genes[i];
            if ((i + 1) != genes.size()) {
                ss << std::endl;
            }
        }
        return ss.str();
    }

    void create_csv(const std::string filename)
    {
        std::ofstream f(filename);
        f.is_open();
        f << "i,value" << std::endl;
        for (size_t i = 0; i < genes.size(); i++) {
            f << i << "," << genes[i] << std::endl;
        }
    }

    void randomize(varga::Settings& s,
                   const std::function<double(void)>& rnd01)
    {
        (void)s;
        for (auto& g : genes) {
            g = rnd01();
        }
        _changed = true;
    }

    double get_value(varga::Settings& s)
    {
        (void)s;
        if (_changed) {
            _value = 0;
            for (auto& g : genes) {
                _value += g;
            }
            _changed = false;
        }
        return _value;
    }

    void crossover(varga::Settings& s,
                   const std::function<double(void)>& rnd01,
                   varga::Solution<std::vector<double> >& parent_a,
                   varga::Solution<std::vector<double> >& parent_b)
    {
        uniform_crossover(s, rnd01, parent_a, parent_b);
    }

    void mutate(varga::Settings& s, const std::function<double(void)>& rnd01)
    {
        assert(genes.size() != 0);
        if (rnd01() < s.p_mutate_gene) {
            size_t mutation_i = rnd01() * genes.size();
            genes[mutation_i] = rnd01();
            _changed = true;
        }
    }
};

int main()
{
    varga::Settings s{g_population_size, g_n_generations};
    s.n_parents_best = g_n_parents_best;
    s.progress_update_period = 10000;
    s.p_mutate_gene = g_p_mutate_gene;
    s.log_filename = g_log_filename;

    varga::StateMachine<MySolution> sm{s};
    sm.init_functions = {varga::init_log<MySolution>,
                         varga::randomize_next_gen<MySolution>};
    sm.state_functions = {
            varga::sort_next_gen_by_value<MySolution>,
            varga::update_log<MySolution>,
            varga::print_progress<MySolution>,
            varga::change_generations<MySolution>,
            varga::select_next_gen_parents<MySolution>,
            varga::add_next_gen_solutions_from_crossover<MySolution>,
            varga::next_gen_mutations<MySolution>};
    sm.closure_functions = {varga::print_stats<MySolution>,
                            varga::create_stats_file<MySolution>,
                            varga::create_best_solution_csv<MySolution>};
    sm.run();
    return 0;
}
