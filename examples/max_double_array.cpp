#include <sstream>
#include <string>

#include "varga.hpp"

const size_t g_n_generations = 1000000;
const size_t g_population_size = 10;
const size_t g_individual_n_genes = 100;
const size_t g_n_parents_best = 3;
const double g_p_mutate_gene = 0.05;

struct MyIndividual : varga::Individual<std::vector<double> > {
    MyIndividual()
    {
        genes.resize(g_individual_n_genes);
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
    }

    double get_fitness(varga::Settings& s)
    {
        (void)s;
        double fitness = 0;
        for (auto& g : genes) {
            fitness += g;
        }
        return fitness;
    }

    void crossover(varga::Settings& s,
                   const std::function<double(void)>& rnd01,
                   varga::Individual<std::vector<double> >& parent_a,
                   varga::Individual<std::vector<double> >& parent_b)
    {
        uniform_crossover(s, rnd01, parent_a, parent_b);
    }

    void mutate(varga::Settings& s, const std::function<double(void)>& rnd01)
    {
        assert(genes.size() != 0);
        if (rnd01() < s.p_mutate_gene) {
            size_t mutation_i = rnd01() * genes.size();
            genes[mutation_i] = rnd01();
        }
    }
};

int main()
{
    varga::Settings s{g_population_size, g_n_generations};
    s.n_parents_best = g_n_parents_best;
    s.progress_update_period = 10000;
    s.p_mutate_gene = g_p_mutate_gene;

    varga::StateMachine<MyIndividual> sm{s};
    sm.init_functions = {varga::init_log<MyIndividual>,
                         varga::randomize_next_gen<MyIndividual>};
    sm.state_functions = {
            varga::evaluate_next_gen<MyIndividual>,
            varga::sort_next_gen_by_fitness<MyIndividual>,
            varga::update_log<MyIndividual>,
            varga::print_progress<MyIndividual>,
            varga::change_generations<MyIndividual>,
            varga::select_next_gen_parents<MyIndividual>,
            varga::add_next_gen_individuals_from_crossover<MyIndividual>,
            varga::next_gen_mutations<MyIndividual>};
    sm.closure_functions = {varga::print_stats<MyIndividual>,
                            varga::create_stats_file<MyIndividual>,
                            varga::create_best_individual_csv<MyIndividual>};
    sm.run();
    return 0;
}
