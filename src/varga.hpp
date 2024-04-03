#include <vector>

template <typename TIndividualValue>
struct Individual
{
    TIndividualValue value;
    double fitness;
};

template <typename TIndividualValue>
struct Population
{
    std::vector<Individual<TIndividualValue> individuals;
};

template <typename TIndividualValue>
struct Context
{
    size_t ngeneration;
    Population<TIndividualValue> new_generation;
    Population<TIndividualValue> prev_generation;
};

template <typename TIndividualValue>
struct Initializer
{
    // TODO: define input/output parameter and type
    // TODO: how exactly it should work
    // TODO: does it need to be a class, or better - a function?
    // TODO: can this become a Context constructor?
    virtual void init(Context<TIndividualValue> &context) {};
};

template <typename TIndividualValue>
struct Evaluator
{
    // TODO: define input/output parameter and type
    // TODO: how exactly it should work
    virtual void evaluate(Context<TIndividualValue> &context) {};
};

template <typename TIndividualValue>
struct Generator
{
    // TODO: define input/output parameter and type
    // TODO: how exactly it should work
    virtual bool generate(Context<TIndividualValue> &context) {};
};

template <typename TIndividualValue>
struct MainLoop {
    Context<TIndividualValue> context;
    Initializer<TIndividualValue> initializer;
    Evaluator<TIndividualValue> evaluator;
    Generator<TIndividualValue> generator;

    // TODO: add constructor
    void run() 
    {
        // TODO: add proper run() logic and parameter passing
        initializer.init(context);
        bool continue_generating = true;
        while (continue_generating) {
            evaluator.evaluate(context);
            continue_generating = generator.generate(context);
        }
    }
};
