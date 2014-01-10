#include <sstream>
#include <boost/lexical_cast.hpp>
#include <ea/evolutionary_algorithm.h>
#include <ea/fitness_function.h>
#include <ea/selection/proportionate.h>
#include <ea/selection/tournament.h>
#include <ea/selection/elitism.h>
#include <ea/selection/random.h>
#include <ea/generational_models/steady_state.h>
#include <ea/datafiles/evaluations.h>
#include <ea/datafiles/fitness.h>
#include <ea/datafiles/population_entropy.h>
#include <ea/representations/realstring.h>
#include <ea/cmdline_interface.h>
#include <ea/algorithm.h>
#include <ea/meta_population.h>
#include <ea/fitness_functions/benchmarks.h>
#include <ea/qhfc.h>


using namespace ealib;

template <typename EA>
struct configuration : public abstract_configuration<EA> {
    //! Called to generate the initial EA population.
    void initial_population(EA& ea) {
        generate_ancestors(ancestors::uniform_real(), get<POPULATION_SIZE>(ea), ea);
    }
    
    //! Called to fill a population to capacity.
    void fill_population(EA& ea) {
        generate_ancestors(ancestors::uniform_real(), get<POPULATION_SIZE>(ea)-ea.size(), ea);
    }
};


//! Evolutionary algorithm definition.
typedef evolutionary_algorithm<
realstring,
mutation::operators::per_site<mutation::site::uniform_real>, // mutation operator
benchmarks,
configuration,
recombination::two_point_crossover, // recombination operator
generational_models::deterministic_crowding< >
> ea_type;

typedef meta_population<
// embedded ea type:
evolutionary_algorithm<bitstring,
mutation::operators::per_site<mutation::site::bitflip>,
all_ones,
configuration,
recombination::two_point_crossover,
generational_models::deterministic_crowding< > >,
// mp types:
mutation::operators::no_mutation,
constant,
qhfc_configuration,
recombination::no_recombination,
generational_models::qhfc,
attr::no_attributes> ea_type;



template <typename EA>
struct mp_configuration : public abstract_configuration<EA> {
    //! Called as the final step of EA initialization.
    virtual void initialize(EA& ea) {
        for(std::size_t i=0; i<get<META_POPULATION_SIZE>(ea); ++i) {
        }
    }
};

typedef meta_population<ea_type, mp_configuration, generational_models::qhfc> mea_type;


/*! Define the EA's command-line interface.
 */
template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    virtual void gather_options() {
        // ea options
        add_option<META_POPULATION_SIZE>(this);
        add_option<REPRESENTATION_SIZE>(this);
        
        add_option<POPULATION_SIZE>(this);
        add_option<REPLACEMENT_RATE_P>(this);
        add_option<TOURNAMENT_SELECTION_N>(this);
        add_option<TOURNAMENT_SELECTION_K>(this);
        add_option<ELITISM_N>(this);
        
        add_option<MUTATION_PER_SITE_P>(this);
        add_option<MUTATION_UNIFORM_REAL_MIN>(this);
        add_option<MUTATION_UNIFORM_REAL_MAX>(this);
        
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<CHECKPOINT_OFF>(this);
        add_option<CHECKPOINT_PREFIX>(this);
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        add_option<ANALYSIS_OUTPUT>(this);
        
        add_option<BENCHMARKS_FUNCTION>(this);
        add_option<QHFC_BREED_TOP_FREQ>(this);
        add_option<QHFC_DETECT_EXPORT_NUM>(this);
        add_option<QHFC_PERCENT_REFILL>(this);
        add_option<QHFC_CATCHUP_GEN>(this);
        add_option<QHFC_NO_PROGRESS_GEN>(this);
    }
    
    virtual void gather_tools() {
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::qhfc>(this, ea);
        add_event<datafiles::meta_population_entropy>(this, ea);
        add_event<datafiles::meta_population_fitness>(this, ea);
        add_event<datafiles::meta_population_fitness_evaluations>(this, ea);
    };
};
LIBEA_CMDLINE_INSTANCE(mea_type, cli);
