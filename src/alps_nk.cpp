#include <sstream>
#include <boost/lexical_cast.hpp>
#include <ea/evolutionary_algorithm.h>
#include <ea/fitness_function.h>
#include <ea/selection/tournament.h>
#include <ea/selection/elitism.h>
#include <ea/selection/random.h>
#include <ea/generational_models/steady_state.h>
#include <ea/datafiles/evaluations.h>
#include <ea/datafiles/fitness.h>
#include <ea/datafiles/population_entropy.h>
#include <ea/representations/bitstring.h>
#include <ea/cmdline_interface.h>
#include <ea/algorithm.h>
#include <ea/alps.h>
#include <ea/meta_population.h>
#include <ea/fitness_functions/nk_model.h>

using namespace ealib;

template <typename EA>
struct configuration : public abstract_configuration<EA> {
    //! Called to generate the initial EA population.
    void initial_population(EA& ea) {
        generate_ancestors(ancestors::random_bitstring(), get<POPULATION_SIZE>(ea), ea);
    }
};


//! Evolutionary algorithm definition.
typedef meta_population<
evolutionary_algorithm<
bitstring,
mutation::operators::per_site<mutation::site::bitflip>, // mutation operator
nk_model< >,
configuration,
recombination::two_point_crossover,
generational_models::steady_state<selection::tournament< >, selection::elitism<selection::random> >
> mea_type;


/*! Define the EA's command-line interface.
 */
template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    virtual void gather_options() {
        // ea options
        add_option<META_POPULATION_SIZE>(this);
        add_option<REPRESENTATION_SIZE>(this);
        
        add_option<GM_AGE>(this);
        add_option<ADMISSION_AGE>(this);
        add_option<GAP_SIZE>(this);
        add_option<ADMISSION_AGING_SCHEME>(this);

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
        
        add_option<NK_MODEL_N>(this);
        add_option<NK_MODEL_K>(this);
    }
    
    virtual void gather_tools() {
    }
    
    virtual void gather_events(EA& ea) {
        add_event<alps>(this, ea);
        add_event<alps_datafile>(this, ea);
        add_event<datafiles::meta_population_entropy>(this, ea);
        add_event<datafiles::meta_population_fitness>(this, ea);
        add_event<datafiles::meta_population_fitness_evaluations>(this, ea);
    };
};
LIBEA_CMDLINE_INSTANCE(mea_type, cli);
