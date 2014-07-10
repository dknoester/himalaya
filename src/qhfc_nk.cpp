#include <ea/qhfc.h>
#include <ea/genome_types/bitstring.h>
#include <ea/fitness_functions/nk_model.h>
#include <ea/datafiles/evaluations.h>
#include <ea/datafiles/metapopulation_fitness.h>
#include <ea/cmdline_interface.h>
using namespace ealib;

typedef qhfc
< direct<bitstring>
, nk_model< >
, mutation::operators::per_site<mutation::site::bitflip>
, recombination::two_point_crossover
, ancestors::random_bitstring
> ea_type;


/*! Define the EA's command-line interface.
 */
template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    virtual void gather_options() {
        // ea options
        add_option<METAPOPULATION_SIZE>(this);
        add_option<REPRESENTATION_SIZE>(this);
        
        add_option<POPULATION_SIZE>(this);
        
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
        
        add_option<ELITISM_N>(this);
        add_option<QHFC_POP_SCALE>(this);
        add_option<QHFC_BREED_TOP_FREQ>(this);
        add_option<QHFC_DETECT_EXPORT_NUM>(this);
        add_option<QHFC_PERCENT_REFILL>(this);
        add_option<QHFC_CATCHUP_GEN>(this);
        add_option<QHFC_NO_PROGRESS_GEN>(this);
        
        add_option<FF_RNG_SEED>(this);
        add_option<NK_MODEL_N>(this);
        add_option<NK_MODEL_K>(this);
    }
    
    virtual void gather_tools() {
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::qhfc_dat>(ea);
//        add_event<datafiles::meta_population_fitness>(ea);
//        add_event<datafiles::meta_population_fitness_evaluations>(ea);
    };
};
LIBEA_CMDLINE_INSTANCE(ea_type, cli);
