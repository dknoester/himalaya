/* delay_bench.cpp
 *
 * This file is part of the Himalaya project.
 *
 * Copyright 2012 David B. Knoester.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <ea/evolutionary_algorithm.h>
#include <ea/representations/bitstring.h>
#include <ea/fitness_functions/nk_model.h>
#include <ea/generational_models/steady_state.h>
#include <ea/generational_models/moran_process.h>
#include <ea/selection/tournament.h>
#include <ea/selection/proportionate.h>
#include <ea/selection/elitism.h>
#include <ea/selection/random.h>
#include <ea/line_of_descent.h>
#include <ea/datafiles/evaluations.h>
#include <ea/datafiles/fitness.h>
#include <ea/cmdline_interface.h>
#include <ea/fitness_functions/all_ones.h>
using namespace ealib;

#include "delay.h"
#include "analysis.h"

typedef evolutionary_algorithm
< individual<bitstring, peak_delay<nk_model< > >, bitstring, directS, default_lod_traits>
, ancestors::random_bitstring
, mutation::operators::per_site<mutation::site::bitflip>
, recombination::two_point_crossover
, generational_models::steady_state<selection::proportionate< >, delayed_elitism<selection::random> >
> ea_type;

/*! Define the EA's command-line interface.  Ealib provides an integrated command-line
 and configuration file parser.  This class specializes that parser for this EA.
 */
template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    //! Define the options that can be parsed.
    virtual void gather_options() {
        add_option<REPRESENTATION_SIZE>(this);
        add_option<POPULATION_SIZE>(this);
        add_option<STEADY_STATE_LAMBDA>(this);
        add_option<TOURNAMENT_SELECTION_N>(this);
        add_option<TOURNAMENT_SELECTION_K>(this);
        add_option<ELITISM_N>(this);
        
        add_option<MUTATION_PER_SITE_P>(this);
        
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<CHECKPOINT_OFF>(this);
        add_option<CHECKPOINT_PREFIX>(this);
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        add_option<ANALYSIS_OUTPUT>(this);
        
        add_option<FF_RNG_SEED>(this);
        add_option<NK_MODEL_N>(this);
        add_option<NK_MODEL_K>(this);
        
        add_option<DELAY_GENERATIONS>(this);
        add_option<DELAY_RANDOM_INSERT>(this);
    }
    
    //! Define events (e.g., datafiles) here.
    virtual void gather_events(EA& ea) {
        add_event<datafiles::fitness>(ea);
        add_event<datafiles::fitness_evaluations>(ea);
        add_event<lod_event>(ea);
        add_event<effective_fitness>(ea);
        add_event<dominant_archive>(ea);
        add_event<random_individuals>(ea);
    };
    
    virtual void gather_tools() {
    }
};

LIBEA_CMDLINE_INSTANCE(ea_type, cli);