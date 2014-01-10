/* himalaya.cpp
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
#include "himalaya2.h"

template <typename EA>
struct configuration : public abstract_configuration<EA> {
    //! Called to generate the initial EA population.
    void initial_population(EA& ea) {
        generate_ancestors(ancestors::uniform_real(), get<POPULATION_SIZE>(ea), ea);
    }
};

typedef evolutionary_algorithm<
realstring,
mutation::operators::per_site<mutation::site::uniform_real>, // mutation operator
benchmarks,
configuration, // user-defined configuration methods
recombination::two_point_crossover, // recombination operator
generational_models::steady_state<selection::proportionate< >, selection::elitism<selection::random> >, // generational model
attr::default_attributes, // individual attributes
individual_lod // using an lod individual automatically turns on LOD tracking.
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

        add_option<REPLACEMENT_RATE_P>(this);
        add_option<MUTATION_GENOMIC_P>(this);
        add_option<MUTATION_PER_SITE_P>(this);
        add_option<MUTATION_UNIFORM_REAL_MIN>(this);
        add_option<MUTATION_UNIFORM_REAL_MAX>(this);

        add_option<TOURNAMENT_SELECTION_N>(this);
        add_option<TOURNAMENT_SELECTION_K>(this);
        add_option<ELITISM_N>(this);

        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<CHECKPOINT_OFF>(this);
        add_option<CHECKPOINT_PREFIX>(this);
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        add_option<ANALYSIS_OUTPUT>(this);
        
        add_option<BENCHMARKS_FUNCTION>(this);
    }
    
    //! Define events (e.g., datafiles) here.
    virtual void gather_events(EA& ea) {
        add_event<himalaya>(this,ea);
        add_event<himalaya_datafile>(this,ea);
        add_event<himalaya_inheritance>(this,ea);
        add_event<datafiles::fitness>(this, ea);
        add_event<datafiles::fitness_evaluations>(this, ea);
        add_event<datafiles::population_entropy>(this, ea);
    };
    
    virtual void gather_tools() {
        add_tool<landscape_data>(this);
    }
};

LIBEA_CMDLINE_INSTANCE(ea_type, cli);