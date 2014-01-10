/* himalaya.h
 *
 * This file is part of EALib Examples.
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
#ifndef _HIMALAYA_H_
#define _HIMALAYA_H_

#include <boost/lexical_cast.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>

#include <ea/evolutionary_algorithm.h>
#include <ea/representations/realstring.h>
#include <ea/representations/bitstring.h>
#include <ea/cmdline_interface.h>
#include <ea/line_of_descent.h>
#include <ea/selection/proportionate.h>
#include <ea/selection/tournament.h>
#include <ea/selection/elitism.h>
#include <ea/selection/random.h>
#include <ea/generational_models/steady_state.h>
#include <ea/datafiles/evaluations.h>
#include <ea/datafiles/fitness.h>
#include <ea/datafiles/population_entropy.h>
#include <ea/fitness_functions/benchmarks.h>
#include <ea/fitness_functions/nk_model.h>
#include <ea/events.h>
#include <ea/algorithm.h>
#include <ea/analysis.h>
#include <ea/datafile.h>
#include <ea/attributes.h>

#include "analysis.h"
using namespace ealib;


LIBEA_MD_DECL(HIMALAYA_MUTATION_RATE, "himalaya.mutation_rate", double);
LIBEA_MD_DECL(HIMALAYA_ACC_OFFSPRING_DELTAW, "himalaya.acc_offspring_deltaw", double);
LIBEA_MD_DECL(HIMALAYA_NUM_OFFSPRING, "himalaya.num_offspring", double);


/*! Per-site mutation.
 */
template <typename MutationType>
struct himalaya_per_site {
    typedef MutationType mutation_type;
    
    //! Iterate through all elements in the given representation, possibly mutating them.
    template <typename EA>
    void operator()(typename EA::individual_type& ind, EA& ea) {
        typename EA::representation_type& repr=ind.repr();
        const double per_site_p=get<MUTATION_PER_SITE_P>(ea);
        for(typename EA::representation_type::iterator i=repr.begin(); i!=repr.end(); ++i){
            if(ea.rng().p(per_site_p)) {
                _mt(repr, i, ea);
            }
        }
    }
    
    mutation_type _mt;
};

template <typename EA>
struct himalaya_inheritance : inheritance_event<EA> {
    himalaya_inheritance(EA& ea) : inheritance_event<EA>(ea) {
    }
    virtual ~himalaya_inheritance() { }
    
    virtual void operator()(typename EA::population_type& parents,
                            typename EA::individual_type& offspring,
                            EA& ea) {
        // calculate w_off:
        double w_off = static_cast<double>(ealib::fitness(offspring,ea));
        
        // accumulate delta_w in parent:
        for(typename EA::population_type::iterator i=parents.begin(); i!=parents.end(); ++i) {
            if((*i)->generation() >= 0.0) {
                double delta_w = w_off - static_cast<double>(ealib::fitness(**i,ea));
                put<HIMALAYA_ACC_OFFSPRING_DELTAW>(get<HIMALAYA_ACC_OFFSPRING_DELTAW>(**i,0.0) + delta_w, **i);
                put<HIMALAYA_NUM_OFFSPRING>(get<HIMALAYA_NUM_OFFSPRING>(**i,0.0)+1.0, **i);
            }
        }
    }
};


template <typename EA>
struct himalaya : public ealib::end_of_update_event<EA> {
    himalaya(EA& ea) : ealib::end_of_update_event<EA>(ea) {
    }
    
    virtual ~himalaya() { }
    
    virtual void operator()(EA& ea) {
        using namespace ealib;

        // fitness peaks are defined relative to a genotype; the only way to identify
        // if an individual lies on a slope vs on the peak is to map it's local
        // neighborhood, which means asexual search / hill climbing.
        
        // recombination, however, serves to teleport the offspring to an entirely
        // different peak.  the genetic relationship between the offspring and
        // its parents doesn't matter; it's on a different peak (unless recombination
        // was degenerate, and only mutated a single locus).
        
        // either way, recombination doesn't inform us of whether an individual
        // is on a fruitful peak or not.  the goal, however, is to identify the
        // set of peaks.  once they've been identified, the trick is to climb them.
        // note that there are many paths up each peak -- some may be easier than
        // others.
        
        // track the delta_w between parent and offspring:
        //  if >0, parent is in valley or on slope
        //  if 0, parent is on a neutral plateau
        //  if <0, parent is on/near peak
        
        // want to:
        // - perform a broad search from a valley or slope (many paths)
        // - explore neutral plateaus (many paths)
        // - search nearby for absolute peak (fine tune)
        
        // can we detect individuals that are good parents?  i.e., those that
        // provide useful building blocks?
        
        // could replace individuals wrt to the entropy that they contribute to
        // the population...
        
        // measure genotypic entropy:
        std::vector<std::string> genotypes;
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            genotypes.push_back(algorithm::vcat(i->repr().begin(), i->repr().end()));
        }
        
        double mu = get<MUTATION_PER_SITE_P>(ea);
        
        // for now, scale mu by a fixed 1% based on genotypic entropy.
        // this is probably a **really** bad set point...
        if(analysis::entropy(genotypes.begin(), genotypes.end()) > 0.8) {
            mu *= 0.99;
        } else {
            mu *= 1.01;
        }
        
        // make sure it stays within a somewhat reasonable range:
        mu = algorithm::clip(mu,
                             1.0/(2.0*get<REPRESENTATION_SIZE>(ea)), // floor mu at ~0.5 genomic
                             0.5); // ceiling mu at 0.5 per site
        
        put<MUTATION_PER_SITE_P>(mu,ea);
        
        //        // adapt the heritable mutation rate in response to an individual's lineage's
        //        // ema of fitness; if it's increasing, don't touch it.  if it's plateaued,
        //        // increase it.  perhaps use semideviation?
        //        for(typename EA::population_type::iterator i=ea.population().begin(); i!=ea.population().end(); ++i) {
        //            int c=0;
        //            typename EA::individual_ptr_type p=*i;
        //            std::vector<double> f;
        //            while(p->has_parents() && c++<20) {
        //                f.push_back(ealib::fitness(**i,ea) - ealib::fitness(*p,ea));
        //                p = p->lod_parent();
        //            }
        //            if(algorithm::exp_mean(f.begin(), f.end(), 5ul) < 5.0) {
        //                put<HIMALAYA_MUTATION_RATE_OFF>(get<HIMALAYA_MUTATION_RATE>(**i)*1.001,**i);
        //            } else {
        //                put<HIMALAYA_MUTATION_RATE_OFF>(1.0,**i);
        //            }
        //        }
    }
};


template <typename EA>
struct himalaya_datafile : ealib::record_statistics_event<EA> {
    himalaya_datafile(EA& ea) : ealib::record_statistics_event<EA>(ea), _df("himalaya.dat") {
        _df.add_field("update")
        .add_field("mu")
        .add_field("avg_delta_w")
        .add_field("avg_num_offspring");
    }
    
    virtual ~himalaya_datafile() {
    }
    
    virtual void operator()(EA& ea) {
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::mean> > delta_w, numoffspring;
        
        for(typename EA::population_type::iterator i=ea.population().begin(); i!=ea.population().end(); ++i) {
            delta_w(get<HIMALAYA_ACC_OFFSPRING_DELTAW>(**i,0.0) / get<HIMALAYA_NUM_OFFSPRING>(**i,1.0));
            numoffspring(get<HIMALAYA_NUM_OFFSPRING>(**i,0.0));
        }
        
        _df.write(ea.current_update())
        .write(get<MUTATION_PER_SITE_P>(ea))
        .write(mean(delta_w))
        .write(mean(numoffspring))
        .endl();
    }
    
    ealib::datafile _df;
};

#endif
