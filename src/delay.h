/* delay.h
 *
 * This file is part of the Himalaya project.
 *
 * Copyright 2014 David B. Knoester.
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
#ifndef _DELAY_H_
#define _DELAY_H_

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include <ea/metadata.h>
#include <ea/selection/elitism.h>

using namespace ealib;

LIBEA_MD_DECL(DELAY_GENERATIONS, "delay.generations", int);
LIBEA_MD_DECL(DELAY_W_REAL, "delay.w_real", double);
LIBEA_MD_DECL(DELAY_W_EFF, "delay.w_eff", double);
LIBEA_MD_DECL(DELAY_RANDOM_INSERT, "delay.random_insert", double);



namespace access {
    struct delayed_priority {
        template <typename EA>
        double operator()(typename EA::individual_type& ind, EA& ea) {
            typedef typename EA::individual_ptr_type individual_ptr_type;
            
            double w = ind.priority();
            put<DELAY_W_REAL>(w, ind);
            
            if(get<IND_GENERATION>(ind) > 0) {
                individual_ptr_type p=ind.traits().lod_parent();
                for(int i=1; (i<=get<DELAY_GENERATIONS>(ea)) && (get<IND_GENERATION>(*p) >= 0); ++i) {
                    w = get<DELAY_W_REAL>(*p);
                    p = p->traits().lod_parent();
                }
            }
            
            put<DELAY_W_EFF>(w, ind);
            return w;
        }
    };
    
}

template <typename EA>
typename EA::individual_ptr_type lod_parent(typename EA::individual_type& ind, EA& ea) {
    if((get<IND_GENERATION>(ind) <= 0) || (!ind.traits().has_parents())) {
        return typename EA::individual_ptr_type();
    }
    
    return ind.traits().lod_parent();
}

/*! Delay the fitness of an individual based on the mean fitness along its
 lineage.
 */
template <typename FitnessFunction>
struct mean_delay : public FitnessFunction {
    typedef FitnessFunction parent;

    //! Common mean delay method.
    template <typename Individual, typename EA>
    double delay(Individual& ind, EA& ea) {
        using namespace boost::accumulators;
        
        accumulator_set<double, stats<tag::mean> > w;
        w(get<DELAY_W_REAL>(ind));
        
        typename EA::individual_ptr_type p=lod_parent(ind,ea);
        for(int n=0; (n<get<DELAY_GENERATIONS>(ea)) && (p!=0); ++n) {
            w(get<DELAY_W_REAL>(*p));
            p = lod_parent(*p,ea);
        }
        
        double w1 = mean(w);
        put<DELAY_W_EFF>(w1,ind);
        return w1;

    }
    
    //! Mean delay a stochastic fitness function.
    template <typename Individual, typename RNG, typename EA>
    double operator()(Individual& ind, RNG& rng, EA& ea) {
        put<DELAY_W_REAL>(static_cast<double>(parent::operator()(ind,rng,ea)), ind);
        return delay(ind,ea);
    }

    //! Mean delay a constant fitness function.
    template <typename Individual, typename EA>
    double operator()(Individual& ind, EA& ea) {
        put<DELAY_W_REAL>(static_cast<double>(parent::operator()(ind,ea)), ind);
        return delay(ind,ea);
    }
};

/*! Delay the fitness of an indivdual by up to DELAY_GENERATIONS number of
 ancestors along its lineage.
 
 For example, the effective fitness w_eff of individual i is the real fitness
 w_real of its n'th ancestor.
 */
template <typename FitnessFunction>
struct generation_delay : public FitnessFunction {
    typedef FitnessFunction parent;
    
    template <typename Individual, typename EA>
    double operator()(Individual& ind, EA& ea) {
        typedef typename EA::individual_ptr_type individual_ptr_type;
        
        double w = parent::operator()(ind,ea);
        put<DELAY_W_REAL>(w, ind);
        
        if(get<IND_GENERATION>(ind) > 0) {
            individual_ptr_type p=ind.traits().lod_parent();
            for(int i=1; (i<=get<DELAY_GENERATIONS>(ea)) && (get<IND_GENERATION>(*p) >= 0); ++i) {
                w = get<DELAY_W_REAL>(*p);
                p = p->traits().lod_parent();
            }
        }
        
        put<DELAY_W_EFF>(w, ind);
        return w;
    }
};

/* Rewrite the fitness of an individual to be the max fitness of its DELAY_GENERATIONS
 ancestors.
 
 For example, w_eff of individual i is the max w_real of its n ancestors (itself
 included).
 */
template <typename FitnessFunction>
struct peak_delay : public FitnessFunction {
    typedef FitnessFunction parent;
    
    //! Compare two fitnesses based on the direction_tag:
    double best(double x, double y, minimizeS) { return std::min(x,y); }
    double best(double x, double y, maximizeS) { return std::max(x,y); }
    
    //! Calculate fitness.
    template <typename Individual, typename EA>
    double operator()(Individual& ind, EA& ea) {
        typedef typename EA::individual_ptr_type individual_ptr_type;
        
        double w = parent::operator()(ind,ea);
        put<DELAY_W_REAL>(w, ind);
        
        individual_ptr_type p=ind.traits().lod_parent();
        for(int i=1; (i<=get<DELAY_GENERATIONS>(ea)) && (get<IND_GENERATION>(*p) >= 0); ++i) {
            w = best(w,get<DELAY_W_REAL>(*p),typename parent::direction_tag());
            p = p->traits().lod_parent();
        }
        
        put<DELAY_W_EFF>(w, ind);
        return w;
    }
};



/*! Stacks elitism on top of another selection strategy.
 
 This selection strategy "stacks" with others; that is, it must be used in
 conjunction with another selection strategy, such as tournament_selection.
 Elitism augments that selection strategy by explicitly preserving N "elite"
 (high-fitness) individuals.  Those selected are *still* maintained as part
 of the source population from which the embedded selection strategy draws its
 own selected individuals.
 
 UPDATE: Keep the N highest-fit, most diverse individuals...
 
 */
template <typename SelectionStrategy>
struct delayed_elitism {
    typedef SelectionStrategy embedded_selection_type;
    
    //! Initializing constructor.
    template <typename Population, typename EA>
    delayed_elitism(std::size_t n, Population& src, EA& ea) : _embedded(n,src,ea) {
    }
    
    /*! Preserve the elite individuals from the src population.
     */
    template <typename Population, typename EA>
    void operator()(Population& src, Population& dst, std::size_t n, EA& ea) {
        std::size_t e = get<ELITISM_N>(ea);
        assert(n > e);
        _embedded(src, dst, n-e, ea);
        
        // now, append the e most-fit individuals:
        if(e > 0) {
            calculate_fitness(src.begin(), src.end(), ea);
            std::sort(src.begin(), src.end(), comparators::metadata<DELAY_W_REAL,EA>());
            typename Population::reverse_iterator rl=src.rbegin();
            std::advance(rl, e);
            dst.insert(dst.end(), src.rbegin(), rl);
        }
    };
    
    embedded_selection_type _embedded; //!< Underlying selection strategy.
};


/*! At the end of each update, insert random individuals into the population.
 */
template <typename EA>
struct random_individuals : end_of_update_event<EA> {
    random_individuals(EA& ea) : end_of_update_event<EA>(ea) {
    }
    
    virtual void operator()(EA& ea) {
        generate_ancestors(typename EA::ancestor_generator_type(), get<DELAY_RANDOM_INSERT>(ea)*get<POPULATION_SIZE>(ea), ea);
    }
};


/*! Store the dominant individual (based on real fitness).
 */
template <typename EA>
struct dominant_archive : fitness_evaluated_event<EA> {
    dominant_archive(EA& ea) : fitness_evaluated_event<EA>(ea), _df("dominant_archive.dat") {
        _df.add_field("update")
        .add_field("dominant_w_real");
    }

    virtual void operator()(typename EA::individual_type& ind, EA& ea) {
        if(_archive.empty()
           || (get<DELAY_W_REAL>(ind) > get<DELAY_W_REAL>(*_archive.back()))) {
            typename EA::individual_ptr_type p=ea.copy_individual(ind);
            p->traits().lod_clear();
            _archive.push_back(p);
            _df.write(ea.current_update()).write(get<DELAY_W_REAL>(*p)).endl();
        }
    }
    typename EA::population_type _archive;
    datafile _df;
};

/*! Datafile for mean generation, and mean & max fitness.
 */
template <typename EA>
struct effective_fitness : record_statistics_event<EA> {
    effective_fitness(EA& ea) : record_statistics_event<EA>(ea), _df("effective_fitness.dat") {
        _df.add_field("update")
        .add_field("mean_w_real")
        .add_field("max_w_real")
        .add_field("mean_w_eff");
    }
    
    virtual ~effective_fitness() {
    }
    
    virtual void operator()(EA& ea) {
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::mean,tag::max> > w_real;
        accumulator_set<double, stats<tag::mean> > w_eff;
        
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            w_real(get<DELAY_W_REAL>(*i));
            w_eff(get<DELAY_W_EFF>(*i));
        }
        
        _df.write(ea.current_update())
        .write(mean(w_real))
        .write(max(w_real))
        .write(mean(w_eff))
        .endl();
    }
    
    datafile _df;
};


#endif
