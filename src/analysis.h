/* analysis.h
 *
 * This file is part of EALib.
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
#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_

///*! Sample the fitness landscape, and write fitness values to a file.
// */
//template <typename EA>
//void sample_landscape(EA& ea) {
//    ea.population().clear();
//    fitness_type_n(ea.population(), get<POPULATION_SIZE>(ea), ea);
//    
//    datafile df(get<ANALYSIS_OUTPUT>(ea));
//    checkpoint_save(ea, df);
//}
//// LIBEA_ANALYSIS_TOOL(sample_landscape, "sample the fitness landscape")

//
//template <typename EA>
//struct landscape_data : public ealib::analysis::unary_function<EA> {
//    static const char* name() { return "landscape_data";}
//    
//    virtual void operator()(EA& ea) {
//        using namespace ealib;
//        using namespace ealib::analysis;
//        
//        datafile df("landscape_data.dat");
//        df.add_field("fitness");
//        
//        for(int i=0; i<get<REPRESENTATION_SIZE>(ea); ++i) {
//            df.add_field("site" + boost::lexical_cast<std::string>(i));
//        }
//        
//        ea.clear();
//        lifecycle::prepare_new(ea);
//        
//        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
//            df.write(static_cast<double>(ealib::fitness(*i,ea))).write_all(i->repr().begin(), i->repr().end()).endl();
//        }
//    }
//};

#endif
