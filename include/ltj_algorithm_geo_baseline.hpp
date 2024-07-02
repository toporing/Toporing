/*
 * ltj_algorithm.hpp
 * Copyright (C) 2020 Author removed for double-blind evaluation
 *
 *
 * This is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */



#ifndef RING_LTJ_ALGORITHM_GEO_BASELINE_HPP
#define RING_LTJ_ALGORITHM_GEO_BASELINE_HPP

#define VERBOSE 1

#include <triple_pattern.hpp>
#include <ring.hpp>
#include <ltj_iterator.hpp>
#include <ltj_iterator_geo.hpp>
#include <veo_adaptive.hpp>
#include <memory>

namespace ring {

    template<class iterator_t = ltj_iterator<ring<>, uint8_t, uint64_t>,
            class geoiterator_t = ltj_iterator_geo<toporing<>, uint8_t, uint64_t>,
             class veo_t = veo::veo_adaptive<iterator_t, util::trait_size> >
    class ltj_algorithm_geo_baseline {

    public:
        typedef uint64_t value_type;
        typedef uint64_t size_type;
        typedef iterator_t ltj_iter_type;
        typedef geoiterator_t ltj_iter_type_geo;
        typedef typename ltj_iter_type::var_type var_type;
        typedef typename ltj_iter_type::ring_type ring_type;
        typedef typename ltj_iter_type::value_type const_type;
        typedef typename ltj_iter_type_geo::georing_type georing_type;
	typedef ltj_iterator_base<var_type, const_type> ltj_iter_base_type;

	typedef std::shared_ptr <ltj_iter_base_type> ltj_iter_base_ptr;
        typedef veo_t veo_type;
        typedef std::unordered_map<var_type, std::vector<ltj_iter_base_ptr>> var_to_iterators_type;
        typedef std::vector<std::pair<var_type, value_type>> tuple_type;
        typedef std::chrono::high_resolution_clock::time_point time_point_type;
    private:
        std::vector<triple_pattern>* m_ptr_triple_patterns;
	std::vector<triple_pattern>* m_ptr_geo_patterns;
	

        veo_type m_veo;
        ring_type* m_ptr_ring;
        georing_type* m_ptr_georing;
        std::vector<ltj_iter_base_ptr> m_iterators;
	std::vector<ltj_iter_base_ptr> m_geo_iterators;
        var_to_iterators_type m_var_to_iterators;
        var_to_iterators_type m_var_to_geo_iterators;
        bool m_is_empty = false;
	size_type m_num_vars;


        void copy(const ltj_algorithm_geo_baseline &o) {
            m_ptr_triple_patterns = o.m_ptr_triple_patterns;
            m_veo = o.m_veo;
            m_ptr_ring = o.m_ptr_ring;
            m_ptr_georing = o.m_ptr_georing;
            m_iterators = o.m_iterators;
            m_var_to_iterators = o.m_var_to_iterators;
            m_geo_iterators = o.m_geo_iterators;
            m_var_to_geo_iterators = o.m_var_to_geo_iterators;
            m_is_empty = o.m_is_empty;
	    m_num_vars = o.m_num_vars;
        }


        inline void add_var_to_iterator(const var_type var, ltj_iter_base_ptr ptr_iterator,
	      var_to_iterators_type &map_var_to_iterators){

            auto it =  map_var_to_iterators.find(var);
            if(it != map_var_to_iterators.end()){
                it->second.push_back(ptr_iterator);
            }else{
                std::vector<ltj_iter_base_ptr> vec = {ptr_iterator};
                map_var_to_iterators.insert({var, vec});
            }
        }


    public:


        ltj_algorithm_geo_baseline() = default;

        ltj_algorithm_geo_baseline(std::vector<triple_pattern>* triple_patterns, ring_type* ring, georing_type*georing){

            m_ptr_triple_patterns = new std::vector<triple_pattern>() ;//triple_patterns;
            m_ptr_geo_patterns = new std::vector<triple_pattern>() ;

            m_ptr_ring = ring;
            m_ptr_georing = georing;
            size_type i = 0;
            m_iterators.reserve(m_ptr_triple_patterns->size());
	    unordered_set<var_type> vars_set;
            for(auto& triple : *triple_patterns){
                //Bulding iterators
		ltj_iter_base_ptr curr_iter;
                if(triple.is_geo()){
                  m_geo_iterators.emplace_back(new ltj_iter_type_geo(&triple, m_ptr_georing));
		  (*m_ptr_geo_patterns).push_back(triple);
		  curr_iter = m_geo_iterators.back();
                }
                else{
                  m_iterators.emplace_back(new ltj_iter_type(&triple, m_ptr_ring));
		  (*m_ptr_triple_patterns).push_back(triple);
		  curr_iter = m_iterators.back();
                }
                if(curr_iter->is_empty()){
                    m_is_empty = true;
                    return;
                }

                //For each variable we add the pointers to its iterators
		var_to_iterators_type &var_to_iter_map = triple.is_geo() ? m_var_to_geo_iterators : m_var_to_iterators;
                if(triple.o_is_variable()){
                    add_var_to_iterator(triple.term_o.value, curr_iter, var_to_iter_map);
		    vars_set.insert(triple.term_o.value);
                }
                if(triple.p_is_variable()){
                    add_var_to_iterator(triple.term_p.value, curr_iter, var_to_iter_map);
		    vars_set.insert(triple.term_p.value);
                }
                if(triple.s_is_variable()){
                    add_var_to_iterator(triple.term_s.value, curr_iter, var_to_iter_map);
		    vars_set.insert(triple.term_s.value);
                }
                ++i;
            }
	    m_num_vars = vars_set.size();
	
            m_veo = veo_type(m_ptr_triple_patterns, &m_iterators, &m_var_to_iterators, m_ptr_ring);

        }

        //! Copy constructor
        ltj_algorithm_geo_baseline(const ltj_algorithm_geo_baseline &o) {
            copy(o);
        }

        //! Move constructor
        ltj_algorithm_geo_baseline(ltj_algorithm_geo_baseline &&o) {
            *this = std::move(o);
        }

        //! Copy Operator=
        ltj_algorithm_geo_baseline &operator=(const ltj_algorithm_geo_baseline &o) {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }

        //! Move Operator=
        ltj_algorithm_geo_baseline &operator=(ltj_algorithm_geo_baseline &&o) {
            if (this != &o) {
                m_ptr_triple_patterns = std::move(o.m_ptr_triple_patterns);
                m_veo = std::move(o.m_veo);
                m_ptr_ring = std::move(o.m_ptr_ring);
                m_ptr_georing = std::move(o.m_ptr_georing);
                m_iterators = std::move(o.m_iterators);
                m_var_to_iterators = std::move(o.m_var_to_iterators);
                m_var_to_geo_iterators = std::move(o.m_var_to_geo_iterators);
                m_is_empty = o.m_is_empty;
		m_num_vars = o.m_num_vars;
            }
            return *this;
        }

        void swap(ltj_algorithm_geo_baseline &o) {
            std::swap(m_ptr_triple_patterns, o.m_ptr_triple_patterns);
            std::swap(m_veo, o.m_veo);
            std::swap(m_ptr_ring, o.m_ptr_ring);
            std::swap(m_ptr_georing, o.m_ptr_georing);
            std::swap(m_iterators, o.m_iterators);
            std::swap(m_var_to_iterators, o.m_var_to_iterators);
            std::swap(m_var_to_geo_iterators, o.m_var_to_geo_iterators);
            std::swap(m_is_empty, o.m_is_empty);
            std::swap(m_num_vars, o.m_num_vars);
        }


        /**
        *
        * @param res               Results
        * @param limit_results     Limit of results
        * @param timeout_seconds   Timeout in seconds
        */
        void join(std::vector<tuple_type> &res,
                  const size_type limit_results = 0, const size_type timeout_seconds = 0){
            if(m_is_empty) return;
            time_point_type start = std::chrono::high_resolution_clock::now();
            tuple_type t(m_num_vars);



            search(0, t, res, start,
		//two_ready, one_ready,sim , bind_value, bind_state, 
		limit_results, timeout_seconds);
        };

        void join_v2(std::vector<tuple_type> &res,
                  const size_type limit_results = 0, const size_type timeout_seconds = 0){
            if(m_is_empty) return;
            time_point_type start = std::chrono::high_resolution_clock::now();
            tuple_type t(m_num_vars);
            search_v2(0, t, res, start, limit_results, timeout_seconds);
        };

        void join_ramas(std::vector<tuple_type> &res, size_type &ramas,
                  const size_type limit_results = 0, const size_type timeout_seconds = 0){
            if(m_is_empty) return;
            time_point_type start = std::chrono::high_resolution_clock::now();
            tuple_type t(m_num_vars);
            search_ramas(0, t, res, start, ramas, limit_results, timeout_seconds);
        };


        /**
         *
         * @param j                 Index of the variable
         * @param tuple             Tuple of the current search
         * @param res               Results
         * @param start             Initial time to check timeout
         * @param limit_results     Limit of results
         * @param timeout_seconds   Timeout in seconds
         */
        bool search(const size_type j, tuple_type &tuple, std::vector<tuple_type> &res,
                    const time_point_type start,
		    /*std::unordered_set<ltj_iter_base_ptr> &two_ready,
		    std::unordered_set<ltj_iter_base_ptr> &one_ready,
		    std::unordered_set<ltj_iter_base_ptr> &sim,
		    std::unordered_map<var_type, value_type> &bind_value,
		    std::unordered_map<ltj_iter_base_ptr, int> &bind_state,*/
                    const size_type limit_results = 0, const size_type timeout_seconds = 0){

            //(Optional) Check timeout
            if(timeout_seconds > 0){
                time_point_type stop = std::chrono::high_resolution_clock::now();
                auto sec = std::chrono::duration_cast<std::chrono::seconds>(stop-start).count();
                if(sec > timeout_seconds){
		  //cout << "TIMED OUT" <<endl;
		  return false;
		}
            }

            //(Optional) Check limit
            if(limit_results > 0 && res.size() == limit_results) return false;

	    //cout << "j = "<<j<<endl;
            if(j == m_veo.size()){
		// initialize map with all variables of all geo iterators
		// to catch variables in sim

	        std::unordered_set<ltj_iter_base_ptr> two_ready, one_ready, sim;

	        std::unordered_map<var_type, value_type> bind_value;
		std::unordered_map<ltj_iter_base_ptr, int> bind_state;

		for(int i=0; i<j; i++){

		  var_type binded_var = tuple[i].first;
		  value_type var_value = tuple[i].second;

		  bind_value[binded_var] = var_value;

		  for(ltj_iter_base_ptr iter : m_var_to_geo_iterators[binded_var] ){
		    iter->down(binded_var, var_value);
		  //  bind_state[iter]++;
		  }

		}

		for(auto &p : m_var_to_geo_iterators){
		  std::vector<ltj_iter_base_ptr> &iters = p.second;
		  for(ltj_iter_base_ptr iter : iters){
		    //dirty hack
		    bind_state[iter] = iter->is_s_binded()+iter->is_o_binded();
		  }
		}
		

		for(auto &p : bind_state){
		  ltj_iter_base_ptr iter = p.first;
		  int state = p.second;

		  //cerr << "state is "<<state<<endl;
		  if(state == 0){
		    sim.insert(iter);
		  }else if(state == 1){
		    one_ready.insert(iter);
		  }else if(state == 2){
		    two_ready.insert(iter);
		  }

		}

		bool ok = filter_and_extend(j, tuple, res,
				  two_ready, one_ready,
				  sim, bind_value,
				  start,
				  limit_results, timeout_seconds);
		
		if ( !ok ) return false;

		//sim.clear();
		//one_ready.clear();
		//two_ready.clear();
		for(int i=0; i<j; i++){

		  var_type binded_var = tuple[i].first;

		  for(ltj_iter_base_ptr iter : m_var_to_geo_iterators[binded_var] ){
		    iter->up(binded_var);
		  //  bind_state[iter]++;
		  }

		}
                //Report results
                //res.emplace_back(tuple);
		/*
		if(true){
                std::cout << "Add result" << std::endl;
                for(const auto &dat : tuple){
                    std::cout << "{" << (uint64_t) dat.first << "=" << dat.second << "} ";
                }
                std::cout << std::endl;
		}
		*/
		
            }else{
                var_type x_j = m_veo.next();
                //std::cout << "Variable: " << (uint64_t) x_j << std::endl;
                std::vector<ltj_iter_base_ptr>& itrs = m_var_to_iterators[x_j];
                bool ok;
                if(itrs.size() == 1 && itrs[0]->in_last_level()) {//Lonely variables
                    //std::cout << "Seeking (last level)" << std::endl;
                    auto results = itrs[0]->seek_all(x_j);
                    //std::cout << "Results: " << results.size() << std::endl;
                    //std::cout << "Seek (last level): (" << (uint64_t) x_j << ": size=" << results.size() << ")" <<std::endl;
                    for (const auto &c : results) {
                        //std::cout << "Seek (bucle): (" << (uint64_t) x_j << ": " << c << ")" <<std::endl;
		        if(c==0) continue;
                        //1. Adding result to tuple
                        tuple[j] = {x_j, c};
                        //2. Going down in the trie by setting x_j = c (\mu(t_i) in paper)
                        itrs[0]->down(x_j, c);
                        m_veo.down();
                        //2. Search with the next variable x_{j+1}
                        ok = search(j + 1, tuple, res, start, 
			    //two_ready, one_ready,sim , bind_value, bind_state, 
			    limit_results, timeout_seconds);
                        if(!ok) return false;
                        //4. Going up in the trie by removing x_j = c
                        itrs[0]->up(x_j);
                        m_veo.up();
                    }
                }else {
                    value_type c = seek(x_j);
                    //std::cout << "Seek (init): (" << (uint64_t) x_j << ": " << c << ")" <<std::endl;
                    while (c != 0) { //If empty c=0
                        //1. Adding result to tuple
                        tuple[j] = {x_j, c};
                        //2. Going down in the tries by setting x_j = c (\mu(t_i) in paper)
                        for (ltj_iter_base_ptr iter : itrs) {
                            iter->down(x_j, c);
                        }
                        m_veo.down();
                        //3. Search with the next variable x_{j+1}
                        ok = search(j + 1, tuple, res, start, 
			    //two_ready, one_ready,sim , bind_value, bind_state, 
			    limit_results, timeout_seconds);
                        if(!ok) return false;
                        //4. Going up in the tries by removing x_j = c
                        for (ltj_iter_base_ptr iter : itrs) {
                            iter->up(x_j);
                        }
                        m_veo.up();
                        //5. Next constant for x_j
                        c = seek(x_j, c + 1);
                        //std::cout << "Seek (bucle): (" << (uint64_t) x_j << ": " << c << ")" <<std::endl;
                    }
                }
                m_veo.done();
            }
            return true;
        };


	void bind_all_geo_iters(var_type var, value_type value){

	  for(auto &iter : m_var_to_geo_iterators[var]){
	    iter->down(var, value);
	  }
	}


	void unbind_all_geo_iters(var_type var){

	  for(auto &iter : m_var_to_geo_iterators[var]){
	    iter->up(var);
	  }
	}

	void update_geo_iters_state(var_type var,std::unordered_set<ltj_iter_base_ptr> &two_ready,
	    std::unordered_set<ltj_iter_base_ptr> &one_ready,std::unordered_set<ltj_iter_base_ptr> &sim){

	  for(auto &iter : m_var_to_geo_iterators[var]){
	    // could be optimized, although set is small 
	    // (stop after one erase is valid?)
	    two_ready.erase(iter);
	    one_ready.erase(iter);
	    sim.erase(iter);

	    int state = iter->is_s_binded() + iter->is_o_binded();

	    if(state == 2){
	      two_ready.insert(iter);
	    }else if (state == 1){
	      one_ready.insert(iter);
	    }else if(state == 0){
	      sim.insert(iter);
	    }
	  }
	}

        bool filter_and_extend(const size_type j, tuple_type &tuple, std::vector<tuple_type> &res,
		    std::unordered_set<ltj_iter_base_ptr> &two_ready, std::unordered_set<ltj_iter_base_ptr> &one_ready,
		    std::unordered_set<ltj_iter_base_ptr> &sim, std::unordered_map<var_type, value_type> &bind_value, 
                    const time_point_type start,
                    const size_type limit_results = 0, const size_type timeout_seconds = 0){
	    //cerr << "filter and extend" << endl;
	    //cerr << "2-ready size "<< two_ready.size()<<endl;
	    //cerr << "ready size "<< one_ready.size()<<endl;
	    //cerr << "sim size "<< sim.size()<<endl;
            if(timeout_seconds > 0){
                time_point_type stop = std::chrono::high_resolution_clock::now();
                auto sec = std::chrono::duration_cast<std::chrono::seconds>(stop-start).count();
                if(sec > timeout_seconds) return false;
            }

            //(Optional) Check limit
            if(limit_results > 0 && res.size() == limit_results) return false;


	    if( !two_ready.empty() ){
	      ltj_iter_base_ptr iter= *two_ready.begin();
	      
	      if(iter->check()){
	        two_ready.erase(iter);
		bool ok = filter_and_extend(j, tuple, res,
				  two_ready, one_ready,
				  sim, bind_value,
				  start,
				  limit_results, timeout_seconds);

		if(!ok) return false;

		two_ready.insert(iter);
	      }

	      return true;

	    }else if( !one_ready.empty() or !sim.empty() ){
	      // get one iterator and assign a value to a variable

	      //get most restrictive iterator instead?
	      //
	      
	      bool iter_from_one = !one_ready.empty();

	      ltj_iter_base_ptr iter= iter_from_one ? *one_ready.begin() : *sim.begin();

	      var_type next_var = (iter->is_s_binded()) ? iter->get_object_variable() :
							  iter->get_subject_variable();

	      auto results = iter->seek_all(next_var);

	      if (iter_from_one) one_ready.erase(iter);
	      else sim.erase(iter);
	      for(const auto &c: results){

		bind_all_geo_iters(next_var, c);
		update_geo_iters_state(next_var, two_ready,
		    one_ready, sim);

		tuple[j] = {next_var, c};

		bool ok = filter_and_extend( j+1, tuple, res,
				  two_ready, one_ready,
				  sim, bind_value,
				  start,
				  limit_results, timeout_seconds);
		
		if(!ok) return false;

		unbind_all_geo_iters(next_var);
		update_geo_iters_state(next_var, two_ready,
		    one_ready, sim);
	      }
	      if(iter_from_one) one_ready.insert(iter);
	      else sim.insert(iter);

	      return true;
	    }
	    /*else if( !sim.empty() ){

	      cerr << "no remaining iters in 2-ready or ready"<< endl;
	      // simular????
	    }*/

	    
	    // no variables left
	    
	    else if(j!=tuple.size()){
	      //cerr << "ERROR: extension finished with unbinded variables"<<endl;
	      return false;
	    }

	    //everything is binded
	    res.push_back(tuple);
	  
	    /*
	    std::cout <<"add result ";
	    for(const auto &dat : tuple){
                    std::cout << "{" << (uint64_t) dat.first << "=" << dat.second << "} ";
            }
	    std::cout <<endl;
	    */

	    return true;

	}
        bool search_v2(const size_type j, tuple_type &tuple, std::vector<tuple_type> &res,
                    const time_point_type start,
                    const size_type limit_results = 0, const size_type timeout_seconds = 0){

            //(Optional) Check timeout
            if(timeout_seconds > 0){
                time_point_type stop = std::chrono::high_resolution_clock::now();
                auto sec = std::chrono::duration_cast<std::chrono::seconds>(stop-start).count();
                if(sec > timeout_seconds) return false;
            }

            //(Optional) Check limit
            if(limit_results > 0 && res.size() == limit_results) return false;

            if(j == m_veo.size()){
	        
	        std::set<ltj_iter_base_ptr> two_ready, one_ready, sim;

	        std::map<var_type, value_type> bind_value;
		std::map<ltj_iter_base_ptr, int> bind_state;

		// initialize map with all variables of all geo iterators
		// to catch variables in sim

		for(int i=0; i<j; i++){

		  var_type binded_var = tuple[i].first;
		  value_type var_value = tuple[i].second;

		  bind_value[binded_var] = var_value;

		  for(ltj_iter_base_ptr iter : m_var_to_geo_iterators[binded_var] ){
		    iter->down(binded_var, var_value);
		  //  bind_state[iter]++;
		  }

		}

		for(auto &p : m_var_to_geo_iterators){
		  std::vector<ltj_iter_base_ptr> &iters = p.second;
		  for(ltj_iter_base_ptr iter : iters){
		    //dirty hack
		    bind_state[iter] = iter->is_s_binded()+iter->is_o_binded();
		  }
		}
		

		for(auto &p : bind_state){
		  ltj_iter_base_ptr iter = p.first;
		  int state = p.second;

		  if(state == 0){
		    sim.insert(iter);
		  }else if(state == 1){
		    one_ready.insert(iter);
		  }else if(state == 2){
		    two_ready.insert(iter);
		  }

		}

		bool ok = filter_and_extend(j, tuple, res,
				  two_ready, one_ready,
				  sim, bind_value,
				  start,
				  limit_results, timeout_seconds);
		
		if ( !ok ) return false;

		for(int i=0; i<j; i++){

		  var_type binded_var = tuple[i].first;

		  for(ltj_iter_base_ptr iter : m_var_to_geo_iterators[binded_var] ){
		    iter->up(binded_var);
		  //  bind_state[iter]++;
		  }

		}
                //res.emplace_back(tuple);
                /*//std::cout << "Add result" << std::endl;
                for(const auto &dat : tuple){
                    //std::cout << "{" << (uint64_t) dat.first << "=" << dat.second << "} ";
                }
                //std::cout << std::endl;*/
            }else{
                var_type x_j = m_veo.next();
                ////std::cout << "Variable: " << (uint64_t) x_j << std::endl;
                std::vector<ltj_iter_base_ptr>& itrs = m_var_to_iterators[x_j];
                bool ok;
                if(itrs.size() == 1 && itrs[0]->in_last_level()) {//Lonely variables
                    ////std::cout << "Seeking (last level)" << std::endl;
                    value_type c = itrs[0]->seek_last(x_j);
                    //auto results = itrs[0]->seek_all(x_j);
                    ////std::cout << "Results: " << results.size() << std::endl;
                    ////std::cout << "Seek (last level): (" << (uint64_t) x_j << ": size=" << results.size() << ")" <<std::endl;
                    while (c != 0) { //If empty c=0
                        //1. Adding result to tuple
                        tuple[j] = {x_j, c};
                        //2. Going down in the trie by setting x_j = c (\mu(t_i) in paper)
                        itrs[0]->down(x_j, c);
                        m_veo.down();
                        //2. Search with the next variable x_{j+1}
                        ok = search_v2(j + 1, tuple, res, start, limit_results, timeout_seconds);
                        if(!ok) return false;
                        //4. Going up in the trie by removing x_j = c
                        itrs[0]->up(x_j);
                        m_veo.up();

                        c = itrs[0]->seek_last_next(x_j);
                    }
                }else {
                    value_type c = seek(x_j);
                    ////std::cout << "Seek (init): (" << (uint64_t) x_j << ": " << c << ")" <<std::endl;
                    while (c != 0) { //If empty c=0
                        //1. Adding result to tuple
                        tuple[j] = {x_j, c};
                        //2. Going down in the tries by setting x_j = c (\mu(t_i) in paper)
                        for (ltj_iter_base_ptr iter : itrs) {
                            iter->down(x_j, c);
                        }
                        m_veo.down();
                        //3. Search with the next variable x_{j+1}
                        ok = search_v2(j + 1, tuple, res, start, limit_results, timeout_seconds);
                        if(!ok) return false;
                        //4. Going up in the tries by removing x_j = c
                        for (ltj_iter_base_ptr iter : itrs) {
                            iter->up(x_j);
                        }
                        m_veo.up();
                        //5. Next constant for x_j
                        c = seek(x_j, c + 1);
                        ////std::cout << "Seek (bucle): (" << (uint64_t) x_j << ": " << c << ")" <<std::endl;
                    }
                }
                m_veo.done();
            }
            return true;
        };


        bool search_ramas(const size_type j, tuple_type &tuple, std::vector<tuple_type> &res,
                    const time_point_type start, size_type &ramas,
                    const size_type limit_results = 0, const size_type timeout_seconds = 0){

            //(Optional) Check timeout
            if(timeout_seconds > 0){
                time_point_type stop = std::chrono::high_resolution_clock::now();
                auto sec = std::chrono::duration_cast<std::chrono::seconds>(stop-start).count();
                if(sec > timeout_seconds) return false;
            }

            //(Optional) Check limit
            if(limit_results > 0 && res.size() == limit_results) return false;

            if(j == m_veo.nolonely_size()){
                ++ramas;
            }
            if(j == m_veo.size()){
                //Report results
                res.emplace_back(tuple);
                /*//std::cout << "Add result" << std::endl;
                for(const auto &dat : tuple){
                    //std::cout << "{" << (uint64_t) dat.first << "=" << dat.second << "} ";
                }
                //std::cout << std::endl;*/
            }else{
                var_type x_j = m_veo.next();
                ////std::cout << "Variable: " << (uint64_t) x_j << std::endl;
                std::vector<ltj_iter_base_ptr>& itrs = m_var_to_iterators[x_j];
                bool ok;
                if(itrs.size() == 1 && itrs[0]->in_last_level()) {//Lonely variables
                    ////std::cout << "Seeking (last level)" << std::endl;
                    auto results = itrs[0]->seek_all(x_j);
                    ////std::cout << "Results: " << results.size() << std::endl;
                    ////std::cout << "Seek (last level): (" << (uint64_t) x_j << ": size=" << results.size() << ")" <<std::endl;
                    for (const auto &c : results) {
                        //1. Adding result to tuple
                        tuple[j] = {x_j, c};
                        //2. Going down in the trie by setting x_j = c (\mu(t_i) in paper)
                        itrs[0]->down(x_j, c);
                        m_veo.down();
                        //2. Search with the next variable x_{j+1}
                        ok = search_ramas(j + 1, tuple, res, start, ramas,limit_results, timeout_seconds);
                        if(!ok) return false;
                        //4. Going up in the trie by removing x_j = c
                        itrs[0]->up(x_j);
                        m_veo.up();
                    }
                }else {
                    value_type c = seek(x_j);
                    ////std::cout << "Seek (init): (" << (uint64_t) x_j << ": " << c << ")" <<std::endl;
                    while (c != 0) { //If empty c=0
                        //1. Adding result to tuple
                        tuple[j] = {x_j, c};
                        //2. Going down in the tries by setting x_j = c (\mu(t_i) in paper)
                        for (ltj_iter_base_ptr iter : itrs) {
                            iter->down(x_j, c);
                        }
                        m_veo.down();
                        //3. Search with the next variable x_{j+1}
                        ok = search_ramas(j + 1, tuple, res, start, ramas,limit_results, timeout_seconds);
                        if(!ok) return false;
                        //4. Going up in the tries by removing x_j = c
                        for (ltj_iter_base_ptr iter : itrs) {
                            iter->up(x_j);
                        }
                        m_veo.up();
                        //5. Next constant for x_j
                        c = seek(x_j, c + 1);
                        ////std::cout << "Seek (bucle): (" << (uint64_t) x_j << ": " << c << ")" <<std::endl;
                    }
                }
                m_veo.done();
            }
            return true;
        };

        /**
         *
         * @param x_j   Variable
         * @param c     Constant. If it is unknown the value is -1
         * @return      The next constant that matches the intersection between the triples of x_j.
         *              If the intersection is empty, it returns 0.
         */

        value_type seek(const var_type x_j, value_type c=-1){
            std::vector<ltj_iter_base_ptr>& itrs = m_var_to_iterators[x_j];
            value_type c_i, c_prev = 0, i = 0, n_ok = 0;
            while (true){
                //Compute leap for each triple that contains x_j
                if(c == -1){
                    c_i = itrs[i]->leap(x_j);
                }else{
                    c_i = itrs[i]->leap(x_j, c);
                }
                if(c_i == 0) return 0; //Empty intersection
                n_ok = (c_i == c_prev) ? n_ok + 1 : 1;
                if(n_ok == itrs.size()) return c_i;
                c = c_prev = c_i;
                i = (i+1 == itrs.size()) ? 0 : i+1;
            }
        }



        void print_veo(std::unordered_map<uint8_t, std::string> &ht){
            std::cout << "veo: " << std::endl;
            for(uint64_t j = 0; j < m_veo.size(); ++j){
                std::cout << "?" << ht[m_veo.next()] << " ";
            }
            std::cout << std::endl;
        }

        void print_query(std::unordered_map<uint8_t, std::string> &ht ){
            std::cout << "Query: " << std::endl;
            for(size_type i = 0; i <  m_ptr_triple_patterns->size(); ++i){
                m_ptr_triple_patterns->at(i).print(ht);
                if(i < m_ptr_triple_patterns->size()-1){
                    std::cout << " . ";
                }
            }
            std::cout << std::endl;
        }

        void print_results(std::vector<tuple_type> &res, std::unordered_map<uint8_t, std::string> &ht, sdsl::int_vector<> &inv_perm){
            std::cout << "Results: " << std::endl;
            uint64_t i = 1;
            for(tuple_type &tuple :  res){
                std::cout << "[" << i << "]: ";
                for(std::pair<var_type, value_type> &pair : tuple){
                    std::cout << "?" << ht[pair.first] << "=" << inv_perm[pair.second] << " ";
                }
                std::cout << std::endl;
                ++i;
            }
        }


    };
}

#endif //RING_LTJ_ALGORITHM_GEO_HPP
