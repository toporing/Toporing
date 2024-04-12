/*
 * ltj_iterator.hpp
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

#ifndef RING_LTJ_ITERATOR_GEO_HPP
#define RING_LTJ_ITERATOR_GEO_HPP

#define VERBOSE 0

#include <configuration.hpp>
#include <ltj_iterator_base.hpp>
#include <toporing.hpp>
#include <triple_pattern.hpp>
namespace ring {

    template<class georing_t, class var_t, class cons_t>
    class ltj_iterator_geo : public ltj_iterator_base< var_t, cons_t>{

    public:
        typedef cons_t value_type;
        typedef var_t var_type;
        typedef georing_t georing_type;
        typedef uint64_t size_type;
        // typedef  wt_range_iterator<typename ring_type::bwt_type::wm_type> wt_so_iterator_type;
        // typedef  wt_range_iterator<typename ring_type::bwt_p_type::wm_type> wt_p_iterator_type;
        //std::vector<value_type> leap_result_type;

    private:
        triple_pattern m_ptr_triple_pattern;
        georing_type *m_ptr_georing; //TODO: should be const
        size_type m_level = 0;
	       size_type binding_s = 0;
	        size_type binding_o = 0;
          bool m_is_empty = false;
          std::array<state_type, 3> m_state;
        //std::stack<state_type> m_states;
    public:

      void copy(const ltj_iterator_geo &o) {
          m_is_empty = o.m_is_empty;
          m_ptr_triple_pattern = o.m_ptr_triple_pattern;
          m_ptr_georing = o.m_ptr_georing;
          // m_intervals = o.m_intervals;
          // m_so_last_iterator = o.m_so_last_iterator;
          // m_p_last_iterator = o.m_p_last_iterator;
          m_state = o.m_state;

          binding_s = o.binding_s;
          binding_o = o.binding_o;
          // m_consts = o.m_consts;
      }

    	ltj_iterator_geo(const triple_pattern *triple, georing_type *georing) :
    		m_ptr_triple_pattern(*triple) , m_ptr_georing(georing){
	    if(m_ptr_triple_pattern.is_geo() && m_ptr_triple_pattern.term_p.value == TOPO_CONTAINED_IN){

	      std::swap(m_ptr_triple_pattern.term_s, m_ptr_triple_pattern.term_o);
	      m_ptr_triple_pattern.term_p.value = TOPO_CONTAINS;
	    }

	    if(m_ptr_triple_pattern.is_geo() && m_ptr_triple_pattern.term_p.value == TOPO_NOT_CONTAINED_IN){

	      std::swap(m_ptr_triple_pattern.term_s, m_ptr_triple_pattern.term_o);
	      m_ptr_triple_pattern.term_p.value = TOPO_NOT_CONTAINS;
	    }


	    if(m_ptr_triple_pattern.is_geo() && m_ptr_triple_pattern.term_p.value == TOPO_TOUCHES &&
		m_ptr_triple_pattern.s_is_variable() && m_ptr_triple_pattern.o_is_variable() &&
		m_ptr_triple_pattern.term_s.value == m_ptr_triple_pattern.term_o.value){
	      //invalid: nothing can touch self
	      m_ptr_triple_pattern.term_p.value = TOPO_INVALID;
	    }
            //Init current values and intervals according to the triple
            if (!m_ptr_triple_pattern.s_is_variable() && !m_ptr_triple_pattern.p_is_variable()
                && !m_ptr_triple_pattern.o_is_variable()) {
                //S->O->P to avoid forward steps

                //Interval in S

                m_state[0] = s;
                binding_s = m_ptr_triple_pattern.term_s.value;

                m_state[1] = o;
                binding_o = m_ptr_triple_pattern.term_o.value;


                m_state[2] = p;

                m_level = 3;

            } else if (!m_ptr_triple_pattern.s_is_variable() && !m_ptr_triple_pattern.p_is_variable()) {
                //P->S to avoid forward steps


                m_state[0] = p;


                m_state[1] = s;
                binding_s = m_ptr_triple_pattern.term_s.value;

                m_level = 2;

            } else if (!m_ptr_triple_pattern.p_is_variable() && !m_ptr_triple_pattern.o_is_variable()) {
                //O->P to avoid forward steps

                m_state[0] = o;
                binding_o = m_ptr_triple_pattern.term_o.value;


                m_state[1] = p;

                m_level = 2;

            } else if (!m_ptr_triple_pattern.s_is_variable() && !m_ptr_triple_pattern.o_is_variable()) {
                //S->O to avoid forward steps

                //Interval in S

                m_state[0] = s;
                binding_s = m_ptr_triple_pattern.term_s.value;


                m_state[1] = o;
                binding_o = m_ptr_triple_pattern.term_o.value;

                m_level = 2;

            } else if (!m_ptr_triple_pattern.s_is_variable()) {

                //Interval in S

                m_state[0] = s;
                binding_s = m_ptr_triple_pattern.term_s.value;


                m_level = 1;

            } else if (!m_ptr_triple_pattern.p_is_variable()) {

                //Interval in P

                m_state[0] = p;

                m_level = 1;

            } else if (!m_ptr_triple_pattern.o_is_variable()) {

                //Interval in O

                m_state[0] = o;
                binding_o = m_ptr_triple_pattern.term_o.value;

                m_level = 1;


            }
        }

        ltj_iterator_geo(const ltj_iterator_geo &o) {
            copy(o);
        }

        //! Move constructor
        ltj_iterator_geo(ltj_iterator_geo &&o) {
            *this = std::move(o);
        }

        ltj_iterator_geo &operator=(const ltj_iterator_geo &o) {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }

        //! Move Operator=
        ltj_iterator_geo &operator=(ltj_iterator_geo &&o) {
            if (this != &o) {
                m_ptr_triple_pattern = std::move(o.m_ptr_triple_pattern);
                m_ptr_georing = std::move(o.m_ptr_georing);
                // m_intervals = std::move(o.m_intervals);
                // m_so_last_iterator = std::move(o.m_so_last_iterator);
                // m_p_last_iterator = std::move(o.m_p_last_iterator);
                // m_consts = std::move(o.m_consts);
                m_state = std::move(o.m_state);
                m_level = o.m_level;
                m_is_empty = o.m_is_empty;
                binding_s = o.binding_s;
                binding_o = o.binding_o;
            }
            return *this;
        }

        void swap(ltj_iterator_geo &o) {
            // m_bp.swap(bp_support.m_bp); use set_vector to set the supported bit_vector
            std::swap(m_ptr_triple_pattern, o.m_ptr_triple_pattern);
            std::swap(m_ptr_georing, o.m_ptr_georing);
            // std::swap(m_intervals, o.m_intervals);
            // std::swap(m_consts, o.m_consts);
            // std::swap(m_so_last_iterator, o.m_so_last_iterator);
            // std::swap(m_p_last_iterator, o.m_p_last_iterator);
            std::swap(m_state, o.m_state);
            std::swap(m_level, o.m_level);
            std::swap(binding_s, o.binding_s);
            std::swap(binding_o, o.binding_o);
            std::swap(m_is_empty, o.m_is_empty);
        }



	bool is_contains(){
		return !m_ptr_triple_pattern.term_p.is_variable && m_ptr_triple_pattern.term_p.value == TOPO_CONTAINS;
	}

	bool is_contained_in(){
		return !m_ptr_triple_pattern.term_p.is_variable && m_ptr_triple_pattern.term_p.value == TOPO_CONTAINED_IN;
	}

	bool is_touches(){
		return !m_ptr_triple_pattern.term_p.is_variable && m_ptr_triple_pattern.term_p.value == TOPO_TOUCHES;
	}

	bool is_empty() override{
		return false;
	}

	bool check_binded(state_type state){

		if( m_level == 0 ) return false;

		if( m_level >= 1 && m_state[0]==state) return true;
		if( m_level >= 2 && m_state[1]==state) return true;
		if( m_level >= 3 && m_state[2]==state) return true;

		return false;
	}

	bool is_s_binded(){
		return check_binded(s);
	}


	bool is_p_binded(){
		return check_binded(p);
	}

	bool is_o_binded(){
		return check_binded(o);
	}


	void down(var_type var, size_type c) override {
		if(m_level>2) return;

		if(is_variable_subject(var) and !is_s_binded()){
      //cout << "binding s to "<<c<<endl;
			binding_s = c;
			m_state[m_level] = s;
		}else if (is_variable_object(var) and !is_o_binded()){
      //cout << "binding o to "<<c<<endl;
			binding_o = c;
			m_state[m_level] = o;
		}else{
			m_state[m_level] = p;
		}
		++m_level;
		return;
	}

	void down(var_type var, size_type c, size_type k) override {
		down(var, c);
		return;
	}

	void up (var_type var) override {
	    //cout << "unbinding"<<endl;
		if ( m_level == 0 ) return;
		--m_level;
	}

	/*
	value_type leap(var_type var) override {

		if( is_contains() ) return leap_contains(var);
		if( is_contained_in() ) return leap_contained_in(var);
		if( is_touches() ) return leap_touches(var);

		return 0;

	}
	*/

	inline value_type min_op_q(value_type binding){
	  switch (m_ptr_triple_pattern.term_p.value) {
	    case TOPO_CONTAINS:
	      return m_ptr_georing->min_contains_q(binding);
	    case TOPO_TOUCHES:
	      return m_ptr_georing->min_touches_q(binding);
	    case TOPO_NOT_CONTAINS:
	      return m_ptr_georing->min_not_contains_q(binding);
	    case TOPO_INVALID:
	      return 0;
	    default:
	      cout << "UNDEFINED TYPE" <<endl;
	      return 0;
	  }
	}

	inline value_type min_inv_op_q(value_type binding){
	  switch (m_ptr_triple_pattern.term_p.value) {
	    case TOPO_CONTAINS:
	      return m_ptr_georing->min_contained_in_q(binding);
	    case TOPO_TOUCHES:
	      return m_ptr_georing->min_touches_q(binding);
	    case TOPO_NOT_CONTAINS:
	      return m_ptr_georing->min_not_contained_in_q(binding);
	    case TOPO_INVALID:
	      return 0;
	    default:
	      cerr << "UNDEFINED TYPE" <<endl;
	      return 0;
	  }
	}


	inline value_type min_op(value_type binding){
	  switch (m_ptr_triple_pattern.term_p.value) {
	    case TOPO_CONTAINS:
	      return m_ptr_georing->min_contains(binding);
	    case TOPO_TOUCHES:
	      return m_ptr_georing->min_touches(binding);
	    case TOPO_NOT_CONTAINS:
	      return m_ptr_georing->min_not_contains(binding);
	    case TOPO_INVALID:
	      return 0;
	    default:
	      cerr << "UNDEFINED TYPE" <<endl;
	      return 0;
	  }
	}

	inline value_type min_inv_op(value_type binding){
	  switch (m_ptr_triple_pattern.term_p.value) {
	    case TOPO_CONTAINS:
	      return m_ptr_georing->min_contained_in(binding);
	    case TOPO_TOUCHES:
	      return m_ptr_georing->min_touches(binding);
	    case TOPO_NOT_CONTAINS:
	      return m_ptr_georing->min_not_contained_in(binding);
	    case TOPO_INVALID:
	      return 0;
	    default:
	      cerr << "UNDEFINED TYPE" <<endl;
	      return 0;
	  }
	}

	
	value_type leap(var_type var) override {
		if (is_s_binded() && is_variable_object(var)){
			return min_op_q(binding_s);
		}else if( is_o_binded() && is_variable_subject(var) ){
			// backwards in
			return min_inv_op_q(binding_o);
		}else if( !is_s_binded() && is_variable_subject(var)){
			// return smallest that op another
			return min_inv_op(1);
		}else if( !is_o_binded() && is_variable_object(var)){
			// return smallest that op another
			return min_op(1);
		}
		return 0;

	}

/*
	value_type leap_contains(var_type var)  {
    cout << "leap contains" <<endl;
		if (is_s_binded() && is_variable_object(var)){
			return m_ptr_georing->min_op_q(binding_s);
		}else if( is_o_binded() && is_variable_subject(var) ){
			// backwards in
			return m_ptr_georing->min_inv_op_q(binding_o);
		}else if( !is_s_binded() && is_variable_subject(var)){
			// return smallest that op another
			return m_ptr_georing->min_inv_op(1);
		}else if( !is_o_binded() && is_variable_object(var)){
			// return smallest that op another
			return m_ptr_georing->min_op(1);
		}
		return 0;

	}

	value_type leap_not_contains(var_type var)  {
    cout << "leap contains" <<endl;
		if (is_s_binded() && is_variable_object(var)){
			return m_ptr_georing->min_not_contains_q(binding_s);
		}else if( is_o_binded() && is_variable_subject(var) ){
			// backwards in
			return m_ptr_georing->min_not_contained_in_q(binding_o);
		}else if( !is_s_binded() && is_variable_subject(var)){
			// return smallest that contains another
			return m_ptr_georing->min_not_contained_in(1);
		}else if( !is_o_binded() && is_variable_object(var)){
			// return smallest that contains another
			return m_ptr_georing->min_not_contains(1);
		}
		return 0;

	}

	value_type leap_contained_in(var_type var)  {
    cout << "leap contained" <<endl;
		if (is_s_binded() && is_variable_object(var)){
			return m_ptr_georing->min_contained_in_q(binding_s);
		}else if( is_o_binded() && is_variable_subject(var) ){
			// backwards in
			return m_ptr_georing->min_contains_q(binding_o);
		}else if( !is_s_binded() && is_variable_subject(var)){
			// return smallest that contains another
			return m_ptr_georing->min_contains(1);
		}else if( !is_o_binded() && is_variable_object(var)){
			// return smallest that contains another
			return m_ptr_georing->min_contained_in(1);
		}
		return 0;

	}

	value_type leap_touches(var_type var)  {

		return 0;

	}
*/
	
	value_type leap(var_type var, size_type c) override{
		if( is_s_binded() && is_variable_object(var) ){
			return op(binding_s, c);
		}else if( is_o_binded() && is_variable_subject(var) ){
			// backwards in
			return inv_op(binding_o, c);
		}	if( !is_s_binded() && is_variable_subject(var) ){
  			return min_inv_op(c);
		}else if( !is_o_binded() && is_variable_object(var) ){
			// backwards in
			return min_op( c);
		}

		return 0;

		// std :: cerr << "WARNING:  undefined path for in\n";

	}
/*
	value_type leap(var_type var, size_type c) override {

		if( is_contains() ) return leap_contains(var, c);
		if( is_contained_in() ) return leap_contained_in(var, c);
		if( is_touches() ) return leap_touches(var, c);

		return 0;

	}
*/
	inline value_type op(value_type binding, size_type c){
	  switch (m_ptr_triple_pattern.term_p.value) {
	    case TOPO_CONTAINS:
	      return m_ptr_georing->contains(binding,c);
	    case TOPO_TOUCHES:
	      return m_ptr_georing->touches(binding,c);
	    case TOPO_NOT_CONTAINS:
	      return m_ptr_georing->not_contains(binding,c);
	    case TOPO_INVALID:
	      return 0;
	    default:
	      cout << "UNDEFINED TYPE" <<endl;
	      return 0;
	  }
	}

	inline bool check(){
	  switch (m_ptr_triple_pattern.term_p.value) {
	    case TOPO_CONTAINS:
	      return m_ptr_georing->contains(binding_s, binding_o)==binding_o;
	    case TOPO_TOUCHES:
	      return m_ptr_georing->touches(binding_s,binding_o)==binding_o;
	    case TOPO_NOT_CONTAINS:
	      return m_ptr_georing->not_contains(binding_s,binding_o)==binding_o;
	    case TOPO_INVALID:
	      return 0;
	    default:
	      cout << "UNDEFINED TYPE" <<endl;
	      return 0;
	  }
	}

	inline value_type inv_op(value_type binding, size_type c){
	  switch (m_ptr_triple_pattern.term_p.value) {
	    case TOPO_CONTAINS:
	      return m_ptr_georing->contained_in(binding, c);
	    case TOPO_TOUCHES:
	      return m_ptr_georing->touches(binding, c);
	    case TOPO_NOT_CONTAINS:
	      return m_ptr_georing->not_contained_in(binding, c);
	    case TOPO_INVALID:
	      return 0;
	    default:
	      cout << "UNDEFINED TYPE" <<endl;
	      return 0;
	  }
	}
/*
	value_type leap_contains(var_type var, size_type c){
		if( is_s_binded() && is_variable_object(var) ){
			return m_ptr_georing->op(binding_s, c);
		}else if( is_o_binded() && is_variable_subject(var) ){
			// backwards in
			return m_ptr_georing->contained_in(binding_o, c);
		}	if( !is_s_binded() && is_variable_subject(var) ){
  			return m_ptr_georing->min_contained_in(c);
		}else if( !is_o_binded() && is_variable_object(var) ){
			// backwards in
			return m_ptr_georing->min_op( c);
		}

		return 0;

		// std :: cerr << "WARNING:  undefined path for in\n";

	}

	value_type leap_contained_in(var_type var, size_type c){

		if( is_s_binded() && is_variable_object(var) ){
			return m_ptr_georing->contained_in(binding_s, c);
		}else if( is_o_binded() && is_variable_subject(var) ){
			// backwards in
			return m_ptr_georing->contains(binding_o, c);
		}if( !is_s_binded() && is_variable_subject(var) ){
  			return m_ptr_georing->min_contains(c);
  		}else if( !is_o_binded() && is_variable_object(var) ){
  			// backwards in
  			return m_ptr_georing->min_contained_in( c);
  		}

		std :: cerr << "WARNING:  undefined path for in\n";
		return 0;
	}

	value_type leap_touches(var_type var, size_type c){

		if( is_s_binded() && is_variable_object(var) ){
			return m_ptr_georing->touches(binding_s, c);
		}else if( is_o_binded() && is_variable_subject(var) ){
			// backwards in
			return m_ptr_georing->touches(binding_o, c);
		}
	}
*/
	bool in_last_level() override {
		return m_level==2;
	}
/*
	std::vector<uint64_t> seek_all(var_type var) override {
    if( is_contains() ) return seek_all_contains(var);
		if( is_contained_in() ) return seek_all_contained_in(var);
		if( is_touches() ) return seek_all_touches(var);
		return std::vector<uint64_t>();
	}
*/

  std::vector<uint64_t> seek_all(var_type var) override {
    std::vector<uint64_t> res;

    //cout << "seeking all"<<endl;
    if( is_s_binded() && is_variable_object(var) ){
      //cout << "on s with binding "<<binding_s<<endl;
			value_type reg = min_op_q(binding_s);
      while( reg!=0 ){

        res.push_back(reg);
        reg = op(binding_s,reg+1);
      }
		}else if( is_o_binded() && is_variable_subject(var) ){
			// backwards in
      //cout << "on o with binding "<<binding_o<<endl;
      value_type reg = min_inv_op_q(binding_o);
      while( reg!=0 ){
        //cout << "appending "<<reg<<endl;
        res.push_back(reg);
	//cout << reg+1<<endl;
        reg = inv_op(binding_o,reg+1);
      }
      //cout << "finished " << endl;
		}else if( !is_s_binded() && is_variable_subject(var) ){
      //cout << "on s without binding "<<endl;
      value_type reg = min_op(0);
      while( reg!=0 ){
        res.push_back(reg);
        reg = min_op(reg+1);
      }
		}else if( !is_o_binded() && is_variable_object(var) ){
			// backwards in
      //cout << "on o without binding "<<endl;
      //cout <<" backwards in seek all"<<endl;
      value_type reg = min_inv_op(0);
      while( reg!=0 ){
        res.push_back(reg);
        reg = min_inv_op(reg+1);
      }
		}
		return res;
	}

  std::vector<uint64_t> seek_all_contains(var_type var)  {
    std::vector<uint64_t> res;

    if( is_s_binded() && is_variable_object(var) ){
      cout << "on s with binding "<<binding_s<<endl;
			value_type reg = m_ptr_georing->min_contains_q(binding_s);
      while( reg!=0 ){

        res.push_back(reg);
        reg = m_ptr_georing->contains(binding_s,reg+1);
      }
		}else if( is_o_binded() && is_variable_subject(var) ){
			// backwards in
      cout << "on o with binding "<<binding_o<<endl;
      value_type reg = m_ptr_georing->min_contained_in_q(binding_o);
      while( reg!=0 ){
        cout << "appending "<<reg<<endl;
        res.push_back(reg);
        reg = m_ptr_georing->contained_in(binding_o,reg+1);
      }
		}else if( !is_s_binded() && is_variable_subject(var) ){
      cout << "on s without binding "<<endl;
      value_type reg = m_ptr_georing->min_contains(0);
      while( reg!=0 ){
        res.push_back(reg);
        reg = m_ptr_georing->min_contains(reg+1);
      }
		}else if( !is_o_binded() && is_variable_object(var) ){
			// backwards in
      cout << "on o without binding "<<endl;
      cout <<" backwards in seek all"<<endl;
      value_type reg = m_ptr_georing->min_contained_in(0);
      while( reg!=0 ){
        res.push_back(reg);
        reg = m_ptr_georing->min_contained_in(reg+1);
      }
		}
		return res;
	}

  std::vector<uint64_t> seek_all_touches(var_type var)  {
    std::cout << "warning: using unimplemented seek_all_contains "<<endl;
    cout <<"on var "<<var<<endl;
    cout << "m_level is "<<m_level<<endl;
    std::vector<uint64_t> res;

    if( is_s_binded() && is_variable_object(var) ){
      cout << "on s with binding "<<binding_s<<endl;
			value_type reg = m_ptr_georing->min_touches_q(binding_s);
      while( reg!=0 ){

        res.push_back(reg);
        reg = m_ptr_georing->touches(binding_s,reg+1);
      }
		}else if( is_o_binded() && is_variable_subject(var) ){
			// backwards in
      cout << "on o with binding "<<binding_o<<endl;
      value_type reg = m_ptr_georing->min_touches_q(binding_o);
      while( reg!=0 ){
        cout << "appending "<<reg<<endl;
        res.push_back(reg);
        reg = m_ptr_georing->touches(binding_o,reg+1);
      }
		}else if( !is_s_binded() && is_variable_subject(var) ){
      cout << "on s without binding "<<endl;
      value_type reg = m_ptr_georing->min_touches(0);
      while( reg!=0 ){
        res.push_back(reg);
        reg = m_ptr_georing->min_touches(reg+1);
      }
		}else if( !is_o_binded() && is_variable_object(var) ){
			// backwards in
      cout << "on o without binding "<<endl;
      cout <<" backwards in seek all"<<endl;
      value_type reg = m_ptr_georing->min_touches(0);
      while( reg!=0 ){
        res.push_back(reg);
        reg = m_ptr_georing->min_touches(reg+1);
      }
		}
		return res;
	}
  std::vector<uint64_t> seek_all_contained_in(var_type var)  {
    std::cout << "warning: using unimplemented seek_all_contained_in "<<endl;
		return std::vector<uint64_t>();
	}

	inline var_type get_subject_variable(){
	  return m_ptr_triple_pattern.term_s.value;
	}

	inline var_type get_object_variable(){
	  return m_ptr_triple_pattern.term_o.value;
	}

	inline var_type get_predicate_variable(){
	  return m_ptr_triple_pattern.term_p.value;
	}

	inline bool is_variable_subject(var_type var) override{
            return m_ptr_triple_pattern.term_s.is_variable && var == m_ptr_triple_pattern.term_s.value;
        }

    inline bool is_variable_predicate(var_type var) override{
        return m_ptr_triple_pattern.term_p.is_variable && var == m_ptr_triple_pattern.term_p.value;
    }

    inline bool is_variable_object(var_type var) override{
        return m_ptr_triple_pattern.term_o.is_variable && var == m_ptr_triple_pattern.term_o.value;
    }




};

}

#endif //RING_LTJ_ITERATOR_GEO_HPP
