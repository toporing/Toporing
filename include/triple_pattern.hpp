/*
 * triple_pattern.hpp
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



#ifndef RING_TRIPLE_PATTERN_HPP
#define RING_TRIPLE_PATTERN_HPP

#include "configuration.hpp"

namespace ring {



    struct term_pattern {
        uint64_t value; //TODO: transform char of variable to uint64_t
        bool is_variable;
    };

    struct triple_pattern {
        term_pattern term_s;
        term_pattern term_p;
        term_pattern term_o;
	bool is_geo_pattern=false;



	void set_geo(string &s){

	  is_geo_pattern = true;
	  term_p.is_variable = false;
	  // TODO: maybe hash these?
	  if(s=="CONTAINS"){
	    term_p.value = TOPO_CONTAINS;
	  }else if(s=="CONTAINED_IN"){
	    term_p.value = TOPO_CONTAINED_IN;
	  }else if(s=="TOUCHES"){
	    term_p.value = TOPO_TOUCHES;
	  }else if(s=="NOT_CONTAINS"){
	    term_p.value = TOPO_NOT_CONTAINS;
	  }else if(s=="NOT_CONTAINED_IN"){
	    term_p.value = TOPO_NOT_CONTAINED_IN;
	  }else if(s=="NOT_TOUCHES"){
			term_p.value = TOPO_NOT_TOUCHES;
		}

	}

	std::string topo_name() const{
	  switch(term_p.value){
	    case TOPO_CONTAINS:
	      return "TOPO_CONTAINS";
	    case TOPO_CONTAINED_IN:
	      return "TOPO_CONTAINED_IN";
	    case TOPO_TOUCHES:
	      return "TOPO_TOUCHES";
	    case TOPO_NOT_CONTAINS:
	      return "TOPO_NOT_CONTAINS";
	    case TOPO_NOT_CONTAINED_IN:
	      return "TOPO_NOT_CONTAINED_IN";
			case TOPO_NOT_TOUCHES:
				return "TOPO_NOT_TOUCHES";
	    default:
	      return "ERROR";
	  }
	}

        void const_s(uint64_t s){
            term_s.is_variable = false;
            term_s.value = s;
        }

        void const_o(uint64_t o){
            term_o.is_variable = false;
            term_o.value = o;
        }

        void const_p(uint64_t p){
            term_p.is_variable = false;
            term_p.value = p;
        }

        void var_s(uint64_t s){
            term_s.is_variable = true;
            term_s.value = s;
        }

        void var_o(uint64_t o){
            term_o.is_variable = true;
            term_o.value = o;
        }

        void var_p(uint64_t p){
            term_p.is_variable = true;
            term_p.value = p;
        }

        bool s_is_variable() const {
            return term_s.is_variable;
        }

        bool p_is_variable() const {
            return term_p.is_variable;
        }

        bool o_is_variable() const {
            return term_o.is_variable;
        }

	bool is_geo() const {
	  return is_geo_pattern;
	}


        void print(std::unordered_map<uint8_t, std::string> &ht) const {
            if(s_is_variable()){
                std::cout << "?" << ht[term_s.value] << " ";
            }else{
                std::cout << term_s.value << " ";
            }

	    if(is_geo()){
                std::cout << topo_name() << " ";
	    }else if(p_is_variable()){
                std::cout << "?" << ht[term_p.value] << " ";
            }else{
                std::cout << term_p.value << " ";
            }

            if(o_is_variable()){
                std::cout << "?" << ht[term_o.value];
            }else{
                std::cout << term_o.value;
            }
        }
    };
}


#endif //RING_TRIPLE_PATTERN_HPP
