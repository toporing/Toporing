/*
 * query-index.cpp
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

#include <iostream>
#include <utility>
#include "ring.hpp"
#include <chrono>
#include <triple_pattern.hpp>
//#include <ltj_algorithm.hpp>
#include <ltj_iterator.hpp>
#include <ltj_algorithm_geo.hpp>
#include <ltj_algorithm_geo_baseline.hpp>
#include <ltj_iterator_geo.hpp>
#include "utils.hpp"
//#include <ring_muthu.hpp>
#include <time.hpp>

using namespace std;

//#include<chrono>
//#include<ctime>

using namespace ::util::time;
typedef ring::ring<> ring_type;
typedef ring::ltj_iterator<ring_type, uint8_t, uint64_t> iterator_type;

ring::triple_pattern create_dummy_triple(int s, int p) {

    ring::triple_pattern triple;
    triple.const_s(s);
    triple.const_p(75); // contains = 31
    triple.var_o(1);
    return triple;
}

std::vector<uint> possible_subjects(ring_type* ring_ptr){
    ring::triple_pattern triple;
    triple.var_s(1);
    triple.const_p(75); // contains
    triple.var_o(2);

    std::vector<uint> res;
    auto aux = iterator_type(&triple, ring_ptr);
    if(aux.is_empty()) return res;
    uint s_v = 0;
    while(true){
      s_v = aux.leap(1, s_v+1);
      if(s_v == 0) break;
      res.push_back(s_v);
    }
    return res;

}


int main(int argc, char* argv[])
{

   

    if(argc < 3){
        std::cout << "Usage: " << argv[0] << " <index> <# of queries>" << std::endl;
        return 0;
    }

    std::string index = argv[1];
    int num_queries = atoi(argv[2]);
    int num_O = 242124917;

    ring_type ring;
    sdsl::load_from_file(ring, index+"-ring");
    auto subjects = possible_subjects(&ring);

    std::cout << "Subjects: " << subjects.size() << std::endl;
    if(subjects.size() == 0) exit(0);

    int *S_ids = new int[num_queries];
    int *O_ids = new int[num_queries];
    std::vector<ring::triple_pattern> triples;
    std::vector<iterator_type> iters;

    triples.reserve(num_queries+1);
    iters.reserve(num_queries+1);
    int qs = 0;
    while(qs < num_queries) {
      S_ids[qs] = subjects[rand() % subjects.size()];
      O_ids[qs] = rand() % num_O;
      auto triple = create_dummy_triple(S_ids[qs], O_ids[qs]);
      auto aux = iterator_type(&triple, &ring);
      if(!aux.is_empty()){
         triples.push_back(triple);
         iters.push_back({&triples.back(), &ring});
         ++qs;
      }
    }
    std::cout << "Queries: " << qs << std::endl;

    ::util::time::usage::usage_type start, stop;
    uint64_t total_elapsed_time;
    uint64_t total_user_time;


    int dump = 0, x;
    // Warm-up
    for(int i=0; i < num_queries; i++){
      // set 1 as var of iterator in O
      //if(iters[i].is_empty()) continue; 
      x = iters[i].leap(1,O_ids[i]);
    }
    dump += x;

    
    cout << "query;total_elapsed_time" << endl;
    

    start = ::util::time::usage::now();

    for(int i=0; i < num_queries; i++){
      // set 1 as var of iterator in O
      //if(iters[i].is_empty()) continue;
      x = iters[i].leap(1,O_ids[i]);
    }

    stop = ::util::time::usage::now();
    dump += x;
    total_elapsed_time = (uint64_t) duration_cast<nanoseconds>(stop.elapsed - start.elapsed);
    cout << "leap_ring;" << total_elapsed_time << endl;

	return 0;
}
