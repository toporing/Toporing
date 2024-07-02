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

ring::triple_pattern create_dummy_triple(int s, int p) {

    ring::triple_pattern triple;
    triple.const_s(s);
    triple.const_p(31); // contains
    triple.var_o(1);
}


int main(int argc, char* argv[])
{

    typedef ring::ring<> ring_type;

    if(argc < 3){
        std::cout << "Usage: " << argv[0] << " <index> <# of queries>" << std::endl;
        return 0;
    }

    std::string index = argv[1];
    int num_queries = atoi(argv[2]);
    int num_S = 107836911;
    int num_O = 242124917;

    ring_type ring;
    sdsl::load_from_file(ring, index+"-ring");

    typedef ring::ltj_iterator<ring_type, uint8_t, uint64_t> iterator_type;
    int *S_ids = new int[num_queries];
    int *O_ids = new int[num_queries];
    std::vector<ring::triple_pattern> triples;
    std::vector<iterator_type> iters;

    triples.reserve(num_queries+1);
    iters.reserve(num_queries+1);
    for(int i=0; i < num_queries; i++) {
      S_ids[i] = rand() % num_S;
      O_ids[i] = rand() % num_O;
      triples.push_back(create_dummy_triple(S_ids[i], O_ids[i]));
      iters.push_back({&triples.back(), &ring});
    }

    ::util::time::usage::usage_type start, stop;
    uint64_t total_elapsed_time;
    uint64_t total_user_time;


    int dump = 0, x;
    // Warm-up
    for(int i=0; i < num_queries; i++){
      // set 1 as var of iterator in O 
      x = iters[i].leap(1,O_ids[i]);
    }
    dump += x;

    
    cout << "query;total_elapsed_time" << endl;
    

    start = ::util::time::usage::now();

    for(int i=0; i < num_queries; i++){
      // set 1 as var of iterator in O
      x = iters[i].leap(1,O_ids[i]);
    }

    stop = ::util::time::usage::now();
    dump += x;
    total_elapsed_time = (uint64_t) duration_cast<nanoseconds>(stop.elapsed - start.elapsed);
    cout << "leap_ring;" << total_elapsed_time << endl;

	return 0;
}
