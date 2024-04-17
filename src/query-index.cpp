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

bool get_file_content(string filename, vector<string> & vector_of_strings)
{
    // Open the File
    ifstream in(filename.c_str());
    // Check if object is valid
    if(!in)
    {
        cerr << "Cannot open the File : " << filename << endl;
        return false;
    }
    string str;
    // Read the next line from File until it reaches the end.
    while (getline(in, str))
    {
        // Line contains string of length > 0 then save it in vector
        if(str.size() > 0)
            vector_of_strings.push_back(str);
    }
    //Close The File
    in.close();
    return true;
}

std::string ltrim(const std::string &s)
{
    size_t start = s.find_first_not_of(' ');
    return (start == std::string::npos) ? "" : s.substr(start);
}

std::string rtrim(const std::string &s)
{
    size_t end = s.find_last_not_of(' ');
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

std::string trim(const std::string &s) {
    return rtrim(ltrim(s));
}

std::vector<std::string> tokenizer(const std::string &input, const char &delimiter){
    std::stringstream stream(input);
    std::string token;
    std::vector<std::string> res;
    while(getline(stream, token, delimiter)){
        res.emplace_back(trim(token));
    }
    return res;
}

bool is_variable(string & s)
{
    return (s.at(0) == '?');
}

uint8_t get_variable(string &s, std::unordered_map<std::string, uint8_t> &hash_table_vars){
    auto var = s.substr(1);
    auto it = hash_table_vars.find(var);
    if(it == hash_table_vars.end()){
        uint8_t id = hash_table_vars.size();
        hash_table_vars.insert({var, id });
        return id;
    }else{
        return it->second;
    }
}

bool is_geo(string &s){

  return s=="CONTAINS" or s=="CONTAINED_IN" or s=="TOUCHES" or s=="NOT_CONTAINS" or s=="NOT_CONTAINED_IN";

}

uint64_t get_constant(string &s){
    return std::stoull(s);
}

ring::triple_pattern get_triple(string & s, std::unordered_map<std::string, uint8_t> &hash_table_vars) {
    vector<string> terms = tokenizer(s, ' ');

    ring::triple_pattern triple;
    if(is_geo(terms[1])){
      triple.is_geo_pattern=true;
    }
    if(is_variable(terms[0])){
        triple.var_s(get_variable(terms[0], hash_table_vars));
    }else{
        triple.const_s(get_constant(terms[0]));
    }
    if(is_geo(terms[1])){
      triple.set_geo(terms[1]);
    }else if(is_variable(terms[1])){
        triple.var_p(get_variable(terms[1], hash_table_vars));
    }else{
        triple.const_p(get_constant(terms[1]));
    }
    if(is_variable(terms[2])){
        triple.var_o(get_variable(terms[2], hash_table_vars));
    }else{
        triple.const_o(get_constant(terms[2]));
    }
    return triple;
}

std::string get_type(const std::string &file){
    auto p = file.find_last_of('.');
    return file.substr(p+1);
}


void permute_triple(sdsl::int_vector<> &perm, ring::triple_pattern &triple){

	if( !triple.term_s.is_variable ){
		triple.term_s.value = perm[triple.term_s.value]; 
	}

	if( !triple.term_o.is_variable ){
		triple.term_o.value = perm[triple.term_o.value]; 
	}

	//cout << "PERM -> "<< triple.term_s.value << " " << triple.term_p.value << " " <<triple.term_o.value <<endl;

}

template<class ring_type, class georing_type, class trait_type>
void queryGeo(const std::string &file, const std::string &queries, const uint64_t limit){
    vector<string> dummy_queries;
    bool result = get_file_content(queries, dummy_queries);

    georing_type graph;

    //cout << " Loading the index..."; fflush(stdout);
    sdsl::load_from_file(graph, file);
    ring_type ring;
    sdsl::load_from_file(ring, file+"-ring");

    sdsl::int_vector<> perm;
    sdsl::load_from_file(perm, file+"-perm");

    sdsl::int_vector<> inv_perm(perm.size());
    for(int i=0; i<perm.size(); i++) inv_perm[perm[i]] = i;

    cerr << "Total toporing size: "<< size_in_mega_bytes(ring) + size_in_mega_bytes(graph) + size_in_mega_bytes(perm) + size_in_mega_bytes(inv_perm) <<"MB" << endl;
    //cout << endl << " Index loaded " << sdsl::size_in_bytes(graph) << " bytes" << endl;

    cerr << "ring;toporing"<<endl;
    cerr << size_in_mega_bytes(ring) << ";"<< size_in_mega_bytes(graph) <<endl;
    ofstream ff("report.html");
    sdsl::write_structure<HTML_FORMAT>(graph, ff);
    std::ifstream ifs;
    uint64_t nQ = 0;

    ::util::time::usage::usage_type start, stop;
    uint64_t total_elapsed_time;
    uint64_t total_user_time;

    if(result)
    {

        int count = 1;
        for (string& query_string : dummy_queries) {

            //vector<Term*> terms_created;
            //vector<Triple*> query;
            std::unordered_map<std::string, uint8_t> hash_table_vars;
            std::vector<ring::triple_pattern> query;
            vector<string> tokens_query = tokenizer(query_string, '.');
            for (string& token : tokens_query) {
                auto triple_pattern = get_triple(token, hash_table_vars);
		permute_triple(perm, triple_pattern);
                query.push_back(triple_pattern);
            }


            typedef ring::ltj_iterator<ring_type, uint8_t, uint64_t> iterator_type;
            typedef ring::ltj_iterator_geo<georing_type, uint8_t, uint64_t> geoiterator_type;
            typedef ring::ltj_algorithm_geo<iterator_type,geoiterator_type,
                                 ring::veo::veo_adaptive<iterator_type, trait_type>> algorithm_type;
            typedef std::vector<typename algorithm_type::tuple_type> results_type;
            results_type res;

            start = ::util::time::usage::now();
            algorithm_type ltj(&query, &ring, &graph);

            std::unordered_map<uint8_t, std::string> ht;
            for(const auto &p : hash_table_vars){
                ht.insert({p.second, p.first});
            }
            //ltj.print_veo(ht);


            ltj.join(res, limit, 600);
            stop = ::util::time::usage::now();

            total_elapsed_time = (uint64_t) duration_cast<nanoseconds>(stop.elapsed - start.elapsed);
            total_user_time = (uint64_t) duration_cast<nanoseconds>(stop.user - start.user);


            //cout << "Query Details:" << endl;
	    //cout << query_string << endl;
            //ltj.print_query(ht, inv_perm);
            //cout << "##########" << endl;
            //ltj.print_results(res, ht, inv_perm);
            cout << nQ <<  ";" << res.size() << ";" << total_elapsed_time << ";" << total_user_time << endl;
            nQ++;
          }
        }

}




template<class ring_type, class georing_type, class trait_type>
void queryGeoBaseline(const std::string &file, const std::string &queries, const uint64_t limit){
    vector<string> dummy_queries;
    bool result = get_file_content(queries, dummy_queries);

    georing_type graph;

    //cout << " Loading the index..."; fflush(stdout);
    sdsl::load_from_file(graph, file);
    ring_type ring;
    sdsl::load_from_file(ring, file+"-ring");
    //cout << endl << " Index loaded " << sdsl::size_in_bytes(graph) << " bytes" << endl;
    sdsl::int_vector<> perm;
    sdsl::load_from_file(perm, file+"-perm");


    sdsl::int_vector<> inv_perm(perm.size());
    for(int i=0; i<perm.size(); i++) inv_perm[perm[i]] = i;

    std::ifstream ifs;
    uint64_t nQ = 0;

    ::util::time::usage::usage_type start, stop;
    uint64_t total_elapsed_time;
    uint64_t total_user_time;

    if(result)
    {

        int count = 1;
        for (string& query_string : dummy_queries) {

            //vector<Term*> terms_created;
            //vector<Triple*> query;
            std::unordered_map<std::string, uint8_t> hash_table_vars;
            std::vector<ring::triple_pattern> query;
            vector<string> tokens_query = tokenizer(query_string, '.');
            for (string& token : tokens_query) {
                auto triple_pattern = get_triple(token, hash_table_vars);
		permute_triple(perm, triple_pattern);
                query.push_back(triple_pattern);
            }


            typedef ring::ltj_iterator<ring_type, uint8_t, uint64_t> iterator_type;
            typedef ring::ltj_iterator_geo<georing_type, uint8_t, uint64_t> geoiterator_type;
            typedef ring::ltj_algorithm_geo_baseline<iterator_type,geoiterator_type,
                                 ring::veo::veo_adaptive<iterator_type, trait_type>> algorithm_type;
// #if ADAPTIVE
//
//             typedef ring::ltj_algorithm<iterator_type,
//                     ring::veo::veo_adaptive<iterator_type, trait_type>> algorithm_type;
//
// #else
// #if VEO_RANDOM
//             typedef ring::ltj_algorithm<iterator_type,
//                     ring::veo::veo_simple_random<iterator_type, trait_type>> algorithm_type;
// #else
//             typedef ring::ltj_algorithm<iterator_type,
//                     ring::veo::veo_simple<iterator_type, trait_type>> algorithm_type;
// #endif
// #endif
            typedef std::vector<typename algorithm_type::tuple_type> results_type;
            results_type res;

            start = ::util::time::usage::now();
            algorithm_type ltj(&query, &ring, &graph);
            ltj.join(res, limit, 600);
            stop = ::util::time::usage::now();

            total_elapsed_time = (uint64_t) duration_cast<nanoseconds>(stop.elapsed - start.elapsed);
            total_user_time = (uint64_t) duration_cast<nanoseconds>(stop.user - start.user);

            std::unordered_map<uint8_t, std::string> ht;
            for(const auto &p : hash_table_vars){
                ht.insert({p.second, p.first});
            }

            //cout << "Query Details:" << endl;
	    
	    //cout << query_string << endl;
            //ltj.print_query(ht, inv_perm);
            //ltj.print_veo(ht);
            //cout << "##########" << endl;
            //ltj.print_results(res, ht, inv_perm);
            cout << nQ <<  ";" << res.size() << ";" << total_elapsed_time << ";" << total_user_time << endl;
            nQ++;
          }
        }

}

int main(int argc, char* argv[])
{

    typedef ring::ring<> ring_type;
    //typedef ring::c_ring ring_type;
    //
    if(argc < 4){
        std::cout << "Usage: " << argv[0] << " <index> <queries> <limit> [baseline]" << std::endl;
        return 0;
    }

    std::string index = argv[1];
    std::string queries = argv[2];
    uint64_t limit = std::atoll(argv[3]);
    int is_baseline = 0;
    if(argc>=5) is_baseline =  atoi(argv[4]);
    std::string type = get_type(index);

    if (is_baseline==0 and type == "georing"){
        queryGeo<ring::ring<>, ring::toporing<>,  ring::util::trait_size>(index, queries, limit);
    }else if (is_baseline==1 and type == "georing"){
        queryGeoBaseline<ring::ring<>, ring::toporing<>,  ring::util::trait_size>(index, queries, limit);
    }else{
        std::cout << "Type of index: " << type << " is not supported." << std::endl;
    }



	return 0;
}
