/*
 * build-index.cpp
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
#include "ring.hpp"
//#include "ring_muthu.hpp"
#include <fstream>
#include <sdsl/construct.hpp>
//#include <ltj_algorithm.hpp>
//#include <georing.hpp>
#include <toporing.hpp>
using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

int IN_ID = -1;
int CONTAINED_ID = 0;
int TOUCHES_ID = 0;




unsigned try_map(unsigned u, unordered_map<unsigned,unsigned> &mapping, 
	unordered_map<unsigned,unsigned> &mapping_inv){

  if(mapping.count(u)!=0){
    return mapping[u];
  }else{
    unsigned new_id = mapping.size() + 1;
    mapping[u] = new_id;
    //mapping_inv[new_id] = u;
    //cout <<"P "<< u << " -> " << new_id << endl;
    return new_id;
  }


}
template<class georing_type, class ring_type>
void build_georing(const std::string &dataset, const std::string &output){
    unordered_map<unsigned, vector<unsigned > > hierarchy;
    vector<pair<unsigned,unsigned> > adj;
    vector<spo_triple> C, T,R;


    std::ifstream ifs(dataset);

    cerr << "reading tuples"<<endl;
    uint64_t s, p , o, maxId=0;
    
    while(ifs >> s >> p >> o){
        //skip non hierarchy relations
        //if(p!=TOPO_IN_ID) continue;

        if(p==IN_ID){
	    if(s==o) continue;
	    //cout << s <<"-"<<IN_ID<<"-"<<o<<endl;
            C.push_back(spo_triple(s, p, o));
            maxId = max({maxId,s,o});
        }else if(p==CONTAINED_ID){
	    if(s==o) continue;
	    //cout << s <<"-"<<CONTAINED_ID<<"-"<<o<<endl;
            C.push_back(spo_triple(o, p, s));
            maxId = max({maxId,s,o});
        }else if(p==TOUCHES_ID){
	    if(s==o) continue;
            T.push_back(spo_triple(s, p, o));
            //T.push_back(spo_triple(o, p, s));
            maxId = max({maxId,s-1,o-1});
        }else{
          R.push_back(spo_triple(s, p, o));
        }

    }
    unordered_set<unsigned> cant_be_root;
    cerr << "finished reading tuples" <<endl;
    for(spo_triple t: C){
        s = get<0>(t);
        o = get<2>(t);
	//can_be_root[o]=false;
	cant_be_root.insert(o);
        hierarchy[s].push_back(o);
	//cerr << s << " contains "<<o <<endl;
	if(hierarchy.count(o) == 0) hierarchy[o] = {}; 
    }

    for(spo_triple t: T){
        s = get<0>(t);
        o = get<2>(t);
        adj.push_back({s,o});
    }
    unordered_map<unsigned, unsigned> mapping;
    unordered_map<unsigned, unsigned> mapping_inv;

    //int root = -1;
    vector<unsigned> roots;
    /*
    for(int i=0; i< can_be_root.size(); i++){
      if(can_be_root[i]) root = i;
    }
    */

    for(auto &p: hierarchy){
      unsigned u = p.first;
      if( cant_be_root.count(u) == 0 ) roots.push_back(u);
    }
    if(roots.empty() ) cout << "error: no possible root found"<< endl;
    cerr << "building georing"<<endl;
    georing_type G(hierarchy, adj, roots, mapping, mapping_inv);

    cout << "mapping ring tuples"<<endl;
    //TODO: complete mapping for arbitrary data
/*
    for(int i=0; i<mapping.size(); i++){
      cout <<i<<" has been mapped to "<<mapping[i]<<endl;
    }
    cout << "num of vertices "<<G.num_vertices()<<endl;
    */
    for(int i=0; i<R.size(); i++){
        auto t = R[i];
        s = get<0>(t);
        p = get<1>(t);
        o = get<2>(t);
	//cout <<"B "<< s <<" "<<p<<" "<<o<<endl;

	s = try_map(s, mapping, mapping_inv);
	o = try_map(o, mapping, mapping_inv);

	//cout <<"A "<< s <<" "<<p<<" "<<o<<endl;
	R[i] = spo_triple({s,p,o});
    }


    cout <<"storing permutation"<<endl;

    sdsl::int_vector<> iv(mapping.size()+1);
    iv[0]=0;

    for(auto p: mapping){
	    if(p.first > iv.size()) iv.resize(p.first+1);
	    iv[p.first] = p.second;
    }

    sdsl::store_to_file(iv, output+"-perm");



    cout <<"building ring"<<endl;
    ring_type A(R);
    cout <<"ring construction finished, storing..."<<endl;
    sdsl::store_to_file(G, output);
    sdsl::store_to_file(A, output+"-ring");
}

template<class ring>
void build_index(const std::string &dataset, const std::string &output){
    vector<spo_triple> D, E;

    std::ifstream ifs(dataset);
    uint64_t s, p , o;
    do {
        ifs >> s >> p >> o;
        D.push_back(spo_triple(s, p, o));
    } while (!ifs.eof());

    D.shrink_to_fit();
    cout << "--Indexing " << D.size() << " triples" << endl;
    memory_monitor::start();
    auto start = timer::now();

    ring A(D);
    auto stop = timer::now();
    memory_monitor::stop();
    cout << "  Index built  " << sdsl::size_in_bytes(A) << " bytes" << endl;

    sdsl::store_to_file(A, output);
    cout << "Index saved" << endl;
    cout << duration_cast<seconds>(stop-start).count() << " seconds." << endl;
    cout << memory_monitor::peak() << " bytes." << endl;

}

int main(int argc, char **argv)
{

    if(argc < 5){
        std::cout << "Usage: " << argv[0] << " <dataset> IN_ID CONTAINED_ID TOUCHES_ID" << std::endl;
        return 0;
    }

    std::string dataset = argv[1];
    //std::string type    = argv[2];
	IN_ID = atoi(argv[2]);
	CONTAINED_ID = atoi(argv[3]);
	TOUCHES_ID = atoi(argv[4]);
        std::string index_name = dataset + ".georing";
        build_georing<ring::toporing<>, ring::ring<>>(dataset, index_name);
   

    return 0;
}
