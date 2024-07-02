#ifndef GEORING_GENERAL_HPP
#define GEORING_GENERAL_HPP

#include "georing_base.hpp" 
#include <sdsl/bp_support_sada.hpp>
#include "configuration.hpp"
#include <omp.h>

using namespace std;

namespace ring{
  template < class wt_t = wm_int<>,
	   class bp_t = bp_support_sada<>,
	   class bit_vector_t = bit_vector,
	   class rank_t = typename bit_vector_t::rank_1_type,
	   class select1_t = typename bit_vector_t::select_1_type ,
	   class select0_t = typename bit_vector_t::select_0_type >
  class toporing : public georing_base{
    private:
      // Topology of the hierarchy
      bit_vector_t m_T_bv;
      bp_t m_T_bp;
      // Touches relation stored for leaves, as matrix
      // stored with WT and BV
      wt_t m_M_wt;
      bit_vector_t m_M_bv;
      select1_t m_M_bv_select1;
      select0_t m_M_bv_select0;
      rank_t m_M_bv_rank;
      size_type m_num_leaves;
      select0_t m_T_select0;

      bit_vector_t m_contains;
      select1_t m_contains_select1;
      select0_t m_contains_select0;
      rank_t m_contains_rank;

      bit_vector_t m_contained;
      select1_t m_contained_select1;
      select0_t m_contained_select0;
      rank_t m_contained_rank;

      bit_vector_t m_not_contains;
      select1_t m_not_contains_select1;
      select0_t m_not_contains_select0;
      rank_t m_not_contains_rank;

      bit_vector_t m_not_contained;
      select1_t m_not_contained_select1;
      select0_t m_not_contained_select0;
      rank_t m_not_contained_rank;

      bit_vector_t m_touches;
      select1_t m_touches_select1;
      select0_t m_touches_select0;
      rank_t m_touches_rank;

      bit_vector_t m_not_touches;
      select1_t m_not_touches_select1;
      select0_t m_not_touches_select0;
      rank_t m_not_touches_rank;

      bit_vector_t m_is_leaf;
      select1_t m_is_leaf_select1;
      rank_t m_is_leaf_rank;

      int_vector<> m_run_inits;
      int_vector<> m_run_endings;
      bit_vector_t m_run_limits;
      select1_t m_run_limits_select;
      rank_t m_run_limits_rank;

      size_type total_leaves;
      
      size_type m_threshold;

    public:

      size_type num_vertices() const{
	return m_T_bv.size()/2;
      }
      
      toporing() = default;

      toporing( unordered_map<unsigned, vector <unsigned> > &hierarchy, const vector<pair<unsigned, unsigned> > &adjacency, vector<unsigned> &roots, unordered_map<unsigned, unsigned> &mapping, unordered_map<unsigned, unsigned> &mapping_inv){
	
	size_type vertices = hierarchy.size();



	bit_vector_t t_local(2*vertices,0);
	//int_vector<> c_local (vertices + 1);
	stack < value_type > s;
	unordered_map<unsigned,unsigned> visited;//(vertices+2, false);
	

	bit_vector_t contains_local(vertices+1,0);
	bit_vector_t contained_local(vertices+1,0);
	bit_vector_t touches_local(vertices+1,0);
	bit_vector_t not_touches_local(vertices+1,0);
	for(unsigned root: roots){
		cerr << root <<" is one root"<<endl;
		s.push(root);
	}
	bit_vector_t not_contains_local(vertices+1,0);
	bit_vector_t not_contained_local(vertices+1,0);

	size_type t_idx = 0;
	size_type c_idx = 0;

	bit_vector is_leaf_local(vertices+1,0);

	total_leaves = 0;

	unsigned curr_postorder = 0;
	unsigned excess = 0;
	while ( !s.empty() ){
	  unsigned u = s.top();
	  //if(!visited[u]){
	  if(visited[u]==0){
	    excess++;
	    //first time visiting, write ( and push children
	    t_local[t_idx++] = 1;
	    // HERE: write mapping to old id
	    //for( unsigned v: hierarchy[u]){
	    for(auto v=hierarchy[u].rbegin(); v!=hierarchy[u].rend(); v++){
	      s.push(*v);
	    }
	    
	    visited[u] = 1;
	    //visited.insert(u);
	  }else if(visited[u]==1){
	    excess--;
	    // second time visiting, remove from stack and write )
	    mapping[u] = ++curr_postorder;
	    mapping_inv[curr_postorder] = u;
	    //cout <<"P "<< u << " -> " << mapping[u] <<endl;
	    s.pop();
	    t_local[t_idx++] = 0;
	    if(hierarchy[u].size()==0){
		    is_leaf_local[curr_postorder] = 1;
		    total_leaves++;
	    }
	    visited[u]=2;
	    /*
	    if(hierarchy[u].size()!=0){
	      //u contains something
	      contains_local[curr_postorder]=1;
	    }
	    // puede estar tocando solo entre descendientes, revisar!
	    if(excess!=0){
	      //u is contained by something
	      contained_local[curr_postorder]=1;
	    }
	    */
	  }
	  
	  else{
		  s.pop();
		  cerr << "Error: node "<<u <<" has already been explored" <<endl;
	  }
	  
	}

	if(visited.size() != hierarchy.size()){
		cout << "Error: hierarchy was not completely traversed" <<endl;
		for(auto &p : hierarchy){
			unsigned id = p.first;
			if(visited.count(id)==0){
				cerr << id << "was not visited!"<<endl;
				for(auto x : p.second ){
					cerr<< "    "<<x <<endl;
				}
			}
		}
	}

	  //cout << "total leaves is "<<total_leaves <<endl;
	m_is_leaf.swap(is_leaf_local);
	sdsl::util::init_support(m_is_leaf_rank, &m_is_leaf);
	sdsl::util::init_support(m_is_leaf_select1, &m_is_leaf);

	m_T_bv.swap(t_local);
	bp_t t_bp_local(&m_T_bv);
	m_T_bp.swap(t_bp_local);
	m_T_bp.set_vector(&m_T_bv);
	sdsl::util::init_support(m_T_select0, &m_T_bv);
	unsigned m = adjacency.size();
	bit_vector b_local(num_vertices() + m+1);


	vector<pair<value_type,value_type> > adjMapped;

	for(unsigned int i=0; i<m; i++){

		if(mapping.count(adjacency[i].first)==0 or mapping.count(adjacency[i].second)==0){
			cerr << "adjacency "<<adjacency[i].first <<" " << adjacency[i].second <<
				"can not be mapped (not in hierarchy)" << endl;
			continue;
		}
		adjMapped.push_back({mapping[adjacency[i].first],mapping[adjacency[i].second]});
		//cout << adjMapped[i].first << " " << adjMapped[i].second <<endl;
		if(adjMapped.back().first > adjMapped.back().second) adjMapped.pop_back();

	}
	sort(adjMapped.begin(),adjMapped.end());

	m = adjMapped.size();
	value_type idxA = 0;
	value_type idxB = 0;
	value_type idxMWT = 0;
	int_vector<> m_wt_local(m+1);

	for(pair<value_type,value_type> p: adjMapped){
		//if(p.first > p.second) continue;
		// continuing with next row
		while(idxA< p.first){
			idxA++;
			b_local[idxB++] = 1;
		}

		b_local[idxB++] = 0;
		m_wt_local[idxMWT++] = p.second;

	}
	while(idxB < b_local.size()) b_local[idxB++]=1;
	//b_local[idxB++]=1;
        m_wt_local[idxMWT++] = 0;
	swap(b_local, m_M_bv);
	construct_im(m_M_wt,m_wt_local);
	sdsl::util::init_support(m_M_bv_select1, &m_M_bv);
	sdsl::util::init_support(m_M_bv_select0, &m_M_bv);
	sdsl::util::init_support(m_M_bv_rank, &m_M_bv);
	// TODO: add adjacency here




	// NOT TOUCHES DS
	//
	

	

	int t = 2;
	m_threshold = t;
	//vector<vector<int> >  vx(num_vertices()+1);
  
	//pos pares para inicio run, impares para fin run
	
	int total_runs = 0;

	vector<vector<int> > runs(num_vertices(), vector<int>());
#ifdef PARALLEL_CONSTRUCTION
#pragma omp parallel for
#endif
	for(int i=1; i<=num_vertices(); i++){

	  auto &this_runs = runs[i-1];
	  //processing leaves here->
	  //


	  //cout << "rl(i) -> "; 


	  int run_start = 0;
	  int curr_run_size = 0;

	  //for(int j=postorder_to_leaf_id(touches(i,1)); j<total_leaves; j++){
	  /*
	  for(int j=1; j<total_leaves; j++){

		  if(j==0) break;
		  int v = leaf_id_to_postorder(j);
		  // i touches v
		  if(touches(i,v) == v){
	//	    cout << "leaf id " << j << " postorder "<<v << " touches " << i <<endl;
			  if(run_start==0){
			    
			      value_type t_leaf = leaf_id_to_postorder(j+t-1);
				  if(touches(i,t_leaf)!=t_leaf){
				    //can skip, as this indicates no run
				    
				    value_type next_postorder = touches(i, t_leaf);
				    value_type leaf_id = postorder_to_leaf_id(next_postorder);
				    j = leaf_id - 1;
				    curr_run_size = 0;
				    continue;

				  }
				  
			  
				  run_start = v;
				  curr_run_size = 0;
			  }
			  curr_run_size++;

		  }else if(run_start!=0 && curr_run_size>=t){
				
			  cout << " " << curr_run_size ;
			  
			  curr_run_size = 0;
		  }else{

			  curr_run_size = 0;
		  }


	  }

	  cout << endl;
	  */
	  // get first id that touches i

	  //cout << "rn(i) -> ";
	  //
	  
	  int a = touches(i, 1);  
	  int lb = a;

	  while(a != 0){
	    //verify end
	    if ( touches(i, a+t-1) != a+t-1 ){
	      a = touches(i, a+t-1);
	      lb = a;
	      continue;
	    }

	    bool complete = true;
	    for( int j = a+t-1; j>lb; j--){
	      // try backwards first
	      if ( touches(i, j) != j ){
		// this node cuts the run, try starting new one from next
		a = touches(i, j);
		// we know that every node in [a,lb=a+t] touches i
		// we dont have to check those in next steps
		lb = a+t;
		complete = false;
		break;

	      }
	    }

	    if(complete){
	      // we know that [a, a+t-1] is a valid run, try to increase
	      // interval to the right
	      int init = a;
	      int ending = a+t-1;
	      //vx[i].push_back(a);
	      //vx[i].push_back(a+t);
	      int j = a+t;
	      while ( touches(i, j)==j ){
		ending = j++;
	      }

	      this_runs.push_back(init);
	      this_runs.push_back(ending);
	      //cout << " "<< ending - init;
	      a = lb = touches(i,j);

	      
	      //total_runs++;
	    }

	  
	  }
	  //cout <<endl;
	  //total_runs += vx[i].size()/2;
	}

	for(auto &v: runs){
		total_runs += v.size()/2;
	}
	cout << "total runs with every node -> "<<total_runs << endl;

	m_run_inits.resize(total_runs);
	m_run_endings.resize(total_runs);

	bit_vector_t run_limits_local(total_runs + num_vertices() + 1);
	int run_idx = 0;
	int bv_idx = 0;

	for(auto &v_runs : runs){
	  run_limits_local[bv_idx++] = 1;
	  for(int i=0; i<v_runs.size(); i+=2){
	    run_limits_local[bv_idx++] = 0;
	    int run_init = v_runs[i];
	    int run_ending = v_runs[i+1];

	    m_run_inits[run_idx] = run_init;
	    m_run_endings[run_idx] = run_ending;

	    run_idx++;
	  }

	}

	run_limits_local[bv_idx++] = 1;

	m_run_limits.swap(run_limits_local);
	sdsl::util::init_support(m_run_limits_select, &m_run_limits);
	sdsl::util::init_support(m_run_limits_rank, &m_run_limits);
	
	//debug_print_runs();

	
	for(int i=1; i<=num_vertices(); i++){
	  if(touches(i,1) !=0){
	    touches_local[i]=1;
	  }
	}
	
	for(int i=1; i<=num_vertices(); i++){
	  if(not_touches(i,1) !=0){
	    not_touches_local[i]=1;
	  }
	}
	

	for(int i=1; i<=num_vertices(); i++){
	  if(not_contains(i,1) !=0){
	    not_contains_local[i]=1;
	  }
	}
	
	for(int i=1; i<=num_vertices(); i++){
	  if(not_contained_in(i,1) !=0){
	    not_contained_local[i]=1;
	  }
	}

	
	for(int i=1; i<=num_vertices(); i++){
	  if(contains(i,1) !=0){
	    contains_local[i]=1;
	  }
	}
	
	for(int i=1; i<=num_vertices(); i++){
	  if(contained_in(i,1) !=0){
	    contained_local[i]=1;
	  }
	}
	m_contains.swap(contains_local);
	sdsl::util::init_support(m_contains_select0, &m_contains);
	sdsl::util::init_support(m_contains_select1, &m_contains);
	sdsl::util::init_support(m_contains_rank, &m_contains);

	m_contained.swap(contained_local);
	sdsl::util::init_support(m_contained_select0, &m_contained);
	sdsl::util::init_support(m_contained_select1, &m_contained);
	sdsl::util::init_support(m_contained_rank, &m_contained);

	
	m_touches.swap(touches_local);
	sdsl::util::init_support(m_touches_select0, &m_touches);
	sdsl::util::init_support(m_touches_select1, &m_touches);
	sdsl::util::init_support(m_touches_rank, &m_touches);

	m_not_touches.swap(not_touches_local);
	sdsl::util::init_support(m_not_touches_select0, &m_not_touches);
	sdsl::util::init_support(m_not_touches_select1, &m_not_touches);
	sdsl::util::init_support(m_not_touches_rank, &m_not_touches);

	m_not_contains.swap(not_contains_local);
	sdsl::util::init_support(m_not_contains_select0, &m_not_contains);
	sdsl::util::init_support(m_not_contains_select1, &m_not_contains);
	sdsl::util::init_support(m_not_contains_rank, &m_not_contains);

	m_not_contained.swap(not_contained_local);
	sdsl::util::init_support(m_not_contained_select0, &m_not_contained);
	sdsl::util::init_support(m_not_contained_select1, &m_not_contained);
	sdsl::util::init_support(m_not_contained_rank, &m_not_contained);
      }
		size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_T_bv.serialize(out, child, "T");
            written_bytes += m_T_bp.serialize(out, child, "T_bp");
            written_bytes += m_T_select0.serialize(out, child, "T_select0");
            written_bytes += m_M_bv.serialize(out, child, "M_bv");
            written_bytes += m_M_wt.serialize(out, child, "M_wt");
            written_bytes += m_M_bv_select1.serialize(out, child, "M_bv_select1");
            written_bytes += m_M_bv_select0.serialize(out, child, "M_bv_select0");
            written_bytes += m_M_bv_rank.serialize(out, child, "M_bv_rank");



            written_bytes += m_contains.serialize(out, child, "contains");
            written_bytes += m_contains_select1.serialize(out, child, "contains_select1");
            written_bytes += m_contains_select0.serialize(out, child, "contains_select0");
            written_bytes += m_contains_rank.serialize(out, child, "contains_rank");

            written_bytes += m_contained.serialize(out, child, "contained");
            written_bytes += m_contained_select1.serialize(out, child, "contained_select1");
            written_bytes += m_contained_select0.serialize(out, child, "contained_select0");
            written_bytes += m_contained_rank.serialize(out, child, "contained_rank");

            written_bytes += m_touches.serialize(out, child, "touches");
            written_bytes += m_touches_select1.serialize(out, child, "touches_select1");
            written_bytes += m_touches_select0.serialize(out, child, "touches_select0");
            written_bytes += m_touches_rank.serialize(out, child, "touches_rank");
	    
            written_bytes += m_not_contains.serialize(out, child, "not_contains");
            written_bytes += m_not_contains_select1.serialize(out, child, "not_contains_select1");
            written_bytes += m_not_contains_select0.serialize(out, child, "not_contains_select0");
            written_bytes += m_not_contains_rank.serialize(out, child, "not_contains_rank");
           
            written_bytes += m_not_contained.serialize(out, child, "not_contained");
            written_bytes += m_not_contained_select1.serialize(out, child, "not_contained_select1");
            written_bytes += m_not_contained_select0.serialize(out, child, "not_contained_select0");
            written_bytes += m_not_contained_rank.serialize(out, child, "not_contained_rank");
           

            written_bytes += m_not_touches.serialize(out, child, "not_touches");
            written_bytes += m_not_touches_select1.serialize(out, child, "not_touches_select1");
            written_bytes += m_not_touches_select0.serialize(out, child, "not_touches_select0");
            written_bytes += m_not_touches_rank.serialize(out, child, "not_touches_rank");


            written_bytes += m_is_leaf.serialize(out, child, "is_leaf");
            written_bytes += m_is_leaf_select1.serialize(out, child, "is_leaf_select1");
            written_bytes += m_is_leaf_rank.serialize(out, child, "is_leaf_rank");


            written_bytes += m_run_inits.serialize(out, child, "run_inits");
            written_bytes += m_run_endings.serialize(out, child, "run_endings");
            written_bytes += m_run_limits.serialize(out, child, "run_limits");
            written_bytes += m_run_limits_select.serialize(out, child, "run_limits_select");
            written_bytes += m_run_limits_rank.serialize(out, child, "run_limits_rank");

	    written_bytes += write_member(total_leaves, out, child, "total_leaves");
	    written_bytes += write_member(m_threshold, out, child, "threshold");
	    sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream &in) {
            m_T_bv.load(in);
            m_T_bp.load(in, &m_T_bv);

	    m_T_select0.load(in, &m_T_bv);
	    m_M_bv.load(in);
	    m_M_wt.load(in);
	    m_M_bv_select1.load(in, &m_M_bv);
	    m_M_bv_select0.load(in, &m_M_bv);
	    m_M_bv_rank.load(in, &m_M_bv);

	    m_contains.load(in);
	    m_contains_select1.load(in, &m_contains);
	    m_contains_select0.load(in, &m_contains);
	    m_contains_rank.load(in, &m_contains);

	    m_contained.load(in);
	    m_contained_select1.load(in, &m_contained);
	    m_contained_select0.load(in, &m_contained);
	    m_contained_rank.load(in, &m_contained);

	    m_touches.load(in);
	    m_touches_select1.load(in, &m_touches);
	    m_touches_select0.load(in, &m_touches);
	    m_touches_rank.load(in, &m_touches);

	    m_not_contains.load(in);
	    m_not_contains_select1.load(in, &m_not_contains);
	    m_not_contains_select0.load(in, &m_not_contains);
	    m_not_contains_rank.load(in, &m_not_contains);

	    m_not_contained.load(in);
	    m_not_contained_select1.load(in, &m_not_contained);
	    m_not_contained_select0.load(in, &m_not_contained);
	    m_not_contained_rank.load(in, &m_not_contained);


	    m_not_touches.load(in);
	    m_not_touches_select1.load(in, &m_not_touches);
	    m_not_touches_select0.load(in, &m_not_touches);
	    m_not_touches_rank.load(in, &m_not_touches);

	    m_is_leaf.load(in);
	    m_is_leaf_select1.load(in, &m_is_leaf);
	    m_is_leaf_rank.load(in, &m_is_leaf);


	    m_run_inits.load(in);
	    m_run_endings.load(in);

	    m_run_limits.load(in);
	    m_run_limits_select.load(in, &m_run_limits);
	    m_run_limits_rank.load(in, &m_run_limits);
	    read_member(total_leaves,in);
	    read_member(m_threshold, in);
	    //debug_print_pars(0,10);
        }

      value_type first(value_type x) const{
	size_type pos_close = m_T_select0(x);
	size_type pos_open = m_T_bp.find_open(pos_close);
	return 2+pos_open-m_T_bp.rank(pos_open);
      }

      size_type num_children(value_type x){
	if(x>num_vertices()) return 0;
	return x - first(x)+1;
      }

      size_type height(value_type x){
	if(x>num_vertices()) return 0;
	size_type pos_open = post_to_pre(x);
	return m_T_bp.excess(pos_open)+1;
      }



      // returns the position of the opening parentheses for node with
      // postorder x
      size_type post_to_pre(value_type x) const{
	return m_T_bp.find_open(m_T_select0(x));
      }

      value_type postorder_to_leaf_id(value_type postorder) const {

	      return m_touches_rank(postorder);
      }

      value_type leaf_id_to_postorder(value_type leaf_id) const {

	      return m_is_leaf_select1(leaf_id);
      }

      value_type node(value_type i) const {
	return m_T_select0(i);
      }

      // returns the postorder id for node with close parentheses at i
      value_type postorder(value_type i) const{
	return i-m_T_bp.rank(i)+1;
      }



      value_type min_contains_q(value_type x) const {
	return contains(x,1);
      }

      value_type min_contains(value_type x) const {
	if(x> num_vertices()+1) return 0;
	value_type ones = m_contains_rank(x);
	if(ones == m_contains_rank(m_contains.size())) return 0;
	return m_contains_select1(ones+1);
      }

      value_type min_contained_in_q(value_type x) const {
	return contained_in(x,1);
      }

      value_type min_contained_in(value_type x) const {

	if(x> num_vertices()+1) return 0;
	value_type ones = m_contained_rank(x);
	if(ones == m_contained_rank(m_contained.size())) return 0;
	return m_contained_select1(ones+1);
      }


      value_type min_not_contains_q(value_type x) const {
	return not_contains(x,1);
      }

      value_type min_not_contains(value_type x) const {

	if(x> num_vertices()+1) return 0;
	value_type ones = m_not_contains_rank(x);
	if(ones == m_not_contains_rank(m_not_contains.size())) return 0;
	return m_not_contains_select1(ones+1);
      }

      value_type min_not_contained_in_q(value_type x) const {
	return not_contained_in(x,1);
      }

      value_type min_not_contained_in(value_type x) const {

	if(x> num_vertices()+1) return 0;
	value_type ones = m_not_contained_rank(x);
	if(ones == m_not_contained_rank(m_not_contained.size())) return 0;
	return m_not_contained_select1(ones+1);
      }

      value_type min_touches_q(value_type x) const {
	return touches(x,1);
      }

      value_type min_touches(value_type x) const {

	if(x> num_vertices()+1) return 0;
	value_type ones = m_touches_rank(x);
	if(ones == m_touches_rank(m_touches.size())) return 0;
	return m_touches_select1(ones+1);
      }



      value_type min_not_touches_q(value_type x) const {
	return not_touches(x,1);
      }

      value_type min_not_touches(value_type x) const {

	if(x> num_vertices()+1) return 0;
	value_type ones = m_not_touches_rank(x);
	if(ones == m_not_touches_rank(m_not_touches.size())) return 0;
	return m_not_touches_select1(ones+1);
      }

      // returns true iff node with postorder x is contained inside postorder y
      bool is_contained(value_type x, value_type y) const{
	return first(y) <= x && x <= y;
      }
      value_type contains(value_type q, value_type x) const {
	if(q>num_vertices() or x>num_vertices()) return 0;

	if( is_contained(x,q) ) return x;
	if ( x<=q ) return first(q);
	
	else return 0;

      }


      value_type not_contains(value_type q, value_type x) const{
	if(q>num_vertices() or x>num_vertices()) return 0;
	if(is_contained(x,q)){
	  return q!=num_vertices() ? q+1: 0;
	}else{
	  return x;
	}
      }
      value_type contained_in(value_type q, value_type x) const{
	if(q>num_vertices() or x>num_vertices()) return 0;
	if(x<=q) return q;
	size_type pos_open_ans = m_T_bp.lca(post_to_pre(q), post_to_pre(x));
	if(pos_open_ans==m_T_bp.size()) return 0;
	return postorder(m_T_bp.find_close(pos_open_ans));
      }
      
      value_type not_contained_in(value_type q, value_type x) const{
	if(q>num_vertices() or x>num_vertices()) return 0;
	if(!is_contained(q,x)) return x;
	size_type pos_close_x = m_T_select0(x);
	value_type preceeding = m_T_bp.rank(pos_close_x);
	if(preceeding==num_vertices()) return 0;
	size_type pos_open_ans = m_T_bp.select(preceeding+1);
	value_type ans = postorder(m_T_bp.find_close(pos_open_ans));
	return first(ans);
      }

      size_type estimate_touches(value_type x) const {

	if(x>num_vertices() or x==0) return 0;
	size_type init = first(x);

	value_type yy1 = m_M_bv_select1(init)-init+1;

	value_type yy2 = m_M_bv_select1(x+1)-x;

	return yy2-yy1;
      }

      value_type matrix_leftmost(value_type x1, value_type x2, value_type y1, value_type y2) const {
	value_type yy1 = m_M_bv_select1(y1)-y1+1;
	value_type yy2 = m_M_bv_select1(y2+1)-y2;
	value_type ans = m_M_wt.range_next_value(x1, yy1, yy2-1);
	return ans>x2 ? 0 : ans;
      }

      value_type matrix_rightmost(value_type x1, value_type x2, value_type y1, value_type y2) const {
	value_type yy1 = m_M_bv_select1(y1)-y1+1;
	value_type yy2 = m_M_bv_select1(y2+1)-y2;
	value_type ans = m_M_wt.range_prev_value(x2, yy1, yy2-1);
	return ans<x1 ? 0 : ans;
      }

      value_type matrix_upmost(value_type x1, value_type x2, value_type y1, value_type y2) const {
	if(x1>x2) return 0;
	value_type xx1 = m_M_bv_select1(x1)-x1+1;
	value_type xx2 = m_M_bv_select1(x2+1)-x2;
	value_type ans = m_M_wt.leftmost(y1,y2, xx1,xx2-1);
	if(ans == m_M_wt.size()) return 0;
	size_type pos_in_array = m_M_bv_select0(ans+1);
	ans = m_M_bv_rank(pos_in_array); 
	return ans>x2 ? 0 : ans;
      }

      value_type matrix_downmost(value_type x1, value_type x2, value_type y1, value_type y2) const {
	if(x1>x2) return 0;
	value_type xx1 = m_M_bv_select1(x1)-x1+1;
	value_type xx2 = m_M_bv_select1(x2+1)-x2;
	value_type ans = m_M_wt.rightmost(y1,y2,xx1,xx2-1);
	if(ans >= m_M_wt.size()) return 0;
	size_type pos_in_array = m_M_bv_select0(ans+1);
	ans = m_M_bv_rank(pos_in_array); 
	return ans<x1 ? 0 : ans;
      }


      value_type left_touches(value_type q, value_type x) const {


	//value_type jj = matrix_leftmost(max(x,q+1),num_vertices()-1,first(q), q);
	//return jj;
	value_type j = matrix_leftmost(x,num_vertices()-1,first(q), q);
	//value_type j = matrix_leftmost(max(x,q+1),num_vertices()-1,first(q), q);

	//either no answer, or answer is contained in q
	if(j==0 or !is_contained(j,q) ) return j;

	//try again after subtree of q
	j = matrix_leftmost(q+1, num_vertices()-1, first(q), q);
	
	return j;
      }


      value_type right_touches(value_type q, value_type x) const {


//	value_type jj = matrix_rightmost(1,min(x-1, first(q)-1) ,first(q), q);

//	value_type jjans = jj==0 ?  0 :contained_in(jj,x);

	//return jjans;

	// OLD VERSION STARTS HERE

	//find rightmost before x
	value_type j = matrix_rightmost(1,x-1,first(q), q);
	if( j==0 ) return j;
	if( !is_contained(j,q) ){
	  value_type anc = contained_in(j,x);
	  if( anc!= 0 and !is_contained(q,anc) )
	    return anc;
	}
	j = matrix_rightmost(1,first(q)-1, first(q), q);
	
	if( j==0 ) return 0;
	if( !is_contained(j,q) ){

	  value_type anc = contained_in(j,x);
	  if( anc!=0 and !is_contained(q,anc) )
	    return anc;
	}
	return 0;
      }


      value_type left_t_touches(value_type q, value_type x) const{


	//value_type jj = matrix_upmost(max(x,q+1),num_vertices()-1,first(q), q);
	//return jj;


	value_type j = matrix_upmost(x,num_vertices()-1,first(q),q);

	//either no answer, or answer is contained in q
	if(j==0 or !is_contained(j,q) ) return j;

	//try again after subtree of q
	j = matrix_upmost(q+1, num_vertices()-1,first(q),q);
	
	return j;
      }


      value_type right_t_touches(value_type q, value_type x) const {

	//find rightmost before x
	value_type j = matrix_downmost(1,x-1,first(q),q);

	if( j==0 ) return j;
	if( !is_contained(j,q) ){
	  value_type anc = contained_in(j,x);
	  if( anc!=0 and !is_contained(q,anc) )
	    return anc;
	}
	j = matrix_downmost(1,first(q)-1,first(q),q);
	
	if( j==0 ) return 0;
	if( !is_contained(j,q) ){

	  value_type anc = contained_in(j,x);
	  if( anc!=0 and !is_contained(q,anc) )
	    return anc;
	}
	return 0;
      }

      value_type touches(value_type q, value_type x) const {

	if(q>num_vertices() or x>num_vertices()) return 0;
	//not_in(q,x);
	//not_contained_in(q,x);
	value_type ans_left = left_touches(q,x);
	ans_left = ans_left==0? num_vertices()+1 : ans_left;
	value_type ans_right = right_touches(q,x);
	ans_right = ans_right==0? num_vertices()+1 : ans_right;

	value_type ans_left_t = left_t_touches(q,x);
	ans_left_t = ans_left_t==0? num_vertices()+1 : ans_left_t;
	value_type ans_right_t = right_t_touches(q,x);
	ans_right_t = ans_right_t==0? num_vertices()+1 : ans_right_t;

	value_type ans = min({ans_left,ans_right,ans_left_t,ans_right_t});
	ans = ans==num_vertices()+1 ? 0 : ans;
	return ans;
      }

      value_type not_touches(value_type q, value_type x) const {
	//cout << "querying not_touches("<<q <<", "<<x<<")\n";
	for(int k=x; k<x+m_threshold; k++){

	  if(touches(q,k) != k){
	    return k;
	  }
	}

	// every node in [k, k+t] touches q, find run and return end+1;
	//
	
	int l = m_run_limits_select(q);
	l = l-m_run_limits_rank(l);

	int r = m_run_limits_select(q+1);
	r = r-m_run_limits_rank(r);
	//cout << "l,r = "<<l <<" " <<r <<endl;

	auto it = lower_bound(m_run_endings.begin()+l, m_run_endings.begin()+r, x);

	int pos = it - m_run_endings.begin();

	return m_run_endings[pos] + 1;

	
	return 0;

      }
      
      void debug_print_pars(int i, int r) const {

	      for(int j=i; j<i+r; j++){
		      cerr << m_T_bv[j];
	      }
	      cerr << endl;
      }

      void debug_print_runs() const {

	
	int r_idx = 0;

	for(int i=0; i<m_run_limits.size(); i++){
	  printf("%3d ", (int)m_run_limits[i]);
	}
	printf("\n");

	for(int i=0; i<m_run_limits.size(); i++){
	  if(m_run_limits[i]==0)
	    printf("%3d ", (int)m_run_inits[r_idx++]);
	  else
	    printf("    ");
	}
	printf("\n");


	r_idx = 0;
	for(int i=0; i<m_run_limits.size(); i++){
	  if(m_run_limits[i]==0)
	    printf("%3d ", (int)m_run_endings[r_idx++]);
	  else
	    printf("    ");
	}
	printf("\n");



      }





  };


}

#endif
