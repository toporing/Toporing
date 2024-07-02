/*
  Stats of the dataset

  107,836,911 subjects
  242,124,917 objects
  5,419 predicates
 */

#include <toporing.hpp>
#include <chrono>
#include <time.hpp>

using namespace ::util::time;
using namespace std;
using namespace ring;


int main(int argc, char* argv[]) {

    if(argc < 3){
        std::cout << "Usage: " << argv[0] << " <index> <# of queries>" << std::endl;
        return 0;
    }

    std::string index = argv[1];
    int num_queries = atoi(argv[2]);
    int num_S = 107836911;
    int num_O = 242124917;
    
    toporing<> graph;
    
    sdsl::load_from_file(graph, index);

    int *S_ids = new int[num_queries];
    int *O_ids = new int[num_queries];

    for(int i=0; i < num_queries; i++) {
      S_ids[i] = rand() % num_S;
      O_ids[i] = rand() % num_O;      
    }

    ::util::time::usage::usage_type start, stop;
    uint64_t total_elapsed_time;
    uint64_t total_user_time;


    int dump = 0, x;
    // Warm-up
    for(int i=0; i < num_queries; i++)
      x = graph.touches(S_ids[i], O_ids[i]);
    dump += x;

    
    cout << "query;total_elapsed_time" << endl;
    
    // Contains
    start = ::util::time::usage::now();
    for(int i=0; i < num_queries; i++)
      x = graph.contains(S_ids[i], O_ids[i]);
    stop = ::util::time::usage::now();
    dump += x;
    total_elapsed_time = (uint64_t) duration_cast<nanoseconds>(stop.elapsed - start.elapsed);
    cout << "contains;" << total_elapsed_time << endl;

    // Contained_in
    start = ::util::time::usage::now();
    for(int i=0; i < num_queries; i++)
      x = graph.contained_in(S_ids[i], O_ids[i]);
    stop = ::util::time::usage::now();
    dump += x;
    total_elapsed_time = (uint64_t) duration_cast<nanoseconds>(stop.elapsed - start.elapsed);
    cout << "contained_in;" << total_elapsed_time << endl;
    
    // Touches
    start = ::util::time::usage::now();
    for(int i=0; i < num_queries; i++)
      x = graph.touches(S_ids[i], O_ids[i]);
    stop = ::util::time::usage::now();
    dump += x;
    total_elapsed_time = (uint64_t) duration_cast<nanoseconds>(stop.elapsed - start.elapsed);
    cout << "touches;" << total_elapsed_time << endl;
    
    // Not contains
    start = ::util::time::usage::now();
    for(int i=0; i < num_queries; i++)
      x = graph.not_contains(S_ids[i], O_ids[i]);
    stop = ::util::time::usage::now();
    dump += x;
    total_elapsed_time = (uint64_t) duration_cast<nanoseconds>(stop.elapsed - start.elapsed);
    cout << "not_contains;" << total_elapsed_time << endl;

    // Not contained_in
    start = ::util::time::usage::now();
    for(int i=0; i < num_queries; i++)
      x = graph.not_contained_in(S_ids[i], O_ids[i]);
    stop = ::util::time::usage::now();
    dump += x;
    total_elapsed_time = (uint64_t) duration_cast<nanoseconds>(stop.elapsed - start.elapsed);
    cout << "not_contained_in;" << total_elapsed_time << endl;
    
    // Not touches
    start = ::util::time::usage::now();
    for(int i=0; i < num_queries; i++)
      x = graph.not_touches(S_ids[i], O_ids[i]);
    stop = ::util::time::usage::now();
    dump += x;
    total_elapsed_time = (uint64_t) duration_cast<nanoseconds>(stop.elapsed - start.elapsed);
    cout << "not_touches;" << total_elapsed_time << endl;
    
    return 0;
}
