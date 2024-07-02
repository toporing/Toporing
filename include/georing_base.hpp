#ifndef GEORING_BASE_HPP
#define GEORING_BASE_HPP

#include "configuration.hpp" 

using namespace std;

namespace ring {

  class georing_base {
    public:

      typedef uint64_t size_type;
      typedef uint64_t value_type;

      virtual size_type num_vertices() const = 0;

      virtual value_type min_contains(value_type x) const { return -1; }
      virtual value_type min_contained_in(value_type x) const { return -1; }
      virtual value_type min_contains_q(value_type x) const { return -1; }
      virtual value_type min_contained_in_q(value_type x) const { return -1; }
      virtual value_type min_touches(value_type x) const { return -1; }
      virtual value_type min_touches_q(value_type x) const { return -1; }
      virtual value_type touches(value_type x) const { return -1; }
      virtual value_type contains(value_type x) const { return -1; }
      virtual value_type contained_in(value_type x) const { return -1; }

  };
}
#endif
