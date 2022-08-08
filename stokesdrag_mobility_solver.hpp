// stokesdrag_mobility_solver.hpp

// =============================================================================
// Include guard
#ifndef MY_STOKESDRAG_MOBILITY_SOLVER_HEADER_INCLUDED
#define MY_STOKESDRAG_MOBILITY_SOLVER_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies

// =============================================================================
// Included dependencies
#include "rpy_mobility_solver.hpp"

class stokesdrag_mobility_solver : public rpy_mobility_solver{

public:

  // =============================================================================
  // Everything we have to define for the base class:

  void evaluate_segment_segment_mobility();
  void evaluate_segment_blob_mobility();
  void evaluate_blob_blob_mobility();
  void evaluate_blob_segment_mobility();

  ~stokesdrag_mobility_solver();
  stokesdrag_mobility_solver();

}; // End of stokesdrag_mobility_solver class definition

#endif // MY_RPY_MOBILITY_SOLVER_HEADER_INCLUDED
