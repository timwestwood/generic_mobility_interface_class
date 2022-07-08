// rpy_mobility_solver.hpp

// =============================================================================
// Include guard
#ifndef MY_RPY_MOBILITY_SOLVER_HEADER_INCLUDED
#define MY_RPY_MOBILITY_SOLVER_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies
class swimmer;

// =============================================================================
// Included dependencies
#include <vector>
#include "matrix.hpp"
#include "config.hpp"

class rpy_mobility_solver{

public:

  // GPU info
  int num_gpus;

  // CUDA variables
  double *v_segs_host;
  double **v_segs_device;

  double *v_blobs_host;
  double **v_blobs_device;

  double *x_segs_host;
  double **x_segs_device;

  double *x_blobs_host;
  double **x_blobs_device;

  double *f_segs_host;
  double **f_segs_device;

  double *f_blobs_host;
  double **f_blobs_device;

  double *f_blobs_repulsion_host;
  double **f_blobs_repulsion_device;

  int *num_segs;
  int *num_blobs;

  ~rpy_mobility_solver();
  rpy_mobility_solver();

  void read_positions_and_forces(std::vector<swimmer>& swimmers);
  void compute_velocities(std::vector<swimmer>& swimmers, int& num_gmres_iterations, const int nt);
  bool compute_errors(matrix& error, const std::vector<swimmer>& swimmers, const int nt);
  void write_data(const int nt, const std::vector<swimmer>& swimmers);

  #if !USE_BROYDEN_FOR_EVERYTHING

    // linear system storage
    matrix rhs;
    matrix v_bodies;
    matrix Q;
    matrix H;
    matrix beta;
    matrix SN;
    matrix CS;

    // Reference matrices for preconditioning
    matrix body_mobility_reference;
    matrix KTKinv_reference;

    void assemble_rhs(const std::vector<swimmer>& swimmers, const int nt);
    matrix apply_preconditioner(const matrix& in, const std::vector<swimmer>& swimmers);
    matrix system_matrix_mult(const matrix& in, const std::vector<swimmer>& swimmers);
    int solve_linear_system(std::vector<swimmer>& swimmers);
    void make_body_reference_matrices(std::vector<swimmer>& swimmers);

  #endif

  #if DYNAMIC_PHASE_EVOLUTION

    std::vector<double> gen_force_refs;

  #endif


}; // End of rpy_mobility_solver class definition

#endif // MY_RPY_MOBILITY_SOLVER_HEADER_INCLUDED
