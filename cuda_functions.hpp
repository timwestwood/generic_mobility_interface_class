// cuda_functions.hpp

// =============================================================================
// Include guard
#ifndef MY_CUDA_FUNCTIONS_HEADER_INCLUDED
#define MY_CUDA_FUNCTIONS_HEADER_INCLUDED

__global__ void Mss_mult(double *V, const double *const F, const double *const X, const int start_seg, const int num_segs);
__global__ void Mbb_mult(double *V, const double *const F, const double *const X, const int start_blob, const int num_blobs);
__global__ void Msb_mult(double *V, const double *const F, const double *const Xs, const double *const Xb, const int start_seg, const int num_segs);
__global__ void Mbs_mult(double *V, const double *const F, const double *const Xb, const double *const Xs, const int start_blob, const int num_blobs);
__global__ void Ms_mult(double *V, const double *const F, const int start_seg, const int num_segs);
__global__ void Mb_mult(double *V, const double *const F, const int start_blob, const int num_blobs);
__global__ void Ms_fill_zero(double *V, const int start_seg, const int num_segs);
__global__ void Mb_fill_zero(double *V, const int start_blob, const int num_blobs);
__global__ void barrier_forces(double *f_segs, double *f_blobs_repulsion, const double *const x_segs, const double *const x_blobs, const int start_seg, const int num_segs, const int start_blob, const int num_blobs);
__global__ void seeding_repulsion_step(double *F, const double *X, const int N);

#endif // MY_CUDA_FUNCTIONS_HEADER_INCLUDED
