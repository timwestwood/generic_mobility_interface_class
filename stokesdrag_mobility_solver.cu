// stokesdrag_mobility_solver.cu

#include <iomanip>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "stokesdrag_mobility_solver.hpp"
#include "cuda_functions.hpp"
#include "swimmer.hpp"

stokesdrag_mobility_solver::~stokesdrag_mobility_solver(){}

stokesdrag_mobility_solver::stokesdrag_mobility_solver(){}

void stokesdrag_mobility_solver::evaluate_segment_segment_mobility(){

  int start_seg = 0;

    for (int n = 0; n < num_gpus; n++){

      cudaSetDevice(n);

      int num_thread_blocks = (num_segs[n] + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

      Ms_mult<<<num_thread_blocks, THREADS_PER_BLOCK>>>(v_segs_device[n], f_segs_device[n], start_seg, num_segs[n]);

      start_seg += num_segs[n];

    }

}

void stokesdrag_mobility_solver::evaluate_segment_blob_mobility(){

  int start_seg = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    int num_thread_blocks = (num_segs[n] + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    Ms_fill_zero<<<num_thread_blocks, THREADS_PER_BLOCK>>>(v_segs_device[n], start_seg, num_segs[n]);

    start_seg += num_segs[n];

  }

}

void stokesdrag_mobility_solver::evaluate_blob_blob_mobility(){

  int start_blob = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    int num_thread_blocks = (num_blobs[n] + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    Mb_mult<<<num_thread_blocks, THREADS_PER_BLOCK>>>(v_blobs_device[n], f_blobs_device[n], start_blob, num_blobs[n]);

    start_blob += num_blobs[n];

  }

}

void stokesdrag_mobility_solver::evaluate_blob_segment_mobility(){

  int start_blob = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    const int num_thread_blocks = (num_blobs[n] + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    Mb_fill_zero<<<num_thread_blocks, THREADS_PER_BLOCK>>>(v_blobs_device[n], start_blob, num_blobs[n]);

    start_blob += num_blobs[n];

  }

}

