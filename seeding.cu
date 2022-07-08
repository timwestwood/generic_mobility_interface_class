// seeding.cu

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <string>
#include "seeding.hpp"
#include "cuda_functions.hpp"
#include "matrix.hpp"
#include "config.hpp"

void repulsive_blob_seeding(double *const blob_references){

  std::cout << std::endl << std::endl << "Input file " << "'sphere"+std::to_string(NBLOB)+".seed'" << " was not found or could not be opened." << std::endl;
  std::cout << "Performing repulsive simulation to seed the blobs..." << std::endl;

  #define SEED_ON_CPU (NBLOB<1000)

  double phi;

  for (int n = 0; n < NBLOB; n++){

    blob_references[3*n + 2] = (NBLOB == 1) ? -1.0 : 2.0*n/double(NBLOB-1) - 1.0;

    const double r = std::sqrt(1.0 - blob_references[3*n + 2]*blob_references[3*n + 2]);

    if ((n == 0) || (n == NBLOB-1)){

      phi = 0.0;

    } else {

      phi += 3.6/(r*std::sqrt(NBLOB));

    }

    blob_references[3*n] = r*std::cos(phi);
    blob_references[3*n + 1] = r*std::sin(phi);

  }

  #if SEED_ON_CPU

    matrix F(3*NBLOB, 1);

  #else

    // Using unified memory as this shouldn't need to be optimised
    double *X, *F;
    cudaMallocManaged(&X, 3*NBLOB*sizeof(double));
    cudaMallocManaged(&F, 3*NBLOB*sizeof(double));

    const int num_thread_blocks = (NBLOB + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

  #endif

  // Make sure it always enters the loop
  double biggest_move_size = 1e4;

  while (biggest_move_size > 1e-4){

    biggest_move_size = 0.0;

    // There are smarter ways of doing this, but it works and was quick to write...

    #if SEED_ON_CPU

      F.zero();

      for (int n = 0; n < NBLOB-1; n++){

        const double Xn[3] = {blob_references[3*n], blob_references[3*n + 1], blob_references[3*n + 2]};

        for (int m = n+1; m < NBLOB; m++){

          const double d[3] = {Xn[0] - blob_references[3*m], Xn[1] - blob_references[3*m + 1], Xn[2] - blob_references[3*m + 2]};
          const double d_norm = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);

          F(3*n) += 5.0*(2.0 - d_norm)*d[0]/d_norm;
          F(3*n + 1) += 5.0*(2.0 - d_norm)*d[1]/d_norm;
          F(3*n + 2) += 5.0*(2.0 - d_norm)*d[2]/d_norm;

          F(3*m) -= 5.0*(2.0 - d_norm)*d[0]/d_norm;
          F(3*m + 1) -= 5.0*(2.0 - d_norm)*d[1]/d_norm;
          F(3*m + 2) -= 5.0*(2.0 - d_norm)*d[2]/d_norm;

        }
      }

      for (int n = 0; n < NBLOB; n++){

        double new_pos[3] = {blob_references[3*n] + F(3*n), blob_references[3*n + 1] + F(3*n + 1), blob_references[3*n + 2] + F(3*n + 2)};
        const double norm = std::sqrt(new_pos[0]*new_pos[0] + new_pos[1]*new_pos[1] + new_pos[2]*new_pos[2]);
        new_pos[0] /= norm;
        new_pos[1] /= norm;
        new_pos[2] /= norm;

        biggest_move_size = std::max<double>(std::abs(new_pos[0] - blob_references[3*n]), biggest_move_size);
        biggest_move_size = std::max<double>(std::abs(new_pos[1] - blob_references[3*n + 1]), biggest_move_size);
        biggest_move_size = std::max<double>(std::abs(new_pos[2] - blob_references[3*n + 2]), biggest_move_size);

        blob_references[3*n] = new_pos[0];
        blob_references[3*n + 1] = new_pos[1];
        blob_references[3*n + 2] = new_pos[2];

      }

    #else

      for (int n = 0; n < 3*NBLOB; n++){

        X[n] = blob_references[n];

      }

      seeding_repulsion_step<<<num_thread_blocks,THREADS_PER_BLOCK>>>(F, X, NBLOB);

      cudaDeviceSynchronize();

      for (int n = 0; n < NBLOB; n++){

        double new_pos[3] = {X[3*n] + F[3*n], X[3*n + 1] + F[3*n + 1], X[3*n + 2] + F[3*n + 2]};
        const double norm = std::sqrt(new_pos[0]*new_pos[0] + new_pos[1]*new_pos[1] + new_pos[2]*new_pos[2]);
        new_pos[0] /= norm;
        new_pos[1] /= norm;
        new_pos[2] /= norm;

        biggest_move_size = std::max<double>(std::abs(new_pos[0] - blob_references[3*n]), biggest_move_size);
        biggest_move_size = std::max<double>(std::abs(new_pos[1] - blob_references[3*n + 1]), biggest_move_size);
        biggest_move_size = std::max<double>(std::abs(new_pos[2] - blob_references[3*n + 2]), biggest_move_size);

        blob_references[3*n] = new_pos[0];
        blob_references[3*n + 1] = new_pos[1];
        blob_references[3*n + 2] = new_pos[2];

      }

    #endif

  }

  #if !SEED_ON_CPU

    cudaFree(X);
    cudaFree(F);

  #endif

  std::ofstream blob_ref_file("sphere"+std::to_string(NBLOB)+".seed");

  for (int n = 0; n < 3*NBLOB; n++){

    blob_ref_file << blob_references[n] << " ";

  }

  blob_ref_file.close();

  std::cout << "...done!" << std::endl;

};

void repulsive_filament_seeding(double *const filament_references){

  #if EQUATORIAL_SEEDING

    std::cout << std::endl << std::endl << "Input file " << "'sphere"+std::to_string(NFIL)+"_equatorial.seed'" << " was not found or could not be opened." << std::endl;

    std::cout << "Seeding locations according to a spiral pattern..." << std::endl;

    const double theta_max = 0.5/double(R_OVER_L); // = 0.5 * 2*PI*L/(2*PI*R), meaning the width of the band is L as measured across the sphere surface.
    const double h = 2.0*std::sin(theta_max);
    const int num_fils_for_sphere = std::ceil(NFIL*2.0/h);
    const int offset = std::round(0.5*(num_fils_for_sphere - NFIL));

    double phi = 0.0;

    for (int n = 0; n < NFIL; n++){

      filament_references[3*n + 2] = 2.0*(n + offset)/double(num_fils_for_sphere - 1) - 1.0;
      const double r = std::sqrt(1.0 - filament_references[3*n + 2]*filament_references[3*n + 2]);
      filament_references[3*n] = r*std::cos(phi);
      filament_references[3*n + 1] = r*std::sin(phi);

      phi += 3.6/(r*std::sqrt(num_fils_for_sphere));

    }

  #elif PLATY_SEEDING

    std::cout << std::endl << std::endl << "Input file " << "'sphere"+std::to_string(NFIL)+"_platy.seed'" << " was not found or could not be opened." << std::endl;

    std::cout << "Seeding..." << std::endl;

    const double theta_max = 0.5/double(R_OVER_L); // = 0.5 * 2*PI*L/(2*PI*R), meaning the width of the band is L as measured across the sphere surface.
    const double h = 2.0*std::sin(theta_max);
    const int num_fils_in_main_band = std::round(0.9*NFIL); // Put 90% in the main band
    const int num_fils_for_sphere = std::ceil(num_fils_in_main_band*2.0/h);
    const int offset = std::round(0.5*(num_fils_for_sphere - num_fils_in_main_band));

    // Seed the main band
    double phi = 0.0;

    for (int n = 0; n < num_fils_in_main_band; n++){

      filament_references[3*n + 2] = 2.0*(n + offset)/double(num_fils_for_sphere - 1) - 1.0;
      const double r = std::sqrt(1.0 - filament_references[3*n + 2]*filament_references[3*n + 2]);
      filament_references[3*n] = r*std::cos(phi);
      filament_references[3*n + 1] = r*std::sin(phi);

      phi += 3.6/(r*std::sqrt(num_fils_for_sphere));

    }

    // Seed the remaining 10% at the bottom of the swimmer
    for (int n = num_fils_in_main_band; n < NFIL; n++){

      filament_references[3*n + 2] = -0.9; // Some fixed location for the bottom band
      const double r = std::sqrt(1.0 - filament_references[3*n + 2]*filament_references[3*n + 2]);
      phi = 2.0*PI*double(n - num_fils_in_main_band)/double(NFIL - num_fils_in_main_band);
      filament_references[3*n] = r*std::cos(phi);
      filament_references[3*n + 1] = r*std::sin(phi);

    }

  #elif UNIFORM_SEEDING

    std::cout << std::endl << std::endl << "Input file " << "'sphere"+std::to_string(NFIL)+".seed'" << " was not found or could not be opened." << std::endl;

    std::cout << "Performing repulsive simulation to seed the filaments..." << std::endl;

    #define SEED_ON_CPU (NFIL<1000)

    double phi;

    for (int n = 0; n < NFIL; n++){

      filament_references[3*n + 2] = (NFIL == 1) ? -1.0 : 2.0*n/double(NFIL-1) - 1.0;

      const double r = std::sqrt(1.0 - filament_references[3*n + 2]*filament_references[3*n + 2]);

      if ((n == 0) || (n == NFIL-1)){

        phi = 0.0;

      } else {

        phi += 3.6/(r*std::sqrt(NFIL));

      }

      filament_references[3*n] = r*std::cos(phi);
      filament_references[3*n + 1] = r*std::sin(phi);

    }

    #if SEED_ON_CPU

      matrix F(3*NFIL, 1);

    #else

      // Using unified memory as this shouldn't need to be optimised
      double *X, *F;
      cudaMallocManaged(&X, 3*NFIL*sizeof(double));
      cudaMallocManaged(&F, 3*NFIL*sizeof(double));

      const int num_thread_blocks = (NFIL + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    #endif

    // Make sure it always enters the loop
    double biggest_move_size = 1e4;

    while (biggest_move_size > 1e-4){

      biggest_move_size = 0.0;

      // There are smarter ways of doing this, but it works and was quick to write...

      #if SEED_ON_CPU

        F.zero();

        for (int n = 0; n < NFIL; n++){

          const double Xn[3] = {filament_references[3*n], filament_references[3*n + 1], filament_references[3*n + 2]};

          for (int m = n+1; m < NFIL; m++){

            const double d[3] = {Xn[0] - filament_references[3*m], Xn[1] - filament_references[3*m + 1], Xn[2] - filament_references[3*m + 2]};
            const double d_norm = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);

            F(3*n) += 5.0*(2.0 - d_norm)*d[0]/d_norm;
            F(3*n + 1) += 5.0*(2.0 - d_norm)*d[1]/d_norm;
            F(3*n + 2) += 5.0*(2.0 - d_norm)*d[2]/d_norm;

            F(3*m) -= 5.0*(2.0 - d_norm)*d[0]/d_norm;
            F(3*m + 1) -= 5.0*(2.0 - d_norm)*d[1]/d_norm;
            F(3*m + 2) -= 5.0*(2.0 - d_norm)*d[2]/d_norm;

          }
        }

        for (int n = 0; n < NFIL; n++){

          double new_pos[3] = {filament_references[3*n] + F(3*n), filament_references[3*n + 1] + F(3*n + 1), filament_references[3*n + 2] + F(3*n + 2)};
          const double norm = sqrt(new_pos[0]*new_pos[0] + new_pos[1]*new_pos[1] + new_pos[2]*new_pos[2]);
          new_pos[0] /= norm;
          new_pos[1] /= norm;
          new_pos[2] /= norm;

          biggest_move_size = std::max<double>(std::abs(new_pos[0] - filament_references[3*n]), biggest_move_size);
          biggest_move_size = std::max<double>(std::abs(new_pos[1] - filament_references[3*n + 1]), biggest_move_size);
          biggest_move_size = std::max<double>(std::abs(new_pos[2] - filament_references[3*n + 2]), biggest_move_size);

          filament_references[3*n] = new_pos[0];
          filament_references[3*n + 1] = new_pos[1];
          filament_references[3*n + 2] = new_pos[2];

        }

      #else

        for (int n = 0; n < 3*NFIL; n++){

          X[n] = filament_references[n];

        }

        seeding_repulsion_step<<<num_thread_blocks,THREADS_PER_BLOCK>>>(F, X, NFIL);

        cudaDeviceSynchronize();

        for (int n = 0; n < NFIL; n++){

          double new_pos[3] = {X[3*n] + F[3*n], X[3*n + 1] + F[3*n + 1], X[3*n + 2] + F[3*n + 2]};
          const double norm = std::sqrt(new_pos[0]*new_pos[0] + new_pos[1]*new_pos[1] + new_pos[2]*new_pos[2]);
          new_pos[0] /= norm;
          new_pos[1] /= norm;
          new_pos[2] /= norm;

          biggest_move_size = std::max<double>(std::abs(new_pos[0] - filament_references[3*n]), biggest_move_size);
          biggest_move_size = std::max<double>(std::abs(new_pos[1] - filament_references[3*n + 1]), biggest_move_size);
          biggest_move_size = std::max<double>(std::abs(new_pos[2] - filament_references[3*n + 2]), biggest_move_size);

          filament_references[3*n] = new_pos[0];
          filament_references[3*n + 1] = new_pos[1];
          filament_references[3*n + 2] = new_pos[2];

        }

      #endif

    }

    #if !SEED_ON_CPU

      cudaFree(X);
      cudaFree(F);

    #endif

  #endif

  std::ofstream fil_ref_file;

  #if EQUATORIAL_SEEDING

    fil_ref_file.open("sphere"+std::to_string(NFIL)+"_equatorial.seed");

  #elif UNIFORM_SEEDING

    fil_ref_file.open("sphere"+std::to_string(NFIL)+".seed");

  #endif

  for (int n = 0; n < 3*NFIL; n++){

    fil_ref_file << filament_references[n] << " " << "\n";

  }

  fil_ref_file.close();

  std::cout << "...done!" << std::endl;

};
