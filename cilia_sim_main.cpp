#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include "swimmer.hpp"
#include "rpy_mobility_solver.hpp"

#if !(PRESCRIBED_CILIA || NO_CILIA_SQUIRMER)

  #include "broyden_solver.hpp"

#endif

#include "config.hpp"

int main(int argc, char** argv){

  double *data_from_file;

  #if READ_INITIAL_CONDITIONS_FROM_BACKUP

    std::cout << std::endl << std::endl << "Reading initial conditions from backup file..." << std::endl;

    std::ifstream input_file(INITIAL_CONDITIONS_FILE_NAME+std::string(".backup"));

    int nt_of_backup;

    #if PRESCRIBED_CILIA

      const int data_per_fil = 10;

    #else

      const int data_per_fil = 9 + 28*NSEG;

    #endif

    #if USE_BROYDEN_FOR_EVERYTHING

      const int data_per_swimmer = 19 + 9*NBLOB + NFIL*data_per_fil;

    #else

      const int data_per_swimmer = 19 + NFIL*data_per_fil;

    #endif

    const int backup_file_length = NSWIM*data_per_swimmer;

    data_from_file = new double[backup_file_length];

    input_file >> nt_of_backup;

    for (int fpos = 0; fpos < backup_file_length; fpos++){

      input_file >> data_from_file[fpos];

    }

    std::cout << std::endl << "Finished reading from file!" << std::endl;

    input_file.close();

  #else

    const int data_per_swimmer = 0;

  #endif

  std::ofstream config_file(SIMULATION_CONFIG_NAME);

  config_file << NFIL << " " << "%% NFIL" << std::endl;
  config_file << NSEG << " " << "%% NSEG" << std::endl;
  config_file << NSWIM << " " << "%% NSWIM" << std::endl;
  config_file << NBLOB << " " << "%% NBLOB" << std::endl;
  config_file << MU << " " << "%% MU" << std::endl;
  config_file << RSEG << " " << "%% RSEG" << std::endl;
  config_file << RBLOB << " " << "%% RBLOB" << std::endl;
  config_file << KB << " " << "%% KB" << std::endl;
  config_file << KT << " " << "%% KT" << std::endl;
  config_file << TOL << " " << "%% TOL" << std::endl;
  config_file << TOTAL_TIME_STEPS << " " << "%% TOTAL_TIME_STEPS" << std::endl;
  config_file << DT << " " << "%% DT" << std::endl;
  config_file << PLOT_FREQUENCY_IN_STEPS << " " << "%% PLOT_FREQUENCY_IN_STEPS" << std::endl;

  #if INSTABILITY_CILIA

    config_file << DIMENSIONLESS_FORCE << " " << "%% DIMENSIONLESS_FORCE" << std::endl;

  #elif GEOMETRIC_CILIA

    config_file << DRIVING_TORQUE_MAGNITUDE_RATIO << " " << "%% DRIVING_TORQUE_MAGNITUDE_RATIO" << std::endl;
    config_file << DEFAULT_BASE_TORQUE_MAGNITUDE << " " << "%% DEFAULT_BASE_TORQUE_MAGNITUDE" << std::endl;
    config_file << CRITICAL_ANGLE << " " << "%% CRITICAL_ANGLE" << std::endl;
    config_file << TRANSITION_STEPS << " " << "%% TRANSITION_STEPS" << std::endl;

  #elif CONSTANT_BASE_ROTATION

    config_file << BASE_ROTATION_RATE << " " << "%% BASE_ROTATION_RATE" << std::endl;

  #elif (PRESCRIBED_CILIA || NO_CILIA_SQUIRMER)

    config_file << STEPS_PER_PERIOD << " " << "%% STEPS_PER_PERIOD" << std::endl;

  #endif

  config_file.close();

  // Initialise the simulation
  rpy_mobility_solver mobility;

  #if !(PRESCRIBED_CILIA || NO_CILIA_SQUIRMER)

    broyden_solver broyden;

  #endif

  std::vector<swimmer> swimmers(NSWIM);

  for (int n = 0; n < NSWIM; n++){

    swimmers[n].initial_setup(n, &data_from_file[n*data_per_swimmer],
                                  &mobility.x_segs_host[3*n*NFIL*NSEG],
                                  &mobility.f_segs_host[6*n*NFIL*NSEG],
                                  &mobility.f_blobs_host[3*n*NBLOB]);

  }

  #if READ_INITIAL_CONDITIONS_FROM_BACKUP

    delete[] data_from_file;

  #endif

  // Swimmers are all identical
  swimmers[0].write_reference_positions();

  #if !(PRESCRIBED_BODY_VELOCITIES || PRESCRIBED_CILIA || USE_BROYDEN_FOR_EVERYTHING)

    mobility.make_body_reference_matrices(swimmers);

    for (int n = 0; n < NSWIM; n++){

      swimmers[n].body_mobility_reference = mobility.body_mobility_reference;

    }

  #endif

  #if INFINITE_PLANE_WALL

    std::cout << std::endl;
    std::cout << "Simulating " << NFIL << " filaments on an infinite no-slip wall, with each filament comprised of " << NSEG << " segments." << std::endl;
    std::cout << std::endl;

  #else

    std::cout << std::endl;
    std::cout << "Simulating " << NSWIM << " swimmers, each having a rigid body resolved using " << NBLOB << " 'blobs'." << std::endl;
    std::cout << "Attached to each rigid body are " << NFIL << " filaments, each comprised of " << NSEG << " segments." << std::endl;
    std::cout << std::endl;

  #endif

  #if WRITE_GENERALISED_FORCES

    std::ofstream generalised_force_file(reference_generalised_force_file_name());

    generalised_force_file << TOTAL_TIME_STEPS << " ";
    generalised_force_file << std::scientific << std::setprecision(15);

  #endif

  // If continuing from backup file...
  #if INITIAL_CONDITIONS_TYPE==1

    const int nt_start = nt_of_backup;

  #else

    const int nt_start = 0;

    //swimmers[0].filaments[1].phase = 0.5*PI;

  #endif

  // Begin time stepping
  for (int nt = nt_start; nt < TOTAL_TIME_STEPS; nt++) {

    for (int i = 0; i < NSWIM; i++) {

      swimmers[i].initial_guess(nt);

      swimmers[i].forces_and_torques(nt);

    }

    int num_gmres_iterations;
    mobility.compute_velocities(swimmers, num_gmres_iterations, nt);

    #if (PRESCRIBED_CILIA || NO_CILIA_SQUIRMER)

      // Explicit time integration for these simulation types.
      for (int n = 0; n < NSWIM; n++){

        swimmers[n].body.x[0] += DT*mobility.v_bodies(6*n);
        swimmers[n].body.x[1] += DT*mobility.v_bodies(6*n + 1);
        swimmers[n].body.x[2] += DT*mobility.v_bodies(6*n + 2);

        swimmers[n].body.u[0] = DT*mobility.v_bodies(6*n + 3);
        swimmers[n].body.u[1] = DT*mobility.v_bodies(6*n + 4);
        swimmers[n].body.u[2] = DT*mobility.v_bodies(6*n + 5);
        swimmers[n].body.q = lie_exp(swimmers[n].body.u)*swimmers[n].body.qm1;

      }

      #if WRITE_GENERALISED_FORCES

        double generalised_force = 0.0;

        for (int n = 0; n < NSEG; n++){

          generalised_force += swimmers[0].filaments[0].vel_dir[3*n]*mobility.f_segs_host[6*n];
          generalised_force += swimmers[0].filaments[0].vel_dir[3*n + 1]*mobility.f_segs_host[6*n + 1];
          generalised_force += swimmers[0].filaments[0].vel_dir[3*n + 2]*mobility.f_segs_host[6*n + 2];

        }

        generalised_force_file << generalised_force << " ";

      #endif

    #else

      bool error_is_too_large = mobility.compute_errors(broyden.error, swimmers, nt);

      if (error_is_too_large){

        for (int i = 0; i < NSWIM; i++) {

          swimmers[i].prepare_jacobian_inv(nt);

        }

      }

      broyden.iter = 0;

      while (error_is_too_large && (broyden.iter < MAX_BROYDEN_ITER)){

        broyden.find_update(swimmers, nt);

        for (int i = 0; i < NSWIM; i++) {

          #if INFINITE_PLANE_WALL

            const int per_body = 6*NFIL*NSEG;

          #elif USE_BROYDEN_FOR_EVERYTHING

            const int per_body = 6*NFIL*NSEG + 3*NBLOB + 6;

          #else

            const int per_body = 6*NFIL*NSEG + 6;

          #endif

          swimmers[i].update(&broyden.update.data[i*per_body]);

          swimmers[i].forces_and_torques(nt);

        }

        mobility.compute_velocities(swimmers, num_gmres_iterations, nt);

        error_is_too_large = mobility.compute_errors(broyden.new_error, swimmers, nt);

        if (!broyden.new_error.is_finite()){

          std::cout << DELETE_CURRENT_LINE << std::flush;
          std::cout << "Broyden's method diverged after " << broyden.iter+1 << " iterations during step " << nt+1 << "." << std::endl;

          return 0;

        }

        broyden.end_of_iter(swimmers, nt, nt_start, error_is_too_large);

        std::cout << DELETE_CURRENT_LINE << std::flush;
        std::cout << "Step " << nt+1 << ": Completed Broyden iteration " << broyden.iter;
        #if !(INFINITE_PLANE_WALL || USE_BROYDEN_FOR_EVERYTHING)
          std::cout << " in " << num_gmres_iterations << " iterations of the linear system solver." << "\r";
        #endif
        std::cout << std::flush;

      }

    #endif

    for (int i = 0; i < NSWIM; i++) {

      swimmers[i].end_of_step(nt);

    }

    #if (PRESCRIBED_CILIA || NO_CILIA_SQUIRMER)

      std::cout << "Completed step " << nt+1 << "/" << TOTAL_TIME_STEPS << " in " << num_gmres_iterations << " iterations of the linear system solver." << std::endl;

    #else

      if (error_is_too_large){

        std::cout << DELETE_CURRENT_LINE << std::flush;
        std::cout << "Broyden's method failed to converge within " << MAX_BROYDEN_ITER << " iterations during step " << nt+1 << "." << std::endl;

        return 0;

      } else {

        std::cout << DELETE_CURRENT_LINE << std::flush;
        std::cout << "Completed step " << nt+1 << "/" << TOTAL_TIME_STEPS << " in " << broyden.iter << " Broyden iterations (max = " << broyden.max_iter << ", avg. = " << broyden.avg_iter << ")." << std::endl;

      }

    #endif

    if (nt%PLOT_FREQUENCY_IN_STEPS == 0){

      std::ofstream seg_state_file(SIMULATION_SEG_STATE_NAME, std::ios::app);
      seg_state_file << nt+1 << " ";
      seg_state_file << std::scientific << std::setprecision(10);

      std::ofstream body_state_file(SIMULATION_BODY_STATE_NAME, std::ios::app);
      body_state_file << nt+1 << " ";
      body_state_file << std::scientific << std::setprecision(10);

      for (int n = 0; n < NSWIM; n++){

        swimmers[n].write_data(seg_state_file, body_state_file);

      }

      seg_state_file << std::endl;
      seg_state_file.close();

      body_state_file << std::endl;
      body_state_file.close();

      mobility.write_data(nt, swimmers); // Writes all velocity and force data.

      std::ofstream backup_file(SIMULATION_BACKUP_NAME);
      backup_file << nt+1 << " ";
      backup_file << std::scientific << std::setprecision(10);

      for (int n = 0; n < NSWIM; n++){

        swimmers[n].write_backup(backup_file);

      }

      backup_file << std::endl;
      backup_file.close();

      #if PRESCRIBED_CILIA

        std::ofstream fil_phase_file(SIMULATION_NAME+std::string("_filament_phases.dat"), std::ios::app);
        fil_phase_file << nt+1 << " ";
        fil_phase_file << std::scientific << std::setprecision(10);

        for (int n = 0; n < NSWIM; n++){

          for (int m = 0; m < NFIL; m++){

            fil_phase_file << swimmers[n].filaments[m].phase << " ";

          }

        }

        fil_phase_file << std::endl;
        fil_phase_file.close();

      #endif

    }

  }

  #if WRITE_GENERALISED_FORCES

    generalised_force_file << std::endl;
    generalised_force_file.close();

  #endif

  return 0;

};
