// config.hpp

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Include guard
#ifndef MY_CONFIG_HEADER_INCLUDED
#define MY_CONFIG_HEADER_INCLUDED

#define SIMULATION_NAME "junk"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Simulation type
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define CILIA_TYPE 4
// Valid options:
// 0 = Instability-driven cilia
// 1 = Geometrically-switching cilia (partially implemented)
// 2 = Constant base rotation (partially implemented)
// 3 = Cilia follow a prescribed sequence of shapes. This choice has some sub-types (see below).
// 4 = Squirmer-type simulation; i.e. there aren't actually any filaments/cilia. The slip velocity can be set in the mobility solver.

#if CILIA_TYPE==3

  #define SHAPE_SEQUENCE 2
  // Valid options:
  // 0 = 'Build-a-beat'. This choice has some parameters to set (see below).
  // 1 = The 'Fulford and Blake' beat pattern for mammalian airway cilia. See the data-fitting description in  "A model for the micro-structure in ciliated organisms", Blake (1972).
  // 2 = Coral larvae beat pattern. Data fitting done by me, in the same way as Blake (1972).

  #if SHAPE_SEQUENCE==0

    #define SCALED_BEAT_AMPLITUDE 2.0 // A value in (0,2), giving the beat amplitude in units of filament length.
    #define RECOVERY_STROKE_WINDOW_LENGTH 0.25 // A value in (0,1), giving the fraction of the beat cycle over which a given point on the filament completes its recovery-stroke tangent angle change.
    #define EFFECTIVE_STROKE_LENGTH 0.25 // A value in (0,1), giving the fraction of the cycle spent in the effective stroke.
    // N.B. We must have RECOVERY_STROKE_WINDOW_LENGTH + EFFECTIVE_STROKE_LENGTH < 1
    #define ZERO_VELOCITY_AVOIDANCE_LENGTH 0.05 // A value in (0,1), giving the maximum fraction of the cycle by which we shift the tangent angle curve to ensure the velocity cannot be zero everywhere along the filament at once.

  #endif

  #define DYNAMIC_PHASE_EVOLUTION true
  // If true, cilia phase speeds are solved for as part of the dynamics.
  // If false, phase_speed = omega0 is constant for each cilium.

  #define WRITE_GENERALISED_FORCES false
  // If true, this simulation will save its generalised forces to file for use as the reference values.
  // NOTE: This will overwrite any existing force reference file unless its name has been changed.

  #define ALLOW_BASE_PIVOTING false // NOT YET IMPLEMENTED!

#endif

#define BODY_OR_SURFACE_TYPE 2
// Valid options:
// 0 = An infinite plane wall at z = 0. This choice has some sub-types (see below).
// 1 = Deformed planes with 2 principal curvatures (partially implemented)
// 2 = Spherical bodies. This choice has some sub-types (see below).
// 3 = Toroidal bodies (partially implemented)

#if BODY_OR_SURFACE_TYPE==0

  #define SEEDING_TYPE 0
  // Valid options:
  // 0 = Filaments are placed on a rectangular grid.
  // 1 = Filaments are placed on a hexagonal grid.

#elif BODY_OR_SURFACE_TYPE==2

  #define SEEDING_TYPE 1
  // Valid options:
  // 0 = Filaments are evenly distributed over the surface.
  // 1 = Filaments are seeded in an equatorial band.
  // 2 = Platynereis-style seeding. Most filament are in an equatorial band but some form a small ring at the rear of the swimmer.

#endif

// Define whether the motion of the rigid bodies is imposed or allowed to evolve dynamically.
#define PRESCRIBED_BODY_VELOCITIES false

#define INITIAL_CONDITIONS_TYPE 0
// Valid options:
// 0 = Default
// 1 = Resume from backup file (only really valid if we ensure DT is the same etc., but this is NOT checked automatically)
// 2 = Fresh start from backup file (i.e. uses saved state but resets t=0. Use for hold-then-release style simulations)
#define INITIAL_CONDITIONS_FILE_NAME SIMULATION_NAME // SIMULATION_NAME or "a_different_sim_name"
// N.B. Simulations using GMRES can resume/start using backup files from Broyden-only simulations, but the opposite is not true.
// N.B. For options 0 and 2, whilst the simulation state will be fresh, all saved data will still be appended to any data from a previous simulation of the same name.

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Physical parameters
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define NFIL 1 // The number of filaments attached to the rigid body/surface in each swimmer.
#define NSEG 20 // The number of segments comprising each filament.
#define NSWIM 1 // The number of swimmers.
#define NBLOB 790 // The number of blobs to use as surface elements in each rigid body.

#define MU 1.0 // Fluid viscosity.

#define RSEG 1.0 // Segment radius.
#define RBLOB 1.4 // Surface blob radius.

#define KB 1800.0 // Bending modulus.
#define KT 1800.0 // Twist modulus.

#if CILIA_TYPE==0

  #define DIMENSIONLESS_FORCE 220.0

#elif CILIA_TYPE==1

  #define DRIVING_TORQUE_MAGNITUDE_RATIO 3.0 // How much stronger is the driving torque in the fast stroke than in the slow
  #define DEFAULT_BASE_TORQUE_MAGNITUDE (0.1*KB) // Driving torque magnitude in the fast stroke
  #define CRITICAL_ANGLE (0.4*PI) // The angle at which the cilia changes stroke

#elif CILIA_TYPE==2

  #define BASE_ROTATION_RATE 0.1

#endif

#if BODY_OR_SURFACE_TYPE==2 // Spherical

  #define R_OVER_L 1 // Aspect ratio of the body.

#endif

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Computational parameters
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Threads per block for kernel execution. Should be a multiple of 32, the warp size.
#define THREADS_PER_BLOCK 64

// This factor avoids over-adjusting during the Broyden's iterations.
// Setting it to 1 will give standard behaviour, and values smaller than 1 should help with convergence problems by having Broyden's method do more of the heavy lifting.
// Choosing 0.4 seems to work well when using GMRES, whilst 0.1 appears to be a good choice for a "Broyden's only" simulation.
#define JACOBIAN_CONFIDENCE_FACTOR 0.1

#define MAX_BROYDEN_ITER 100 // Maximum number of Broyden's method iterations per time-step.
#define TOL 1e-4 // Tolerance to be reached by the Broyden's method solve.

#define SOLVER_TYPE 1
// Valid options:
// 0: Use Broyden's method for absolutely everything. When there is a rigid body with forces (and velocities if they're not prescribed) to solve for,
//    the associated linear system is embedded in the wider Broyden's solve, rather than being solved for the current iterate at each iteration.
// 1: Use GMRES to solve the linear system at each iteration of Broyden's method.

#if SOLVER_TYPE==1

  #define MAX_LINEAR_SYSTEM_ITER 150 // Maximum number of iterations used to solve the linear system in each mobility solve.
  #define LINEAR_SYSTEM_TOL 1e-10 // Relative tolerance in the linear system solves.

  // GMRES preconditioner type.
  // Uses left preconditioning if set to false; if you don't want a preconditioner,
  // simply write "return in;" at the start of rpy_mobility_solver::apply_preconditioner(...) and it should be optimised away.
  // Right preconditioning seems to result in fewer GMRES iterations and also means that the error in GMRES
  // is the same as the error in the original system we want to solve, so this is the default option.
  // This ought to be checked whenever the preconditioner is changed though.
  #define USE_RIGHT_PRECON true

#endif

#define TOTAL_TIME_STEPS 20//(1*STEPS_PER_PERIOD) // Total number of time-steps in the simulation.
#define NUM_EULER_STEPS 1 // Number of time-steps to use backwards-Euler before switching to BDF2.

#if CILIA_TYPE==1

  #define TRANSITION_STEPS 21 // The number of time-steps over which to change the driving torque between strokes, must be >= 1
  #define DT 2.0
  #define PLOT_FREQUENCY_IN_STEPS 10

#else

  #define STEPS_PER_PERIOD 300
  #define SAVES_PER_PERIOD 20

#endif

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Derived/redefined parameters (these should be left alone)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define GEOMETRIC_CILIA (CILIA_TYPE==1)
#define INSTABILITY_CILIA (CILIA_TYPE==0)
#define CONSTANT_BASE_ROTATION (CILIA_TYPE==2)
#define PRESCRIBED_CILIA (CILIA_TYPE==3)
#define NO_CILIA_SQUIRMER (CILIA_TYPE==4)

#define INFINITE_PLANE_WALL (BODY_OR_SURFACE_TYPE==0)
#define SADDLE_BODIES (BODY_OR_SURFACE_TYPE==1)
#define SPHEROID_BODIES (BODY_OR_SURFACE_TYPE==2)
#define TORUS_BODIES (BODY_OR_SURFACE_TYPE==3)

#define USE_BROYDEN_FOR_EVERYTHING (SOLVER_TYPE==0)
#define USE_GMRES_FOR_LINEAR_SYSTEM (SOLVER_TYPE==1)

#define READ_INITIAL_CONDITIONS_FROM_BACKUP (INITIAL_CONDITIONS_TYPE != 0)

#define DL (2.2*RSEG) // Inter-segment distance.

#define PI 3.14159265358979323846264338327950288

#define SIMULATION_CONFIG_NAME SIMULATION_NAME ".par"
#define SIMULATION_BACKUP_NAME SIMULATION_NAME ".backup"
#define SIMULATION_BODY_STATE_NAME SIMULATION_NAME "_body_states.dat" // Blob states are recoverable from body states.
#define SIMULATION_SEG_STATE_NAME SIMULATION_NAME "_seg_states.dat"
#define SIMULATION_BODY_VEL_NAME SIMULATION_NAME "_body_vels.dat" // Blob velocities are recoverable from body velocities.
#define SIMULATION_SEG_VEL_NAME SIMULATION_NAME "_seg_vels.dat"
#define SIMULATION_BLOB_FORCES_NAME SIMULATION_NAME "_blob_forces.dat" // Body forces are recoverable from blob forces.
#define SIMULATION_SEG_FORCES_NAME SIMULATION_NAME "_seg_forces.dat"

#define DELETE_CURRENT_LINE "                                                                                                               " << "\r"

#define NTOTAL (NSWIM*(NFIL*NSEG + NBLOB))

#if USE_BROYDEN_FOR_EVERYTHING

  #define NBROY (3*NSWIM*(NBLOB + 2*(NFIL*NSEG + 1)))

#else

  #define NBROY (6*(NSWIM*NFIL*NSEG + NSWIM))

#endif

#if !GEOMETRIC_CILIA

  #define PLOT_FREQUENCY_IN_STEPS (STEPS_PER_PERIOD/SAVES_PER_PERIOD)

#endif

#if INSTABILITY_CILIA

  #define END_FORCE_MAGNITUDE (DIMENSIONLESS_FORCE*KB/(DL*DL*NSEG*NSEG))
  #define REPULSIVE_FORCE_FACTOR 2.0 // How much stronger is the barrier force than the driving force.
  #define DT (36.3833/STEPS_PER_PERIOD) // Based on the period of a single DIMENSIONLESS_FORCE = 220.0 filament above a no-slip wall.

#elif CONSTANT_BASE_ROTATION

  #define DT (2.0*PI/(BASE_ROTATION_RATE*STEPS_PER_PERIOD))
  #define REPULSIVE_FORCE_FACTOR 1.0
  #define END_FORCE_MAGNITUDE (0.1*KB) // There isn't really an end force, it's just used to define the repulsion.

#elif PRESCRIBED_CILIA

  #define DT (1.0/STEPS_PER_PERIOD) // Pick T = 1.

  #if USE_BROYDEN_FOR_EVERYTHING

    #error "The fully-prescribed cilia motion option currently doesn't support using Broyden's method to solve the system."

  #endif

  #if WRITE_GENERALISED_FORCES

    #define DYNAMIC_PHASE_EVOLUTION false
    #define NSWIM 1
    #define NFIL 1
    #define TOTAL_TIME_STEPS STEPS_PER_PERIOD

  #endif

  #define BUILD_A_BEAT (SHAPE_SEQUENCE==0)

  #define FIT_TO_DATA_BEAT (SHAPE_SEQUENCE != 0)
  #define FULFORD_AND_BLAKE_BEAT (SHAPE_SEQUENCE==1)
  #define CORAL_LARVAE_BEAT (SHAPE_SEQUENCE==2)

#endif

#if INFINITE_PLANE_WALL

  #define NSWIM 1
  #define NBLOB 0
  #define NTOTAL (NFIL*NSEG)
  #define NBROY (6*NFIL*NSEG)
  #define PRESCRIBED_BODY_VELOCITIES true
  #define RSEG 1.0
  #define MU 1.0

  #define RECTANGULAR_SEEDING (SEEDING_TYPE==0)
  #define HEXAGONAL_SEEDING (SEEDING_TYPE==1)

#endif

#if NO_CILIA_SQUIRMER

  #define NFIL 0
  #define NTOTAL (NSWIM*NBLOB)
  #define DT (1.0/STEPS_PER_PERIOD) // Pick T = 1.

#endif

#if SPHEROID_BODIES

  #define UNIFORM_SEEDING (SEEDING_TYPE==0)
  #define EQUATORIAL_SEEDING (SEEDING_TYPE==1)
  #define PLATY_SEEDING (SEEDING_TYPE==2)

#endif

#endif // MY_CONFIG_HEADER_INCLUDED
