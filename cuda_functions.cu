// cuda_functions.cu

#include "cuda_functions.hpp"
#include "config.hpp"

// RPY kernels

#if INFINITE_PLANE_WALL

  __global__ void Mss_mult(double * __restrict__ V, const double *const __restrict__ F, const double *const __restrict__ X, const int start_seg, const int num_segs){

    // Calculates the velocities of filament segments, given the forces and torques
    // on the segments, in a domain bounded by an infinite no-slip wall at z = 0.
    // It does so using the expressions given by Swan and Brady (2007).
    // NOTE: There seems to be a typo/sign error/unusual convention in the alternating
    // tensors that appear in the appendices of Swan and Brady's work. To get the correct
    // behaviour, it suffices to change the sign of all rot-trans (and hence trans-rot) terms.
    // This is highlighted in the code by evaluating these terms as they appear in the paper
    // but using -= rather than += to evaluate their contributions (except for the self-mobility
    // corrections, which are given on one line and have the signs changed in the expression itself).

    // N.B. This function assumes that RSEG = MU = 1!

    const double ForceFactor = 0.053051647697298; // = 1/(6*PI)
    const double TorqueFactor = 0.039788735772974; // = 1/(8*PI)

    const int index = threadIdx.x + blockIdx.x*blockDim.x;
    const int stride = blockDim.x*gridDim.x;

    // Declare the shared memory for this thread block
    __shared__ double x_shared[THREADS_PER_BLOCK];
    __shared__ double y_shared[THREADS_PER_BLOCK];
    __shared__ double z_shared[THREADS_PER_BLOCK];
    __shared__ double fx_shared[THREADS_PER_BLOCK];
    __shared__ double fy_shared[THREADS_PER_BLOCK];
    __shared__ double fz_shared[THREADS_PER_BLOCK];
    __shared__ double taux_shared[THREADS_PER_BLOCK];
    __shared__ double tauy_shared[THREADS_PER_BLOCK];
    __shared__ double tauz_shared[THREADS_PER_BLOCK];

    double vx, vy, vz, wx, wy, wz;
    double xi, yi, zi;

    // Stay in the loop as long as any thread in the block still needs to compute velocities.
    for (int i = (start_seg + index); (i-threadIdx.x) < (start_seg + num_segs); i+=stride){

      if (i < (start_seg + num_segs)){

        vx = 0.0; vy = 0.0; vz = 0.0; wx = 0.0; wy = 0.0; wz = 0.0;
        xi = X[3*i]; yi = X[3*i + 1]; zi = X[3*i + 2];

      }

      for (int j_start = 0; j_start < NFIL*NSEG; j_start += THREADS_PER_BLOCK){

        const int j_to_read = j_start + threadIdx.x;

        if (j_to_read < NFIL*NSEG){

          x_shared[threadIdx.x] = X[3*j_to_read];
          y_shared[threadIdx.x] = X[3*j_to_read + 1];
          z_shared[threadIdx.x] = X[3*j_to_read + 2];
          fx_shared[threadIdx.x] = F[6*j_to_read];
          fy_shared[threadIdx.x] = F[6*j_to_read + 1];
          fz_shared[threadIdx.x] = F[6*j_to_read + 2];
          taux_shared[threadIdx.x] = F[6*j_to_read + 3];
          tauy_shared[threadIdx.x] = F[6*j_to_read + 4];
          tauz_shared[threadIdx.x] = F[6*j_to_read + 5];

        }

        __syncthreads();

        if (i < (start_seg + num_segs)){

          for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < NFIL*NSEG); j++){

            if (i == (j + j_start)){

              // Self-mobility
              const double hm1 = 1.0/zi;
              const double hm3 = hm1*hm1*hm1;
              const double hm4 = hm3*hm1;
              const double hm5 = hm4*hm1;

              double a0 = ForceFactor - (9.0*hm1 - 2.0*hm3 + hm5)*0.003315727981081; // 0.003315727981081 = 1/(6*PI*16)
              double a1 = ForceFactor - (9.0*hm1 - 4.0*hm3 + hm5)*0.006631455962162; // 0.006631455962162 = 1/(6*PI*8)

              #if PRESCRIBED_CILIA

                vx += a0*fx_shared[j];
                vy += a0*fy_shared[j];
                vz += a1*fz_shared[j];

              #else

                const double a2 = 0.004973591971622*hm4; // 0.004973591971622 = 3/(6*PI*32)

                vx += a0*fx_shared[j] + a2*tauy_shared[j];
                vy += a0*fy_shared[j] - a2*taux_shared[j];
                vz += a1*fz_shared[j];

                a0 = TorqueFactor - 0.012433979929054*hm3; // 0.012433979929054 = 15/(6*PI*64)
                a1 = TorqueFactor - 0.004973591971622*hm3; // 0.004973591971622 = 3/(6*PI*32)

                wx += a0*taux_shared[j] - a2*fy_shared[j];
                wy += a0*tauy_shared[j] + a2*fx_shared[j];
                wz += a1*tauz_shared[j];

              #endif

            } else {

              double rm1 = 1.0/sqrt((xi - x_shared[j])*(xi - x_shared[j]) + (yi - y_shared[j])*(yi - y_shared[j]) + (zi - z_shared[j])*(zi - z_shared[j]));

              // Begin with the unbounded terms
              double rhatx = (xi - x_shared[j])*rm1;
              double rhaty = (yi - y_shared[j])*rm1;
              double rhatz = (zi - z_shared[j])*rm1;
              double rhat_dot_force = rhatx*fx_shared[j] + rhaty*fy_shared[j] + rhatz*fz_shared[j];
              double rhat_dot_torque = rhatx*taux_shared[j] + rhaty*tauy_shared[j] + rhatz*tauz_shared[j];

              double a0 = 0.039788735772974*rm1; // 0.039788735772974 = 1/(8*PI)
              double a1 = 1.0 + 0.666666666666667*rm1*rm1;
              double a2 = 1.0 - 2.0*rm1*rm1;

              vx += a0*(a1*fx_shared[j] + a2*rhat_dot_force*rhatx);
              vy += a0*(a1*fy_shared[j] + a2*rhat_dot_force*rhaty);
              vz += a0*(a1*fz_shared[j] + a2*rhat_dot_force*rhatz);

              #if !PRESCRIBED_CILIA

                a0 = 0.019894367886487*rm1*rm1*rm1; // 0.019894367886487 = 1/(16*PI)

                wx += a0*(3.0*rhat_dot_torque*rhatx - taux_shared[j]);
                wy += a0*(3.0*rhat_dot_torque*rhaty - tauy_shared[j]);
                wz += a0*(3.0*rhat_dot_torque*rhatz - tauz_shared[j]);

                a0 = 0.039788735772974*rm1*rm1; // 0.039788735772974 = 1/(8*PI)

                vx += a0*(tauy_shared[j]*rhatz - tauz_shared[j]*rhaty);
                vy += a0*(tauz_shared[j]*rhatx - taux_shared[j]*rhatz);
                vz += a0*(taux_shared[j]*rhaty - tauy_shared[j]*rhatx);

                wx += a0*(fy_shared[j]*rhatz - fz_shared[j]*rhaty);
                wy += a0*(fz_shared[j]*rhatx - fx_shared[j]*rhatz);
                wz += a0*(fx_shared[j]*rhaty - fy_shared[j]*rhatx);

              #endif

              // Now for the wall-induced corrections
              rm1 = 1.0/sqrt((xi - x_shared[j])*(xi - x_shared[j]) + (yi - y_shared[j])*(yi - y_shared[j]) + (zi + z_shared[j])*(zi + z_shared[j]));
              double rm2 = rm1*rm1;
              double rm3 = rm2*rm1;
              double rm4 = rm3*rm1;
              double rm5 = rm4*rm1;

              double ex = (xi - x_shared[j])*rm1;
              double ey = (yi - y_shared[j])*rm1;
              double ez = (zi + z_shared[j])*rm1;

              double h = z_shared[j]/(zi + z_shared[j]);

              double e_dot_force = ex*fx_shared[j] + ey*fy_shared[j] + ez*fz_shared[j];
              double e_dot_torque = ex*taux_shared[j] + ey*tauy_shared[j] + ez*tauz_shared[j];

              double ez2 = ez*ez;

              a0 = -(0.75*rm1*(1.0 + 2.0*h*(1.0-h)*ez2) + 0.5*(rm3*(1.0 - 3.0*ez2) - rm5*(1.0 - 5.0*ez2)));
              a1 = -(0.75*rm1*(1.0 - 6.0*h*(1.0-h)*ez2) - 1.5*rm3*(1.0 - 5.0*ez2) + 2.5*rm5*(1.0 - 7.0*ez2))*e_dot_force;
              a2 = ez*fz_shared[j]*(1.5*h*rm1*(1.0 - 6.0*(1.0-h)*ez2) - 3.0*rm3*(1.0 - 5.0*ez2) + 5.0*rm5*(2.0 - 7.0*ez2));
              double a3 = (1.5*h*rm1 - 5.0*rm5)*e_dot_force;
              double a4 = -((3.0*h*h*rm1 + 3.0*rm3)*ez2 + (2.0 - 15.0*ez2)*rm5)*fz_shared[j];

              vx += ForceFactor*(a0*fx_shared[j] + (a1 + a2)*ex);
              vy += ForceFactor*(a0*fy_shared[j] + (a1 + a2)*ey);
              vz += ForceFactor*(a0*fz_shared[j] + (a1 + a2 + a3)*ez + a4);

              #if !PRESCRIBED_CILIA

                a0 = rm3*(0.375 - 2.25*ez2);
                a1 = -1.125*rm3*e_dot_torque;
                a2 = 2.25*ez*rm3*e_dot_torque;
                a3 = 2.25*rm3*(ex*tauy_shared[j] - ey*taux_shared[j]);

                wx += ForceFactor*(a0*taux_shared[j] + a1*ex - a3*ey);
                wy += ForceFactor*(a0*tauy_shared[j] + a1*ey + a3*ex);
                wz += ForceFactor*(a0*tauz_shared[j] + a1*ez + a2);

                a0 = 0.75*rm2;
                a1 = fz_shared[j]*(9.0*h*ez2*rm2 + (1.5 - 15.0*ez2)*rm4);
                a2 = -ez*e_dot_force*(4.5*h*rm2 - 7.5*rm4);
                a3 = -1.5*ez*(h*rm2 - rm4);

                wx -= ForceFactor*(a0*(fy_shared[j]*ez - fz_shared[j]*ey) - (a1+a2)*ey + a3*fy_shared[j]);
                wy -= ForceFactor*(a0*(fz_shared[j]*ex - fx_shared[j]*ez) + (a1+a2)*ex - a3*fx_shared[j]);
                wz -= ForceFactor*a0*(fx_shared[j]*ey - fy_shared[j]*ex);

                // We use symmetry of the grand mobility matrix to calculate this contribution as
                // the torque multiplied by the transpose of the sub-block relating the angular velocity
                // of particle j to the force on particle i.

                h = 1.0 - h;
                ex *= -1.0;
                ey *= -1.0;

                const double e_cross_T_3 = ex*tauy_shared[j] - ey*taux_shared[j];
                a1 = (9.0*h*ez2*rm2 + (1.5 - 15.0*ez2)*rm4)*e_cross_T_3;
                a2 = -ez*(4.5*h*rm2 - 7.5*rm4)*e_cross_T_3;
                a3 = -1.5*ez*(h*rm2 - rm4);

                vx -= ForceFactor*(a0*(tauz_shared[j]*ey - tauy_shared[j]*ez) + a2*ex - a3*tauy_shared[j]);
                vy -= ForceFactor*(a0*(taux_shared[j]*ez - tauz_shared[j]*ex) + a2*ey + a3*taux_shared[j]);
                vz -= ForceFactor*(a0*(tauy_shared[j]*ex - taux_shared[j]*ey) + a2*ez + a1);

              #endif

            }
        }

      }

        __syncthreads();

      } // End of loop over filament segment forces and torques.

      if (i < (start_seg + num_segs)){

        const int p = 6*(i - start_seg);

        V[p] = vx;
        V[p + 1] = vy;
        V[p + 2] = vz;
        V[p + 3] = wx;
        V[p + 4] = wy;
        V[p + 5] = wz;

      }

    } // End of striding loop over filament segment velocities.

  } // End of Mss_mult kernel.

#else

  __global__ void Mss_mult(double * __restrict__ V, const double *const __restrict__ F, const double *const __restrict__ X, const int start_seg, const int num_segs){

    // Calculates the velocities of filament segments given the forces and torques
    // on the segments.

    const int index = threadIdx.x + blockIdx.x*blockDim.x;
    const int stride = blockDim.x*gridDim.x;

    // Declare the shared memory for this thread block
    __shared__ double x_shared[THREADS_PER_BLOCK];
    __shared__ double y_shared[THREADS_PER_BLOCK];
    __shared__ double z_shared[THREADS_PER_BLOCK];
    __shared__ double fx_shared[THREADS_PER_BLOCK];
    __shared__ double fy_shared[THREADS_PER_BLOCK];
    __shared__ double fz_shared[THREADS_PER_BLOCK];
    __shared__ double taux_shared[THREADS_PER_BLOCK];
    __shared__ double tauy_shared[THREADS_PER_BLOCK];
    __shared__ double tauz_shared[THREADS_PER_BLOCK];

    double vx, vy, vz, wx, wy, wz;
    double xi, yi, zi;

    // Stay in the loop as long as any thread in the block still needs to compute velocities.
    for (int i = (start_seg + index); (i-threadIdx.x) < (start_seg + num_segs); i+=stride){

      if (i < (start_seg + num_segs)){

        vx = 0.0; vy = 0.0; vz = 0.0; wx = 0.0; wy = 0.0; wz = 0.0;
        xi = X[3*i]; yi = X[3*i + 1]; zi = X[3*i + 2];

      }

      for (int j_start = 0; j_start < NSWIM*NFIL*NSEG; j_start += THREADS_PER_BLOCK){

        const int j_to_read = j_start + threadIdx.x;

        if (j_to_read < NSWIM*NFIL*NSEG){

          x_shared[threadIdx.x] = X[3*j_to_read];
          y_shared[threadIdx.x] = X[3*j_to_read + 1];
          z_shared[threadIdx.x] = X[3*j_to_read + 2];
          fx_shared[threadIdx.x] = F[6*j_to_read];
          fy_shared[threadIdx.x] = F[6*j_to_read + 1];
          fz_shared[threadIdx.x] = F[6*j_to_read + 2];
          taux_shared[threadIdx.x] = F[6*j_to_read + 3];
          tauy_shared[threadIdx.x] = F[6*j_to_read + 4];
          tauz_shared[threadIdx.x] = F[6*j_to_read + 5];

        }

        __syncthreads();

        if (i < (start_seg + num_segs)){

          for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < NSWIM*NFIL*NSEG); j++){

            double rx = xi - x_shared[j];
            double ry = yi - y_shared[j];
            double rz = zi - z_shared[j];

            const double r = sqrt(rx*rx + ry*ry + rz*rz);
            const double rm1 = 1.0/(r + 1e-20);

            rx *= rm1; ry *= rm1; rz *= rm1;

            double A, B;

            // translation-translation

            double r_dot = rx*fx_shared[j] + ry*fy_shared[j] + rz*fz_shared[j];

            if (r > 2.0*RSEG){

              double temp = 2.0*RSEG*RSEG*rm1*rm1;

              A = 1.0 + temp/3.0;
              B = 1.0 - temp;

              temp = rm1/(8.0*PI*MU);

              A *= temp;
              B *= temp*r_dot;

            } else {

              B = 3.0*r/(32.0*RSEG);
              A = 1.0 - 3.0*B;

              double temp = 1.0/(6.0*PI*MU*RSEG);

              A *= temp;
              B *= temp*r_dot;

            }

            vx += A*fx_shared[j] + B*rx;
            vy += A*fy_shared[j] + B*ry;
            vz += A*fz_shared[j] + B*rz;

            #if !PRESCRIBED_CILIA

              // rotation-rotation

              r_dot = rx*taux_shared[j] + ry*tauy_shared[j] + rz*tauz_shared[j];

              if (r > 2.0*RSEG){

                A = -rm1*rm1*rm1/(16.0*PI*MU);
                B = -3.0*r_dot*A;

              } else {

                double temp = r/RSEG;

                B = 9.0*temp/32.0;
                A = -3.0*B;

                temp *= temp*temp;

                A += 1.0 + 5.0*temp/64.0;
                B -= 3.0*temp/64.0;

                temp = 1.0/(8.0*PI*MU*RSEG*RSEG*RSEG);

                A *= temp;
                B *= temp*r_dot;

              }

              wx += A*taux_shared[j] + B*rx;
              wy += A*tauy_shared[j] + B*ry;
              wz += A*tauz_shared[j] + B*rz;

              // translation-rotation

              if (r > 2.0*RSEG){

                A = rm1*rm1/(8.0*PI*MU);

              } else {

                A = r*(1.0 - 3.0*r/(8.0*RSEG));
                A /= 16.0*PI*MU*RSEG*RSEG*RSEG;

              }

              wx += A*(fy_shared[j]*rz - fz_shared[j]*ry);
              wy += A*(fz_shared[j]*rx - fx_shared[j]*rz);
              wz += A*(fx_shared[j]*ry - fy_shared[j]*rx);

              vx += A*(tauy_shared[j]*rz - tauz_shared[j]*ry);
              vy += A*(tauz_shared[j]*rx - taux_shared[j]*rz);
              vz += A*(taux_shared[j]*ry - tauy_shared[j]*rx);

            #endif

          }

        }

        __syncthreads();

      } // End of loop over filament segment forces and torques.

      if (i < (start_seg + num_segs)){

        const int p = 6*(i - start_seg);

        V[p] = vx;
        V[p + 1] = vy;
        V[p + 2] = vz;
        V[p + 3] = wx;
        V[p + 4] = wy;
        V[p + 5] = wz;

      }

    } // End of striding loop over filament segment velocities.

  } // End of Mss_mult kernel.

#endif

__global__ void Mbb_mult(double * __restrict__ V, const double *const __restrict__ F, const double *const __restrict__ X, const int start_blob, const int num_blobs){

  // Calculates the velocities of rigid-body blobs given the forces they experience.

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  // Declare the shared memory for this thread block
  __shared__ double x_shared[THREADS_PER_BLOCK];
  __shared__ double y_shared[THREADS_PER_BLOCK];
  __shared__ double z_shared[THREADS_PER_BLOCK];
  __shared__ double fx_shared[THREADS_PER_BLOCK];
  __shared__ double fy_shared[THREADS_PER_BLOCK];
  __shared__ double fz_shared[THREADS_PER_BLOCK];

  double vx, vy, vz;
  double xi, yi, zi;

  // Stay in the loop as long as any thread in the block still needs to compute velocities.
  for (int i = (start_blob + index); (i-threadIdx.x) < (start_blob + num_blobs); i+=stride){

    if (i < (start_blob + num_blobs)){

      vx = 0.0; vy = 0.0; vz = 0.0;
      xi = X[3*i]; yi = X[3*i + 1]; zi = X[3*i + 2];

    }

    for (int j_start = 0; j_start < NSWIM*NBLOB; j_start += THREADS_PER_BLOCK){

      const int j_to_read = j_start + threadIdx.x;

      if (j_to_read < NSWIM*NBLOB){

        x_shared[threadIdx.x] = X[3*j_to_read];
        y_shared[threadIdx.x] = X[3*j_to_read + 1];
        z_shared[threadIdx.x] = X[3*j_to_read + 2];
        fx_shared[threadIdx.x] = F[3*j_to_read];
        fy_shared[threadIdx.x] = F[3*j_to_read + 1];
        fz_shared[threadIdx.x] = F[3*j_to_read + 2];

      }

      __syncthreads();

      if (i < (start_blob + num_blobs)){

        for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < NSWIM*NBLOB); j++){

          double rx = xi - x_shared[j];
          double ry = yi - y_shared[j];
          double rz = zi - z_shared[j];

          const double r = sqrt(rx*rx + ry*ry + rz*rz);
          const double rm1 = 1.0/(r + 1e-20);

          rx *= rm1; ry *= rm1; rz *= rm1;

          double A, B;

          const double r_dot = rx*fx_shared[j] + ry*fy_shared[j] + rz*fz_shared[j];

          if (r > 2.0*RBLOB){

            double temp = 2.0*RBLOB*RBLOB*rm1*rm1;

            A = 1.0 + temp/3.0;
            B = 1.0 - temp;

            temp = rm1/(8.0*PI*MU);

            A *= temp;
            B *= temp*r_dot;

          } else {

            B = 3.0*r/(32.0*RBLOB);
            A = 1.0 - 3.0*B;

            double temp = 1.0/(6.0*PI*MU*RBLOB);

            A *= temp;
            B *= temp*r_dot;

          }

          vx += A*fx_shared[j] + B*rx;
          vy += A*fy_shared[j] + B*ry;
          vz += A*fz_shared[j] + B*rz;

        }

      }

      __syncthreads();

    } // End of loop over blob forces.

    if (i < (start_blob + num_blobs)){

      const int p = 3*(i - start_blob);

      V[p] = vx;
      V[p + 1] = vy;
      V[p + 2] = vz;

    }

  } // End of striding loop over blob velocities.

} // End of Mbb_mult kernel.

__global__ void Msb_mult(double * __restrict__ V, const double *const __restrict__ F, const double *const __restrict__ Xs, const double *const __restrict__ Xb, const int start_seg, const int num_segs){

  // Calculates the velocities of filament segments given the forces on
  // rigid-body blobs.

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  // Declare the shared memory for this thread block
  __shared__ double x_shared[THREADS_PER_BLOCK];
  __shared__ double y_shared[THREADS_PER_BLOCK];
  __shared__ double z_shared[THREADS_PER_BLOCK];
  __shared__ double fx_shared[THREADS_PER_BLOCK];
  __shared__ double fy_shared[THREADS_PER_BLOCK];
  __shared__ double fz_shared[THREADS_PER_BLOCK];

  double vx, vy, vz, wx, wy, wz;
  double xi, yi, zi;

  const double RMAX_MINUS_RMIN = abs(RSEG - RBLOB);

  // Stay in the loop as long as any thread in the block still needs to compute velocities.
  for (int i = (start_seg + index); (i-threadIdx.x) < (start_seg + num_segs); i+=stride){

    if (i < (start_seg + num_segs)){

      vx = 0.0; vy = 0.0; vz = 0.0; wx = 0.0; wy = 0.0; wz = 0.0;
      xi = Xs[3*i]; yi = Xs[3*i + 1]; zi = Xs[3*i + 2];

    }

    for (int j_start = 0; j_start < NSWIM*NBLOB; j_start += THREADS_PER_BLOCK){

      const int j_to_read = j_start + threadIdx.x;

      if (j_to_read < NSWIM*NBLOB){

        x_shared[threadIdx.x] = Xb[3*j_to_read];
        y_shared[threadIdx.x] = Xb[3*j_to_read + 1];
        z_shared[threadIdx.x] = Xb[3*j_to_read + 2];
        fx_shared[threadIdx.x] = F[3*j_to_read];
        fy_shared[threadIdx.x] = F[3*j_to_read + 1];
        fz_shared[threadIdx.x] = F[3*j_to_read + 2];

      }

      __syncthreads();

      if (i < (start_seg + num_segs)){

        for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < NSWIM*NBLOB); j++){

          double rx = xi - x_shared[j];
          double ry = yi - y_shared[j];
          double rz = zi - z_shared[j];

          const double r = sqrt(rx*rx + ry*ry + rz*rz);
          const double rm1 = 1.0/(r + 1e-20);

          rx *= rm1; ry *= rm1; rz *= rm1;

          double A, B;

          // translation-translation

          const double r_dot = rx*fx_shared[j] + ry*fy_shared[j] + rz*fz_shared[j];

          if (r > RSEG + RBLOB){

            A = rm1*rm1*(RSEG*RSEG + RBLOB*RBLOB)/3.0;
            B = 1.0 - 3.0*A;
            A += 1.0;

            double temp = rm1/(8.0*PI*MU);

            A *= temp;
            B *= temp*r_dot;

          } else if (r > RMAX_MINUS_RMIN){

            double temp = 32.0*r*r*r;

            A = RMAX_MINUS_RMIN*RMAX_MINUS_RMIN + 3.0*r*r;
            A *= -A/temp;
            A += 0.5*(RSEG + RBLOB);

            B = RMAX_MINUS_RMIN*RMAX_MINUS_RMIN - r*r;
            B *= 3.0*B/temp;

            temp = 1.0/(6.0*PI*MU*RSEG*RBLOB);

            A *= temp;
            B *= temp*r_dot;

          } else {

            A = 1.0/(6.0*PI*MU*((RSEG > RBLOB) ? RSEG : RBLOB)); // ternary operator gives largest radius
            B = 0.0;

          }

          vx += A*fx_shared[j] + B*rx;
          vy += A*fy_shared[j] + B*ry;
          vz += A*fz_shared[j] + B*rz;

          #if !PRESCRIBED_CILIA

            // translation-rotation

            if (r > RSEG + RBLOB){

              A = rm1*rm1/(8.0*PI*MU);

            } else if (r > RMAX_MINUS_RMIN){

              A = RSEG - RBLOB + r;
              A *= rm1*rm1*A*(RBLOB*RBLOB + 2.0*RBLOB*(r + RSEG) - 3.0*(RSEG - r)*(RSEG - r))/(128.0*PI*MU*RBLOB*RSEG*RSEG*RSEG);

            } else {

              A = (RSEG > RBLOB)*r/(8.0*PI*MU*RSEG*RSEG*RSEG);

            }

            wx += A*(fy_shared[j]*rz - fz_shared[j]*ry);
            wy += A*(fz_shared[j]*rx - fx_shared[j]*rz);
            wz += A*(fx_shared[j]*ry - fy_shared[j]*rx);

          #endif

        }

      }

      __syncthreads();

    } // End of loop over blob forces.

    if (i < (start_seg + num_segs)){

      const int p = 6*(i - start_seg);

      V[p] += vx;
      V[p + 1] += vy;
      V[p + 2] += vz;
      V[p + 3] += wx;
      V[p + 4] += wy;
      V[p + 5] += wz;

    }

  } // End of striding loop over filament segment velocities.

} // End of Msb_mult kernel.


__global__ void Mbs_mult(double * __restrict__ V, const double *const __restrict__ F, const double *const __restrict__ Xb, const double *const __restrict__ Xs, const int start_blob, const int num_blobs){

  // Calculates the velocities of blobs given the forces and torques on
  // filament segments.

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  // Declare the shared memory for this thread block
  __shared__ double x_shared[THREADS_PER_BLOCK];
  __shared__ double y_shared[THREADS_PER_BLOCK];
  __shared__ double z_shared[THREADS_PER_BLOCK];
  __shared__ double fx_shared[THREADS_PER_BLOCK];
  __shared__ double fy_shared[THREADS_PER_BLOCK];
  __shared__ double fz_shared[THREADS_PER_BLOCK];
  __shared__ double taux_shared[THREADS_PER_BLOCK];
  __shared__ double tauy_shared[THREADS_PER_BLOCK];
  __shared__ double tauz_shared[THREADS_PER_BLOCK];

  double vx, vy, vz;
  double xi, yi, zi;

  const double RMAX_MINUS_RMIN = abs(RSEG - RBLOB);

  // Stay in the loop as long as any thread in the block still needs to compute velocities.
  for (int i = (start_blob + index); (i-threadIdx.x) < (start_blob + num_blobs); i+=stride){

    if (i < (start_blob + num_blobs)){

      vx = 0.0; vy = 0.0; vz = 0.0;
      xi = Xb[3*i]; yi = Xb[3*i + 1]; zi = Xb[3*i + 2];

    }

    for (int j_start = 0; j_start < NSWIM*NFIL*NSEG; j_start += THREADS_PER_BLOCK){

      const int j_to_read = j_start + threadIdx.x;

      if (j_to_read < NSWIM*NFIL*NSEG){

        x_shared[threadIdx.x] = Xs[3*j_to_read];
        y_shared[threadIdx.x] = Xs[3*j_to_read + 1];
        z_shared[threadIdx.x] = Xs[3*j_to_read + 2];
        fx_shared[threadIdx.x] = F[6*j_to_read];
        fy_shared[threadIdx.x] = F[6*j_to_read + 1];
        fz_shared[threadIdx.x] = F[6*j_to_read + 2];
        taux_shared[threadIdx.x] = F[6*j_to_read + 3];
        tauy_shared[threadIdx.x] = F[6*j_to_read + 4];
        tauz_shared[threadIdx.x] = F[6*j_to_read + 5];

      }

      __syncthreads();

      if (i < (start_blob + num_blobs)){

        for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < NSWIM*NFIL*NSEG); j++){

          double rx = xi - x_shared[j];
          double ry = yi - y_shared[j];
          double rz = zi - z_shared[j];

          const double r = sqrt(rx*rx + ry*ry + rz*rz);
          const double rm1 = 1.0/(r + 1e-20);

          rx *= rm1; ry *= rm1; rz *= rm1;

          double A, B;

          // translation-translation

          const double r_dot = rx*fx_shared[j] + ry*fy_shared[j] + rz*fz_shared[j];

          if (r > RSEG + RBLOB){

            A = rm1*rm1*(RSEG*RSEG + RBLOB*RBLOB)/3.0;
            B = 1.0 - 3.0*A;
            A += 1.0;

            double temp = rm1/(8.0*PI*MU);

            A *= temp;
            B *= temp*r_dot;

          } else if (r > RMAX_MINUS_RMIN){

            double temp = 32.0*r*r*r;

            A = RMAX_MINUS_RMIN*RMAX_MINUS_RMIN + 3.0*r*r;
            A *= -A/temp;
            A += 0.5*(RSEG + RBLOB);

            B = RMAX_MINUS_RMIN*RMAX_MINUS_RMIN - r*r;
            B *= 3.0*B/temp;

            temp = 1.0/(6.0*PI*MU*RSEG*RBLOB);

            A *= temp;
            B *= temp*r_dot;

          } else {

            A = 1.0/(6.0*PI*MU*((RSEG > RBLOB) ? RSEG : RBLOB)); // ternary operator gives largest radius
            B = 0.0;

          }

          vx += A*fx_shared[j] + B*rx;
          vy += A*fy_shared[j] + B*ry;
          vz += A*fz_shared[j] + B*rz;

          #if !PRESCRIBED_CILIA

            // translation-rotation

            if (r > RSEG + RBLOB){

              A = rm1*rm1/(8.0*PI*MU);

            } else if (r > RMAX_MINUS_RMIN){

              A = RSEG - RBLOB + r;
              A *= rm1*rm1*A*(RBLOB*RBLOB + 2.0*RBLOB*(r + RSEG) - 3.0*(RSEG - r)*(RSEG - r))/(128.0*PI*MU*RBLOB*RSEG*RSEG*RSEG);

            } else {

              A = (RSEG > RBLOB)*r/(8.0*PI*MU*RSEG*RSEG*RSEG);

            }

            vx += A*(tauy_shared[j]*rz - tauz_shared[j]*ry);
            vy += A*(tauz_shared[j]*rx - taux_shared[j]*rz);
            vz += A*(taux_shared[j]*ry - tauy_shared[j]*rx);

          #endif

        }

      }

      __syncthreads();

    } // End of loop over segment forces and torques.

    if (i < (start_blob + num_blobs)){

      const int p = 3*(i - start_blob);

      #if USE_BROYDEN_FOR_EVERYTHING

        V[p] += vx;
        V[p + 1] += vy;
        V[p + 2] += vz;

      #else

        V[p] = vx;
        V[p + 1] = vy;
        V[p + 2] = vz;

      #endif

    }

  } // End of striding loop over blob velocities.

} // End of Mbs_mult kernel.







// Stokes drag kernels

__global__ void Ms_mult(double * __restrict__ V, const double *const __restrict__ F, const int start_seg, const int num_segs){

    // Calculates the velocities of filament segments given the forces and torques
    // on the segments.

    const int index = threadIdx.x + blockIdx.x*blockDim.x;
    const int stride = blockDim.x*gridDim.x;

    const double temp = 1.0/(6.0*PI*MU*RSEG);
    const double temp2 = 1.0/(8.0*PI*MU*RSEG*RSEG*RSEG);

    // Stay in the loop as long as any thread in the block still needs to compute velocities.
    for (int i = (start_seg + index); i < (start_seg + num_segs); i+=stride){

      const int p = 6*(i - start_seg);

      V[p] = temp*F[6*i];
      V[p + 1] = temp*F[6*i + 1];
      V[p + 2] = temp*F[6*i + 2];
      V[p + 3] = temp2*F[6*i + 3];
      V[p + 4] = temp2*F[6*i + 4];
      V[p + 5] = temp2*F[6*i + 5];

    } // End of striding loop over filament segment velocities.

  } // End of Ms_mult kernel.

__global__ void Mb_mult(double * __restrict__ V, const double *const __restrict__ F, const int start_blob, const int num_blobs){

  // Calculates the velocities of rigid-body blobs given the forces they experience.

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  const double temp = 1.0/(6.0*PI*MU*RBLOB);

  // Stay in the loop as long as any thread in the block still needs to compute velocities.
  for (int i = (start_blob + index); i < (start_blob + num_blobs); i+=stride){

    const int p = 3*(i - start_blob);

    V[p] = temp*F[3*i];
    V[p + 1] = temp*F[3*i + 1];
    V[p + 2] = temp*F[3*i + 2];

  } // End of striding loop over blob velocities.

} // End of Mb_mult kernel.

__global__ void Mb_fill_zero(double * __restrict__ V, const int start_blob, const int num_blobs){

  // Fill zero velocity arrays

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  // Stay in the loop as long as any thread in the block still needs to fill zeros.
  for (int i = (start_blob + index); i < (start_blob + num_blobs); i+=stride){

    const int p = 3*(i - start_blob);

    #if USE_BROYDEN_FOR_EVERYTHING

      V[p] += 0.0;
      V[p + 1] += 0.0;
      V[p + 2] += 0.0;

    #else

      V[p] = 0.0;
      V[p + 1] = 0.0;
      V[p + 2] = 0.0;

    #endif

  } // End of striding loop over blob velocities.

} // End of Mb_fill_zero kernel.








// Generic interaction kernels

__global__ void barrier_forces(double * __restrict__ f_segs, double * __restrict__ f_blobs_repulsion, const double *const __restrict__ x_segs, const double *const __restrict__ x_blobs, const int start_seg, const int num_segs, const int start_blob, const int num_blobs){

  #if !(PRESCRIBED_CILIA || NO_CILIA_SQUIRMER)

    // Work out which particle(s) this thread will compute the force for
    const int index = threadIdx.x + blockIdx.x*blockDim.x;
    const int stride = blockDim.x*gridDim.x;

    // Declare the shared memory for this thread block
    __shared__ double x_shared[THREADS_PER_BLOCK];
    __shared__ double y_shared[THREADS_PER_BLOCK];
    __shared__ double z_shared[THREADS_PER_BLOCK];

    double fx, fy, fz;
    double xi, yi, zi;
    int fili;

    for (int i = (start_seg + index); (i-threadIdx.x) < (start_seg + num_segs); i+=stride){

      if (i < (start_seg + num_segs)){

        fx = 0.0;
        fy = 0.0;
        fz = 0.0;

        xi = x_segs[3*i];
        yi = x_segs[3*i + 1];
        zi = x_segs[3*i + 2];

        fili = i/NSEG;

        #if INFINITE_PLANE_WALL

          if (zi < 0.5*DL){ // "const double height_of_base_above_surface = 0.5*DL;" in rigid_body.cu

            fz = fmin(1.0, 1.0 - (zi - RSEG)/(0.5*DL - RSEG)); // max magnitude one radius from wall
            fz *= REPULSIVE_FORCE_FACTOR*END_FORCE_MAGNITUDE*fz*fz*fz; // 4th power

          }

        #endif

      }

      for (int j_start = 0; j_start < NTOTAL; j_start += THREADS_PER_BLOCK){

        const int j_to_read = j_start + threadIdx.x;

        if (j_to_read < NSWIM*NFIL*NSEG){

          x_shared[threadIdx.x] = x_segs[3*j_to_read];
          y_shared[threadIdx.x] = x_segs[3*j_to_read + 1];
          z_shared[threadIdx.x] = x_segs[3*j_to_read + 2];

        } else if (j_to_read < NTOTAL){

          x_shared[threadIdx.x] = x_blobs[3*(j_to_read - NSWIM*NFIL*NSEG)];
          y_shared[threadIdx.x] = x_blobs[3*(j_to_read - NSWIM*NFIL*NSEG) + 1];
          z_shared[threadIdx.x] = x_blobs[3*(j_to_read - NSWIM*NFIL*NSEG) + 2];

        }

        __syncthreads();

        if (i < (start_seg + num_segs)){

          for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < NTOTAL); j++){

            const double a_sum = (j_start + j < NSWIM*NFIL*NSEG) ? 2.0*RSEG : RSEG + RBLOB;
            const double chi_fac = 10.0/a_sum; // 1.0/(1.1*a_sum - a_sum)

            const double dx = xi - x_shared[j];
            const double dy = yi - y_shared[j];
            const double dz = zi - z_shared[j];

            const double dist = sqrt(dx*dx + dy*dy + dz*dz);

            int filj = (j_start + j)/NSEG;

            if (!(fili==filj && abs(i -(j_start + j))<=1) && (dist < 1.1*a_sum)){

              double fac = fmin(1.0, 1.0 - chi_fac*(dist - a_sum));
              fac *= REPULSIVE_FORCE_FACTOR*END_FORCE_MAGNITUDE*fac*fac*fac;

              const double dm1 = 1.0/dist;

              fx += fac*dx*dm1;
              fy += fac*dy*dm1;
              fz += fac*dz*dm1;

            }

          }

        }

        __syncthreads();

      }

      if (i < (start_seg + num_segs)){

        f_segs[6*i] += fx; // Don't shift index as f_segs has global size.
        f_segs[6*i + 1] += fy;
        f_segs[6*i + 2] += fz;

      }

    }

    #if !PRESCRIBED_BODY_VELOCITIES

      for (int i = (start_blob + index); (i-threadIdx.x) < (start_blob + num_blobs); i+=stride){

        if (i < (start_blob + num_blobs)){

          fx = 0.0;
          fy = 0.0;
          fz = 0.0;

          xi = x_blobs[3*i];
          yi = x_blobs[3*i + 1];
          zi = x_blobs[3*i + 2];

        }

        for (int j_start = 0; j_start < NTOTAL; j_start += THREADS_PER_BLOCK){

          const int j_to_read = j_start + threadIdx.x;

          if (j_to_read < NSWIM*NFIL*NSEG){

            x_shared[threadIdx.x] = x_segs[3*j_to_read];
            y_shared[threadIdx.x] = x_segs[3*j_to_read + 1];
            z_shared[threadIdx.x] = x_segs[3*j_to_read + 2];

          } else if (j_to_read < NTOTAL){

            x_shared[threadIdx.x] = x_blobs[3*(j_to_read - NSWIM*NFIL*NSEG)];
            y_shared[threadIdx.x] = x_blobs[3*(j_to_read - NSWIM*NFIL*NSEG) + 1];
            z_shared[threadIdx.x] = x_blobs[3*(j_to_read - NSWIM*NFIL*NSEG) + 2];

          }

          __syncthreads();

          if (i < (start_blob + num_blobs)){

            for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < NTOTAL); j++){

              const double a_sum = (j_start + j < NSWIM*NFIL*NSEG) ? RBLOB + RSEG : 2.0*RBLOB;
              const double chi_fac = 10.0/a_sum; // 1.0/(1.1*a_sum - a_sum)

              const double dx = xi - x_shared[j];
              const double dy = yi - y_shared[j];
              const double dz = zi - z_shared[j];

              const double dist = sqrt(dx*dx + dy*dy + dz*dz);

              if (((i + NSWIM*NFIL*NSEG) != (j_start + j)) && (dist < 1.1*a_sum)){

                double fac = fmin(1.0, 1.0 - chi_fac*(dist - a_sum));
                fac *= REPULSIVE_FORCE_FACTOR*END_FORCE_MAGNITUDE*fac*fac*fac;

                const double dm1 = 1.0/dist;

                fx += fac*dx*dm1;
                fy += fac*dy*dm1;
                fz += fac*dz*dm1;

              }

            }

          }

          __syncthreads();

        }

        if (i < (start_blob + num_blobs)){

          const int p = 3*(i - start_blob);

          f_blobs_repulsion[p] = fx;
          f_blobs_repulsion[p + 1] = fy;
          f_blobs_repulsion[p + 2] = fz;

        }

      }

    #endif

  #endif

} // End of barrier forces kernel.

__global__ void seeding_repulsion_step(double *F, const double *X, const int N){

  // Work out which particle(s) this thread will compute the force for
  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  // Declare the shared memory for this thread block
  __shared__ double x_shared[THREADS_PER_BLOCK];
  __shared__ double y_shared[THREADS_PER_BLOCK];
  __shared__ double z_shared[THREADS_PER_BLOCK];

  double fx, fy, fz;
  double xi, yi, zi;

  for (int i = index; (i-threadIdx.x) < N; i+=stride){

    if (i < N){

      fx = 0.0;
      fy = 0.0;
      fz = 0.0;

      xi = X[3*i];
      yi = X[3*i + 1];
      zi = X[3*i + 2];

      #if PRESCRIBED_CILIA

        if (N == NFIL){ // This does what we want (i.e. ensures we only apply this force when seeding filaments) unless this is a weird simulation with NFIL=NBLOB.

          fz -= 0.5*zi; // Wants to be at z = 0.

        }

      #endif

    }

    for (int j_start = 0; j_start < N; j_start += THREADS_PER_BLOCK){

      const int j_to_read = j_start + threadIdx.x;

      if (j_to_read < N){

        x_shared[threadIdx.x] = X[3*j_to_read];
        y_shared[threadIdx.x] = X[3*j_to_read + 1];
        z_shared[threadIdx.x] = X[3*j_to_read + 2];

      }

      __syncthreads();

      if (i < N){

        for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < N); j++){

          if (j_start + j != i){

            const double dx = xi - x_shared[j];
            const double dy = yi - y_shared[j];
            const double dz = zi - z_shared[j];

            const double dist = sqrt(dx*dx + dy*dy + dz*dz);
            const double dm1 = 1.0/dist;

            const double fac = 5.0*(2.0 - dist);

            fx += fac*dx*dm1;
            fy += fac*dy*dm1;
            fz += fac*dz*dm1;

          }

        }

      }

      __syncthreads();

    }

    if (i < N){

      F[3*i] = fx;
      F[3*i + 1] = fy;
      F[3*i + 2] = fz;

    }

  }

} // End of seeding kernel.
