// SPDX-License-Identifier: CC0-1.0
#include <math.h>
#include <omp.h>
#include <stdio.h>

extern "C" void hotspot3d_target_c(float *tin, float *pin, float *tout,
                                   int numCols, int numRows, int layers, int iterations,
                                   float ce, float cw, float cn, float cs, float ct, float cb,
                                   float cc, float stepDivCap, double *kernel_time)
{
  int size = numCols * numRows * layers;
  float *tIn = tin;
  float *tOut = tout;

#pragma omp target data map(to: tIn[0:size], pin[0:size]) map(alloc: tOut[0:size])
  {
    double kstart = omp_get_wtime();

    for (int iter = 0; iter < iterations; iter++) {
#pragma omp target teams distribute parallel for collapse(2) thread_limit(256)
      for (int j = 0; j < numRows; j++) {
        for (int i = 0; i < numCols; i++) {
          float amb_temp = 80.0;

          int c = i + j * numCols;
          int xy = numCols * numRows;

          int W = (i == 0) ? c : c - 1;
          int E = (i == numCols - 1) ? c : c + 1;
          int N = (j == 0) ? c : c - numCols;
          int S = (j == numRows - 1) ? c : c + numCols;

          float temp1, temp2, temp3;
          temp1 = temp2 = tIn[c];
          temp3 = tIn[c + xy];
          tOut[c] = cc * temp2 + cw * tIn[W] + ce * tIn[E] + cs * tIn[S] +
                    cn * tIn[N] + cb * temp1 + ct * temp3 + stepDivCap * pin[c] + ct * amb_temp;
          c += xy;
          W += xy;
          E += xy;
          N += xy;
          S += xy;

          for (int k = 1; k < layers - 1; ++k) {
            temp1 = temp2;
            temp2 = temp3;
            temp3 = tIn[c + xy];
            tOut[c] = cc * temp2 + cw * tIn[W] + ce * tIn[E] + cs * tIn[S] +
                      cn * tIn[N] + cb * temp1 + ct * temp3 + stepDivCap * pin[c] + ct * amb_temp;
            c += xy;
            W += xy;
            E += xy;
            N += xy;
            S += xy;
          }
          temp1 = temp2;
          temp2 = temp3;
          tOut[c] = cc * temp2 + cw * tIn[W] + ce * tIn[E] + cs * tIn[S] +
                    cn * tIn[N] + cb * temp1 + ct * temp3 + stepDivCap * pin[c] + ct * amb_temp;
        }
      }
      auto temp = tIn;
      tIn = tOut;
      tOut = temp;
    }

    double kend = omp_get_wtime();
    *kernel_time = kend - kstart;

    if (iterations & 01) {
#pragma omp target update from(tIn[0:size])
    } else {
#pragma omp target update from(tOut[0:size])
    }
  }
}

extern "C" void compute_temp_cpu_c(float *pIn, float *tIn, float *tOut,
                                   int nx, int ny, int nz, float Cap,
                                   float Rx, float Ry, float Rz,
                                   float dt, float amb_temp, int numiter)
{
  float ce, cw, cn, cs, ct, cb, cc;
  float stepDivCap = dt / Cap;
  ce = cw = stepDivCap / Rx;
  cn = cs = stepDivCap / Ry;
  ct = cb = stepDivCap / Rz;

  cc = 1.0 - (2.0 * ce + 2.0 * cn + 3.0 * ct);

  int c, w, e, n, s, b, t;
  int x, y, z;
  int i = 0;
  do {
    for (z = 0; z < nz; z++)
      for (y = 0; y < ny; y++)
        for (x = 0; x < nx; x++) {
          c = x + y * nx + z * nx * ny;

          w = (x == 0) ? c : c - 1;
          e = (x == nx - 1) ? c : c + 1;
          n = (y == 0) ? c : c - nx;
          s = (y == ny - 1) ? c : c + nx;
          b = (z == 0) ? c : c - nx * ny;
          t = (z == nz - 1) ? c : c + nx * ny;

          tOut[c] = tIn[c] * cc + tIn[n] * cn + tIn[s] * cs + tIn[e] * ce + tIn[w] * cw +
                    tIn[t] * ct + tIn[b] * cb + (dt / Cap) * pIn[c] + ct * amb_temp;
        }
    float *temp = tIn;
    tIn = tOut;
    tOut = temp;
    i++;
  } while (i < numiter);
}

extern "C" float accuracy_c(float *arr1, float *arr2, int len)
{
  float err = 0.0;
  int i;
  for (i = 0; i < len; i++) {
    err += (arr1[i] - arr2[i]) * (arr1[i] - arr2[i]);
  }

  return (float)sqrt(err / len);
}

extern "C" void print_rms_c(float rms)
{
  printf("Root-mean-square error: %e\n", rms);
}
