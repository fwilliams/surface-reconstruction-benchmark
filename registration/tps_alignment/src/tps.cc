/*
Copyright 2007 Benedict Brown
Princeton University

tps.cc
Routines for doing tps alignment.

-----

This file is part of tps_alignment.

tps_alignment is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.

tps_alignment is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <vector>
#include <queue>

#include "tnt_array2d.h"
#include "jama_lu.h"

#include "tps.h"

using namespace TNT;


#define SQ(x) ((x) * (x))

// #define PRINT_DATA

#if 1
static void print_matrix(const farr &M, FILE *f = stdout, const char *fmt = "%g ") {
  for (int row = 0; row < M.dim1(); row++) {
    for (int col = 0; col < M.dim2(); col++) {
      fprintf(f, fmt, M[row][col]);
    }
    putc('\n', f);
  }
  putc('\n', f);
}
#endif

// return a green function of |a - b|, but pass in |a - b|^2
// that's because it will save a multiple in the 2D case
static inline matrix_t green(matrix_t ab2, int dim) {
  switch (dim) {
  case 3: return sqrt(ab2); break;
  case 2: return ab2 > 0.00001f ? ab2 * log(sqrt(ab2)) : 0; break;
  }

  fprintf(stderr, "Green's function is only implemented for dimensions 2 and 3.  You asked for dim %d.\n", dim);
  abort();
  return 0;
}

// return Kw, generating K line-by-line
static farr gen_Kw(const farr &x, const farr &z, const farr &w) {
  int nx = x.dim1(), nz = z.dim1();
  assert(nz == w.dim1());

  farr out(nx, w.dim2());
  matrix_t tmp[nz];

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < nz; j++) {
      tmp[j] = 0;
      for (int k = 0; k < x.dim2(); k++)
        tmp[j] += SQ(x[i][k] - z[j][k]);

      tmp[j] = green(tmp[j], x.dim2());
    }

    for (int k = 0; k < w.dim2(); k++) {
      out[i][k] = 0;
      for (int j = 0; j < nz; j++) {
        out[i][k] += tmp[j] * w[j][k];
      }
    }
  }

  return out;
}

// convert to homogeneous coordinates:
// (x, y, z) -> (1, x, y, z)
static farr augment(const farr &x) {
  farr aug_x(x.dim1(), x.dim2() + 1);
  for (int row = 0; row < x.dim1(); row++) {
    for (int i = 0; i < x.dim2(); i++) aug_x[row][i] = x[row][i];
    aug_x[row][x.dim2()] = 1.0f;
  }

  return aug_x;
}

static farr my_matmult(const farr &a, const farr &b) {
  // assert(a.dim2() == b.dim1());
  if (a.dim2() != b.dim1()) {
    printf("a: (%d x %d), b: (%d x %d)\n",
           a.dim1(), a.dim2(), b.dim1(), b.dim2());
    assert(0);
  }

  int rows = a.dim1(), cols = b.dim2();

  farr out(rows, cols);

  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      register matrix_t sum = 0.0f;
      for (int i = 0; i < a.dim2(); i++)
        sum += a[row][i] * b[i][col];
      out[row][col] = sum;
    }
  }

  return out;
}

// compute the TPS warping parameters required to map point set X to Y
void update_transform_mat(const farr &x, const farr &y, matrix_t lambda,
                          farr &w, farr &A) {
  int nx  = x.dim1();
  int dim = x.dim2() + 1;
  farr K(nx + dim, nx + dim);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < nx; j++) {
      matrix_t sqsum = 0;
      for (int k = 0; k < dim - 1; k++)
        sqsum += SQ(x[i][k] - x[j][k]);
      K[i][j] = green(sqsum, dim - 1);
    }

    K[i][i] += lambda;

    for (int j = 0; j < dim - 1; j++) K[i][nx + j] = x[i][j];
    K[i][nx + dim - 1] = 1;
  }

  for (int i = nx; i < nx + dim - 1; i++) {
    for (int j = 0; j < nx; j++) K[i][j] = x[j][i - nx];
    for (int j = nx; j < nx + dim; K[i][j++] = 0);
  }

  for (int j = 0; j < nx; K[nx + dim - 1][j++] = 1);
  for (int j = nx; j < nx + dim; K[nx + dim - 1][j++] = 0);

  farr aug_y = farr(nx + dim, dim);
  for (int i = nx; i < nx + dim; i++)
    for (int j = 0; j < dim; aug_y[i][j++] = 0);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < dim - 1; j++) aug_y[i][j] = y[i][j];
    aug_y[i][dim - 1] = 1;
  }

  JAMA::LU<matrix_t> lu(K);
  farr out = lu.solve(aug_y);
  w = farr(nx,  dim);
  A = farr(dim, dim);

  for (int i = 0; i < nx; i++)
    for (int j = 0; j < dim; j++)
      w[i][j] = out[i][j];
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      A[i][j] = out[i + nx][j];

  printf("Ky =\n");
  print_matrix(aug_y, stdout, "%.5f ");
  printf("w =\n");
  print_matrix(w, stdout, "%.5f ");
  printf("A =\n");
  print_matrix(A, stdout, "%.5f ");
}

farr warp_points(const farr &pts, const farr &x,
                 const farr &w, const farr &A) {
  int n = pts.dim1();

  farr aug_pts = augment(pts);
  aug_pts = my_matmult(aug_pts, A);

  farr Kx = gen_Kw(pts, x, w);
  aug_pts += Kx;

  farr out(n, 3);
  for (int i = 0; i < n; i++) {
    out[i][0] = aug_pts[i][0];
    out[i][1] = aug_pts[i][1];
    out[i][2] = aug_pts[i][2];
  }

  return out;
}
