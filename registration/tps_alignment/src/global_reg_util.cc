/*
Copyright 2007 Benedict Brown
Princeton University

global_reg_util.cc
Utility functions for global ragne scan registration.

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

#include "global_reg.h"
#include "jama_svd.h"

// grab the next file name from a list stored in a file
void read_fname(char *fname, FILE *flist) {
  int cpos = 0, c = 0;

  while (!feof(flist) && isspace(c = getc(flist)));
  if (!feof(flist)) fname[cpos++] = c;
  while (!feof(flist) && !isspace(c = getc(flist))) fname[cpos++] = (char) c;
  fname[cpos] = '\0';
}

// extract the base file name of a mesh
const char *get_mesh_name(const char *plyfilename) {
  int path_comp = strlen(plyfilename);
  while (path_comp && plyfilename[path_comp] != '/') path_comp--;
  return plyfilename + path_comp + 1;
}

// compute the spring error on a single point
// this is just the weighted sum of the errors on all the incident springs
//
// skip_unstable tells whether or not springs marked as unstable should be included
float point_error(int pnum, const spring_vector &springs, bool skip_unstable,
                  const vector<point> &targets, const vector<bool> &use_points) {
  float err = 0;
  float total_weight = 0;
  for (unsigned int i = 0; i < springs.size(); i++) {
    if (!use_points[springs[i].tgt]) continue;
    if (skip_unstable && !springs[i].stable) continue;
    float d = dist(targets[pnum], targets[springs[i].tgt]);
    err += springs[i].wgt * sqr(springs[i].len - d);
    total_weight += springs[i].wgt;
  }
  err = (total_weight > 0) ? err / total_weight : -1;

  return err;
}

// set the errors on all points
void set_errors(vector<float> &errors, const vector<spring_vector> &springs,
                const vector<point> &targets, const vector<bool> &use_points) {
  unsigned int num_points = errors.size();
  assert(springs.size() == num_points);

  for (unsigned int i = 0; i < num_points; i++)
    errors[i] = use_points[i] ? point_error(i, springs[i], true, targets, use_points) : -2;
}

// compute the best rigid-body mapping of x to y
void compute_rigid_alignment(const farr &x, const farr &y, farr &R, point &T) {
  assert(x.dim1() == y.dim1());
  assert(x.dim2() == 3);
  assert(y.dim2() == 3);

  // do rigid body warp only
  farr covar(3, 3, 0.0);
  point mx(0, 0, 0);
  point my(0, 0, 0);

  // compute mean
  for (int j = 0; j < x.dim1(); j++) {
    mx[0] += x[j][0];
    mx[1] += x[j][1];
    mx[2] += x[j][2];
    my[0] += y[j][0];
    my[1] += y[j][1];
    my[2] += y[j][2];
  }
  mx /= x.dim1();
  my /= x.dim1();

  // compute covariance
  for (int j = 0; j < x.dim1(); j++) {
    covar[0][0] += (y[j][0] - my[0]) * (x[j][0] - mx[0]);
    covar[0][1] += (y[j][0] - my[0]) * (x[j][1] - mx[1]);
    covar[0][2] += (y[j][0] - my[0]) * (x[j][2] - mx[2]);
    covar[1][0] += (y[j][1] - my[1]) * (x[j][0] - mx[0]);
    covar[1][1] += (y[j][1] - my[1]) * (x[j][1] - mx[1]);
    covar[1][2] += (y[j][1] - my[1]) * (x[j][2] - mx[2]);
    covar[2][0] += (y[j][2] - my[2]) * (x[j][0] - mx[0]);
    covar[2][1] += (y[j][2] - my[2]) * (x[j][1] - mx[1]);
    covar[2][2] += (y[j][2] - my[2]) * (x[j][2] - mx[2]);
  }

  float det = covar[0][0] * covar[1][1] * covar[2][2] +
    covar[0][1] * covar[1][2] * covar[2][0] +
    covar[0][2] * covar[1][0] * covar[2][1] -
    covar[0][2] * covar[1][1] * covar[2][0] -
    covar[0][1] * covar[1][0] * covar[2][2] -
    covar[0][0] * covar[1][2] * covar[2][1];

  // now do SVD of covar
  JAMA::SVD<matrix_t> svd(covar);
  farr U, V;
  svd.getU(U);
  svd.getV(V);
  R = farr(3, 3);

  fprintf(stderr, "Covariance det: %g\n", det);
  if (det < 0) {
    // assert(0);
    U[0][2] = -U[0][2];
    U[1][2] = -U[1][2];
    U[2][2] = -U[2][2];
  }

  // R = UV'
  for (int row = 0; row < 3; row++)
    for (int col = 0; col < 3; col++)
      R[row][col] = U[row][0] * V[col][0] + U[row][1] * V[col][1] + U[row][2] * V[col][2];

  point Rmx(R[0][0] * mx[0] + R[0][1] * mx[1] + R[0][2] * mx[2],
            R[1][0] * mx[0] + R[1][1] * mx[1] + R[1][2] * mx[2],
            R[2][0] * mx[0] + R[2][1] * mx[1] + R[2][2] * mx[2]);
  T = point(my[0] - Rmx[0], my[1] - Rmx[1], my[2] - Rmx[2]);
}

// write the vertices for a cube in a vertex array v
// the cube is centered at o and has width 2w
void write_cube_vertex(point *v, point &o, float w) {
  for (int i = 0; i < 8; i++) {
    v[i] = o;
    v[i][0] += (i & 4) ? -w : w;
    v[i][1] += (i & 2) ? -w : w;
    v[i][2] += (i & 1) ? -w : w;
  }
}

// Write a cube into the face array faces, starting at face index f
// and vertex index v.  It is assumed that the vertices are already
// created at the proper positions.
void write_cube_face(vector<TriMesh::Face> &faces, int f, int v) {
  faces[f + 0][0] = v + 0;
  faces[f + 0][1] = v + 3;
  faces[f + 0][2] = v + 1;
  faces[f + 1][0] = v + 0;
  faces[f + 1][1] = v + 2;
  faces[f + 1][2] = v + 3;

  faces[f + 2][0] = v + 1;
  faces[f + 2][1] = v + 3;
  faces[f + 2][2] = v + 7;
  faces[f + 3][0] = v + 1;
  faces[f + 3][1] = v + 7;
  faces[f + 3][2] = v + 5;

  faces[f + 4][0] = v + 4;
  faces[f + 4][1] = v + 0;
  faces[f + 4][2] = v + 1;
  faces[f + 5][0] = v + 4;
  faces[f + 5][1] = v + 1;
  faces[f + 5][2] = v + 5;

  faces[f + 6][0] = v + 2;
  faces[f + 6][1] = v + 0;
  faces[f + 6][2] = v + 4;
  faces[f + 7][0] = v + 2;
  faces[f + 7][1] = v + 4;
  faces[f + 7][2] = v + 6;

  faces[f + 8][0] = v + 7;
  faces[f + 8][1] = v + 6;
  faces[f + 8][2] = v + 4;
  faces[f + 9][0] = v + 7;
  faces[f + 9][1] = v + 4;
  faces[f + 9][2] = v + 5;

  faces[f + 10][0] = v + 3;
  faces[f + 10][1] = v + 6;
  faces[f + 10][2] = v + 7;
  faces[f + 11][0] = v + 3;
  faces[f + 11][1] = v + 2;
  faces[f + 11][2] = v + 6;
}

#ifdef COMPUTE_RMS_ERR
// Determine whether point i is an edge point
static inline bool edge(const TriMesh *s, int i) {
  return s->adjacentfaces[i].size() != s->neighbors[i].size();
}

// compute the RMS alignment error between a pair of meshes
// this is a useful number, but the computation is very slow
double rms_err(TriMesh *src, TriMesh *tgt, KDtree *tgt_kd) {
  float max_kd_dist = 5;
  unsigned int tested_points = 0;
  double err = 0;

  vector<bool> used_vertex(src->vertices.size());
  for (int i = 0; i < 1000; i++) {
    int j = 0;
    const float *p = NULL;
    point tmp_vert;
    do {
      tested_points++;

      // get_random_point(meshx, cum_prob, tmp_vert);
      size_t vnum = (size_t) (drand48() * src->vertices.size());
      if (vnum == src->vertices.size()) { tested_points--; continue; }
      tmp_vert = src->vertices[vnum];

      if (used_vertex[vnum]) {
        // if (tested_points < src->vertices.size())
        //   fprintf(stderr, "Selected already seen vertex %d at try %u\n", vnum, tested_points);
        p = NULL;
        continue;
      } else { used_vertex[vnum] = true; }

#if 1
      // disallow boundary points in source mesh
      if (edge(src, vnum)) {
        p = NULL;
        continue;
      }
#endif

      if (!(p = tgt_kd->closest_to_pt(tmp_vert, sqr(max_kd_dist)))) { continue; }

#if 1
      // disallow boundary points in target mesh
      j = (p - (const float *) &tgt->vertices[0][0]) / 3;
      if (edge(tgt, j)) {
        p = NULL;
        continue;
      }
#endif
    } while (!p && tested_points < BREAK_FACTOR * tgt->vertices.size());

    if (tested_points == BREAK_FACTOR * tgt->vertices.size()) break;

    err += sqr((tmp_vert - point(p)) DOT src->normals[j]);
  }

  return sqrt(err / tested_points);
}
#endif

