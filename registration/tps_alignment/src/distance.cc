/*
Copyright 2007 Benedict Brown
Princeton University

Compute histgram of distances of range scan points to final mesh

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

#include <assert.h>
#include "XForm.h"
#include "TriMesh.h"
#include "KDtree.h"

#if 0
#define MAX_KD_DIST (10 * 10)
#define HIST_SIZE 1500
#define HIST_CELL .002
#else
// for sea floor
#define MAX_KD_DIST (5 * 5)
#define HIST_SIZE 250
#define HIST_CELL .004
#endif

// Find closest point to p on segment from v0 to v1
point closest_on_segment(const point &v0, const point &v1, const point &p) {
  vec v01 = v1 - v0;
  float d = (p - v0) DOT v01;
  d /= len2(v01);
  if (d < 0.0f)
    d = 0.0f;
  else if (d > 1.0f)
    d = 1.0f;
  return v0 + d * v01;
}


// Find closest point to p on face i of mesh
point closest_on_face(const TriMesh *mesh, int i, const point &p) {
  const TriMesh::Face &f = mesh->faces[i];
  const point &v0 = mesh->vertices[f[0]];
  const point &v1 = mesh->vertices[f[1]];
  const point &v2 = mesh->vertices[f[2]];
  vec a = v1 - v0, b = v2 - v0, p1 = p - v0, n = a CROSS b;

  float A[3][3] = { { a[0], b[0], n[0] },
                    { a[1], b[1], n[1] },
                    { a[2], b[2], n[2] } };
  float x[3] = { p1[0], p1[1], p1[2] };
  int indx[3];
  ludcmp(A, indx);
  lubksb(A, indx, x);

  if (x[0] >= 0.0f && x[1] >= 0.0f && x[0] + x[1] <= 1.0f)
    return v0 + x[0] * a + x[1] * b;

  point c01 = closest_on_segment(v0, v1, p);
  point c12 = closest_on_segment(v1, v2, p);
  point c20 = closest_on_segment(v2, v0, p);
  float d01 = dist2(c01, p);
  float d12 = dist2(c12, p);
  float d20 = dist2(c20, p);
  if (d01 < d12) {
    if (d01 < d20) return c01; else return c20;
  } else {
    if (d12 < d20) return c12; else return c20;
  }
}


// Find (good approximation to) closest point on mesh to p.
// Finds closest vertex, then checks all faces that touch it.
// If the mesh has no faces, returns the closest vertex
// Returns nearest vertex num, with nearest point in out
// Return < 0 for no match
int closest_pt(const TriMesh *mesh, const KDtree *kd, float mdist, const point &p, point &out) {
  // float maxdist2 = sqr(dist(p, mesh->bsphere.center) + mesh->bsphere.r);
  const float *match = kd->closest_to_pt(p, mdist);
  if (!match) return -1;

  int ind = (match - (const float *) &(mesh->vertices[0][0])) / 3;
  if (ind < 0 || ind >= (int) mesh->vertices.size()) {
    fprintf(stderr, "This can't happen - bogus index\n");
    return -1;
  }

  if (mesh->faces.size() == 0)  {
    out = mesh->vertices[ind];
    return ind;
  } else assert(mesh->adjacentfaces.size());

  const vector<int> &a = mesh->adjacentfaces[ind];
  if (a.empty()) {
    fprintf(stderr, "Found match to unconnected point\n");
    out =  mesh->vertices[ind];
    return ind;
  }

  float closest_dist2 = 3.3e33;
  for (unsigned int i = 0; i < a.size(); i++) {
    point c1 = closest_on_face(mesh, a[i], p);
    float this_dist2 = dist2(c1, p);
    if (this_dist2 < closest_dist2) {
      closest_dist2 = this_dist2;
      out = c1;
    }
  }
  return ind;
}

int main(int argc, char **argv) {
  long num_points = 0;
  double total_dist = 0.0;
  double total_sq_dist = 0.0;

  double mean, stddev;

  vector<int> hist(HIST_SIZE);

  // read in master mesh
  TriMesh *mesh = TriMesh::read(argv[1]);
  mesh->need_neighbors();
  mesh->need_adjacentfaces();

  KDtree *kd = new KDtree(mesh->vertices);

  while (!feof(stdin)) {
    char fname[1024];
    scanf(" %s ", fname);
    TriMesh *m = TriMesh::read(fname);

    for (unsigned int i = 0; i < m->vertices.size(); i++) {
#if 0
      const float *p = kd->closest_to_pt(m->vertices[i], MAX_KD_DIST);
      if (!p) continue; // skip outliers

      int j = (p - (const float *) &mesh->vertices[0]) / 3;
#else
      point p;
      int j = closest_pt(mesh, kd, MAX_KD_DIST, m->vertices[i], p);    
      if (j < 0) continue;
#endif

      // skip boundary points
      if (mesh->neighbors[j].size() != mesh->adjacentfaces[j].size()) continue;
      num_points++;
      double d2 = dist2(m->vertices[i], mesh->vertices[j]);
      double d = sqrt(d2);
      total_sq_dist += d2;
      total_dist += d;

      int bin = (int) (d / HIST_CELL);
      if (bin < HIST_SIZE) hist[bin]++;
    }

    mean = total_dist / num_points;
    stddev = sqrt(total_sq_dist / num_points - sqr(mean));
    printf("Analyzed %ld points.\nMean: %f StdDev: %f\n", num_points, mean, stddev);

    delete m;
  }

  // write histogram
  FILE *h = fopen(argv[2], "w");
  for (int i = 0; i < HIST_SIZE; i++)
    fprintf(h, "%d\n", hist[i]);

#if 0
  mean = total_dist / num_points;
  var = sqrt(total_sq_dist / num_points - sqr(mean));
  printf("Analyzed %ld points.\nMean: %f StdDev: %f\n", num_points, mean, var);
#endif

  return 0;
}
