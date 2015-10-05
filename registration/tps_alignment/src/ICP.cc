/*
Copyright 2007 Szymon Rusinkiewicz and Benedict Brown
Princeton University

ICP.cc
Routines for doing ICP.

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

#include <cmath>
#include <algorithm>
#include <stdio.h>
#include "ICP.h"
#include "KDtree.h"
#include "timestamp.h"
#include "lineqn.h"
using namespace std;

#define MAX_ITERS 100
#define MIN_PAIRS 25
#define DESIRED_PAIRS 500
#define COMPAT_THRESH 0.9f
#define CURVE_THRESH  0.99f // higher is more restrictive
// #define CURVE_THRESH  0.0f // higher is more restrictive
#define TERM_THRESH 3
#define TERM_HIST 7
#define EIG_THRESH 0.01f

// One or both of the following can be #defined
#define USE_GRID_FOR_OVERLAPS
#define USE_KD_FOR_OVERLAPS


// Determine whether point i is an edge point
static inline bool edge(const EdgeMesh *s, int i) {
  // return s->adjacentfaces[i].size() != s->neighbors[i].size();
  return s->isedge(i);
}


// A pair of points, with an associated normal
struct PtPair {
  point p1, p2;
  vec norm1, norm2;
  PtPair(const point &p1_, const point &p2_,
         const vec &norm1_, const vec &norm2_) :
    p1(p1_), p2(p2_), norm1(norm1_), norm2(norm2_) {}
};


// A class for evaluating compatibility of normals during KDtree searches
class NormCompat : public KDtree::CompatFunc {
private:
  const vec n;
  float dot;
  const EdgeMesh *s;

public:
  NormCompat(const vec &n_, float dot_, const EdgeMesh *s_):
    n(n_), dot(dot_), s(s_) {}

  virtual bool operator () (const float *p) const {
    int idx = (const point *)p - (const point *)&(s->vertices[0]);
    if (edge(s, idx))
      return true;
    else
      return (n DOT s->normals[idx]) > dot;
  }
};

// compute compatibility of curvature and normals both
// BUG: If curvatures are very close to 0, this could fail
class CurveAndNormCompat : public KDtree::CompatFunc {
private:
  NormCompat nc;
  const EdgeMesh *s;
  float curvature, icurvature;
  float curve_ratio, icurve_ratio;

public:
  CurveAndNormCompat(float curv1, float curv2, float curve_ratio_,
                     const vec &n_, float dot_, const EdgeMesh *s_) :
    nc(n_, dot_, s_), s(s_), curvature(0.5f * (curv1 + curv2)) {
    icurvature = 1 / curvature;
    if (curve_ratio_ > 1) {
      curve_ratio = curve_ratio_;
      icurve_ratio = 1 / curve_ratio_;
    } else {
      curve_ratio = 1 / curve_ratio_;
      icurve_ratio = curve_ratio_;
    }
  }

  virtual bool operator() (const float *p) const {
    int idx = (const point *) p - (const point *) &(s->vertices[0]);
    if (edge(s, idx)) return true;
    else {
      bool rc = nc(p);
      if (!rc) return false;
      float ratio = 0.5f * (s->curv1[idx] + s->curv2[idx]) * icurvature;
      rc = (ratio < curve_ratio) && (ratio > icurve_ratio);
      return isnan(ratio) || rc;
    }

    // return false;
  }
};

// Find the median squared distance between points
static float median_dist2(const vector<PtPair> &pairs) {
  size_t n = pairs.size();
  if (!n)
    return 0.0f;

  vector<float> distances2;
  distances2.reserve(n);
  for (size_t i = 0; i < n; i++)
    distances2.push_back(dist2(pairs[i].p1, pairs[i].p2));

  size_t pos = n / 2;
  nth_element(distances2.begin(),
              distances2.begin() + pos,
              distances2.end());
  return distances2[pos];
}


// A spatial grid datastructure for fast overlap computation
class Grid {
public:
  enum { GRID_SHIFT = 4, GRID_MAX = (1 << GRID_SHIFT) - 1 };
  float xmin, xmax, ymin, ymax, zmin, zmax, scale;
  char g[1 << 3*GRID_SHIFT];
  bool valid(const point &p) {
    return p[0] > xmin && p[1] > ymin && p[2] > zmin &&
           p[0] < xmax && p[1] < ymax && p[2] < zmax;
  }

  int ind(const point &p) {
    int x = min(int(scale * (p[0] - xmin)), (int) GRID_MAX);
    int y = min(int(scale * (p[1] - ymin)), (int) GRID_MAX);
    int z = min(int(scale * (p[2] - zmin)), (int) GRID_MAX);
    return (x << (2*GRID_SHIFT)) + (y << GRID_SHIFT) + z;
  }

  bool overlaps(const point &p) { return valid(p) && g[ind(p)]; }
  Grid(const vector<point> &pts);
};


// Compute a Grid from a list of points
Grid::Grid(const vector<point> &pts)
{
  memset(g, 0, sizeof(g));
  xmin = xmax = pts[0][0];
  ymin = ymax = pts[0][1];
  zmin = zmax = pts[0][2];
  for (size_t i = 1; i < pts.size(); i++) {
    if (pts[i][0] < xmin)  xmin = pts[i][0];
    if (pts[i][0] > xmax)  xmax = pts[i][0];
    if (pts[i][1] < ymin)  ymin = pts[i][1];
    if (pts[i][1] > ymax)  ymax = pts[i][1];
    if (pts[i][2] < zmin)  zmin = pts[i][2];
    if (pts[i][2] > zmax)  zmax = pts[i][2];
  }
  scale = 1.0f / max(max(xmax-xmin, ymax-ymin), zmax-zmin);
  scale *= float(1 << GRID_SHIFT);
  for (size_t i = 0; i < pts.size(); i++) {
    int idx = ind(pts[i]);
    g[idx] = 1;
  }
}


// Determine which points on s1 and s2 overlap the other, filling in o1 and o2
// Also fills in maxdist, if it is <= 0 on input
void compute_overlaps(ICP_struct &m0, ICP_struct &m1, float &maxdist, int verbose) {
  ICP_struct *m[2] = { &m0, &m1 };
  size_t nv[2];
  xform  xf[2];
  Grid *g[2];

  timestamp t = now();

  for (int idx = 0; idx < 2; idx++) {
    m[idx]->m->need_normals();
    m[idx]->m->need_bbox();
    nv[idx] = m[idx]->m->vertices.size();
    g[idx] = new Grid(m[idx]->m->vertices);

    // compute transform from other scan to this one
    xf[idx] = inv(m[1 - idx]->xf) * m[idx]->xf;
  }

  if (maxdist <= 0.0f) {
    maxdist = min(1.0f / g[0]->scale, 1.0f / g[1]->scale);
    fprintf(stderr, "new maxidst: %f\n", maxdist);
  }
  float maxdist2 = sqr(maxdist);

  for (int idx = 0; idx < 2; idx++) {
    int odx = 1 - idx;

    m[idx]->weights.resize(nv[idx]);

    for (size_t i = 0; i < nv[idx]; i++) {
      m[idx]->weights[i] = 0;
      if (edge(m[idx]->m, i)) continue;
      point p = xf[idx] * m[idx]->m->vertices[i];
#ifdef USE_GRID_FOR_OVERLAPS
      if (!g[odx]->overlaps(p))
        continue;
#endif
#ifdef USE_KD_FOR_OVERLAPS
      const float *match = m[odx]->kd->closest_to_pt(p, maxdist2);
      if (!match)
        continue;
      if (edge(m[odx]->m, ((match - (const float *) &m[odx]->m->vertices[0][0])) / 3))
        continue;
#endif
      m[idx]->weights[i] = 1;
    }
  }

  delete g[0];
  delete g[1];

  if (verbose > 1) {
    fprintf(stderr, "Computed overlaps in %.2f msec.\n",
            (now() - t) * 1000.0f);
  }
}


#if 1
// XXX - HACK - record point pairs
static void save_pairs(const vector<PtPair> &pairs, vector<point> &p) {
  p.resize(pairs.size());
  for (size_t i = 0; i < pairs.size(); i++) {
    p[i][0] = pairs[i].p1[0];
    p[i][1] = pairs[i].p1[1];
    p[i][2] = pairs[i].p1[2];
  }
}
#else
static void save_pairs(const vector<PtPair> &pairs) {}
#endif


// Select a number of points and find correspondences 
static void select_and_match_cdf(ICP_struct &m0, ICP_struct &m1, float incr, float maxdist, int verbose,
                                 vector<PtPair> &pairs, vector< pair<int, int> > &ipairs, bool flip) {
  xform xf0r  = norm_xf(m0.xf);
  xform xf1r  = norm_xf(m1.xf);
  xform xf01  = inv(m1.xf) * m0.xf;
  xform xf01r = norm_xf(xf01);

  assert(m0.cdf.size() == m0.m->vertices.size());
  assert(m1.cdf.size() == m1.m->vertices.size());
  float maxdist2 = sqr(maxdist);

  size_t i = 0;
  float cval = 0.0f;

  // ipairs.clear();
  while (1) {
    cval += incr * drand48();
    if (cval >= 1.0f) break;
    while (m0.cdf[i] <= cval) i++;
    cval = m0.cdf[i];
    // assert(!edge(s1, i));
    if (edge(m0.m, i)) continue;

    point p = xf01  * m0.m->vertices[i];
    vec   n = xf01r * m0.m->normals[i];

    // Do the matching
    // NormCompat nc(n, COMPAT_THRESH, m1.m);
    CurveAndNormCompat cc(m0.m->curv1[i], m0.m->curv2[i], CURVE_THRESH, n, COMPAT_THRESH, m1.m);

    // const float *match = kd2->closest_to_pt(p, maxdist2, &nc);
    const float *match = m1.kd->closest_to_pt(p, maxdist2, &cc);
    if (!match) continue;
    int imatch = (match - (const float *) &(m1.m->vertices[0][0])) / 3;
    if (edge(m1.m, imatch)) continue;

    assert(i < m0.m->vertices.size());
    assert(imatch < (int) m1.m->vertices.size());

    // Project both points into world coords and save 
    if (flip) {
      ipairs.push_back(pair<int, int>(imatch, (int) i));
      pairs.push_back(PtPair(m1.xf * m1.m->vertices[imatch],
                             m0.xf * m0.m->vertices[i],
                             xf1r  * m1.m->normals[imatch],
                             xf0r  * m0.m->normals[i]));
    } else {
      ipairs.push_back(pair<int, int>((int) i, imatch));
      pairs.push_back(PtPair(m0.xf * m0.m->vertices[i],
                             m1.xf * m1.m->vertices[imatch],
                             xf0r  * m0.m->normals[i],
                             xf1r  * m1.m->normals[imatch]));
    }
  }

  // fprintf(stderr, "Computed %u pairs (non-probability).\n", ipairs.size());
}


// Select a number of points and find correspondences 
static void select_and_match_idx(ICP_struct &m0, ICP_struct &m1, size_t desired_pairs,
                                 float maxdist, int verbose, vector<PtPair> &pairs,
                                 vector< pair<int, int> > &ipairs, bool flip) {
  xform xf0r  = norm_xf(m0.xf);
  xform xf1r  = norm_xf(m1.xf);
  xform xf01  = inv(m1.xf) * m0.xf;
  xform xf01r = norm_xf(xf01);
  float maxdist2 = sqr(maxdist);

  // assert(indexes.size() == 6);
  // assert(indexes[0].size() == flip ? s2->vertices.size() : s1->vertices.size());

  // ipairs.clear();

  float eval_total = m0.evals[0] + m0.evals[1] + m0.evals[2] +
                     m0.evals[3] + m0.evals[4] + m0.evals[5];

  size_t num_pairs[6];
  for (int i = 0; i < 6; i++)
    num_pairs[i] = (int) ((desired_pairs - pairs.size()) * m0.evals[i] / eval_total + 0.5f);
#if 0
  fprintf(stderr, "%d + %d + %d + %d + %d + %d = %d\n",
          num_pairs[0], num_pairs[1], num_pairs[2],
          num_pairs[3], num_pairs[4], num_pairs[5],
          num_pairs[0]+ num_pairs[1]+ num_pairs[2] +
          num_pairs[3]+ num_pairs[4]+ num_pairs[5]);
#endif

  vector<bool> taken(m0.m->vertices.size());
  for (int ilist = 0; ilist < 6; ilist++) {
    size_t end_pairs = pairs.size() + num_pairs[ilist];
    for (int idx = 0; pairs.size() < end_pairs; idx++) {
      if ((int) m0.indexes[ilist].size() <= idx) break;

      int i = m0.indexes[ilist][idx];

      if (taken[i]) {
        // fprintf(stderr, "taken: %d %d %d\n", ilist, idx, indexes[ilist][idx]);
        continue;
      }

      // assert(!edge(s1, i));
      if (edge(m0.m, i)) {
        // fprintf(stderr, "edge: %d %d %d\n", ilist, idx, indexes[ilist][idx]);
        continue;
      }

      point p = xf01  * m0.m->vertices[i];
      vec   n = xf01r * m0.m->normals[i];

      // Do the matching
      // NormCompat nc(n, COMPAT_THRESH, m1.m);
      CurveAndNormCompat cc(m0.m->curv1[i], m0.m->curv2[i], CURVE_THRESH, n, COMPAT_THRESH, m1.m);

      // const float *match = kd2->closest_to_pt(p, maxdist2, &nc);
      const float *match = m1.kd->closest_to_pt(p, maxdist2, &cc);
      if (!match) {
        continue;
      }
      int imatch = (match - (const float *) &(m1.m->vertices[0][0])) / 3;
      if (edge(m1.m, imatch)) {
        continue;
      }

      assert(i < (int) m0.m->vertices.size());
      assert(imatch < (int) m1.m->vertices.size());

      taken[i] = true;

      // Project both points into world coords and save 
      if (flip) {
        ipairs.push_back(pair<int, int>(imatch, (int) i));
        pairs.push_back(PtPair(m1.xf * m1.m->vertices[imatch],
                               m0.xf * m0.m->vertices[i],
                               xf1r  * m1.m->normals[imatch],
                               xf0r  * m0.m->normals[i]));
      } else {
        ipairs.push_back(pair<int, int>((int) i, imatch));
        pairs.push_back(PtPair(m0.xf * m0.m->vertices[i],
                               m1.xf * m1.m->vertices[imatch],
                               xf0r  * m0.m->normals[i],
                               xf1r  * m1.m->normals[imatch]));
      }
    }
  }

  // fprintf(stderr, "Computed %u pairs (probability).\n", ipairs.size());
}


// Compute stability probabilities for each mesh (for sampling)
static void compute_stable_probs(const vector<PtPair> &pairs,
                                 float evec1[6][6], float eval1[6], point &centroid1, float &scale1,
                                 float evec2[6][6], float eval2[6], point &centroid2, float &scale2) {
  size_t n = pairs.size();

  centroid1 = centroid2 = point(0,0,0);
  for (size_t i = 0; i < n; i++) {
    centroid1 += pairs[i].p1;
    centroid2 += pairs[2].p1;
  }

  centroid1 /= float(n);
  centroid2 /= float(n);

  scale1 = scale2 = 0.0f;
  for (size_t i = 0; i < n; i++) {
    scale1 += dist2(pairs[i].p1, centroid1);
    scale2 += dist2(pairs[i].p2, centroid2);
  }
  scale1 /= float(n);
  scale2 /= float(n);

  scale1 = 1.0f / sqrt(scale1);
  scale2 = 1.0f / sqrt(scale2);

  memset(&evec1[0][0], 0, 6*6*sizeof(float));
  memset(&evec2[0][0], 0, 6*6*sizeof(float));

  for (size_t i = 0; i < n; i++) {
    const point &p1 = pairs[i].p1;
    const point &p2 = pairs[i].p2;
    const vec &n1 = pairs[i].norm1;
    const vec &n2 = pairs[i].norm2;

    vec p1c = p1 - centroid1;
    p1c *= scale1;
    vec c1 = p1c CROSS n1;

    vec p2c = p2 - centroid2;
    p2c *= scale2;
    vec c2 = p2c CROSS n2;

    float x1[6] = { c1[0], c1[1], c1[2], n1[0], n1[1], n1[2] };
    float x2[6] = { c2[0], c2[1], c2[2], n2[0], n2[1], n2[2] };
    for (int j = 0; j < 6; j++) {
      for (int k = 0; k < 6; k++) {
        evec1[j][k] += x1[j] * x1[k];
        evec2[j][k] += x2[j] * x2[k];
      }
    }
  }

  eigdc(evec1, eval1);
  eigdc(evec2, eval2);
}


// Compute ICP alignment matrix, including eigenvector decomposition
static void compute_ICPmatrix(const vector<PtPair> &pairs,
                              float evec[6][6], float eval[6], float b[6],
                              point &centroid, float &scale, float &err) {
  size_t n = pairs.size();

  centroid = point(0,0,0);
  for (size_t i = 0; i < n; i++)
    centroid += pairs[i].p2;
  centroid /= float(n);

  scale = 0.0f;
  for (size_t i = 0; i < n; i++)
    scale += dist2(pairs[i].p2, centroid);
  scale /= float(n);
  scale = 1.0f / sqrt(scale);

  memset(&evec[0][0], 0, 6*6*sizeof(float));
  memset(&b[0], 0, 6*sizeof(float));

  err = 0.0f;
  for (size_t i = 0; i < n; i++) {
    const point &p1 = pairs[i].p1;
    const point &p2 = pairs[i].p2;
    const vec &n = pairs[i].norm1;

    float d = (p1 - p2) DOT n;
    d *= scale;
    vec p2c = p2 - centroid;
    p2c *= scale;
    vec c = p2c CROSS n;

    err += d * d;
    float x[6] = { c[0], c[1], c[2], n[0], n[1], n[2] };
    for (int j = 0; j < 6; j++) {
      b[j] += d * x[j];
      for (int k = 0; k < 6; k++)
        evec[j][k] += x[j] * x[k];
    }
  }

  err /= float(n);
  err = sqrt(err) / scale;
  eigdc(evec, eval);
}

// Compute ICP alignment, given matrix computed by compute_ICPmatrix
static void compute_alignxf(float evec[6][6], float eval[6], float b[6],
                            point &centroid, float scale, xform &alignxf)
{
        float einv[6];
        for (int i = 0; i < 6; i++) {
                if (eval[i] < EIG_THRESH * eval[5])
                        einv[i] = 0.0;
                else
                        einv[i] = 1.0 / eval[i];
        }
        float x[6];
        eigmult(evec, einv, b, x);

        // Interpret results
        float sx = min(max(x[0], -1.0f), 1.0f);
        float sy = min(max(x[1], -1.0f), 1.0f);
        float sz = min(max(x[2], -1.0f), 1.0f);
        float cx = sqrt(1.0f - sx*sx);
        float cy = sqrt(1.0f - sy*sy);
        float cz = sqrt(1.0f - sz*sz);

        alignxf[0]  = cy*cz;
        alignxf[1]  = sx*sy*cz + cx*sz;
        alignxf[2]  = -cx*sy*cz + sx*sz;
        alignxf[3]  = 0;
        alignxf[4]  = -cy*sz;
        alignxf[5]  = -sx*sy*sz + cx*cz;
        alignxf[6]  = cx*sy*sz + sx*cz;
        alignxf[7]  = 0;
        alignxf[8]  = sy;
        alignxf[9]  = -sx*cy;
        alignxf[10] = cx*cy;
        alignxf[11] = 0;
        alignxf[12] = x[3] / scale + centroid[0] - alignxf[0]*centroid[0] -
                      alignxf[4]*centroid[1] - alignxf[8]*centroid[2];
        alignxf[13] = x[4] / scale + centroid[1] - alignxf[1]*centroid[0] -
                      alignxf[5]*centroid[1] - alignxf[9]*centroid[2];
        alignxf[14] = x[5] / scale + centroid[2] - alignxf[2]*centroid[0] -
                      alignxf[6]*centroid[1] - alignxf[10]*centroid[2];
        alignxf[15] = 1;
}


// Compute anisotropic scale
void compute_scale(const vector<PtPair> &pairs, xform &alignxf, int verbose)
{
        int n = pairs.size();

        // Compute COM
        point centroid;
        for (int i = 0; i < n; i++)
                centroid += pairs[i].p1 + pairs[i].p2;
        centroid /= 2.0f * n;

        // Compute covariance matrices
        double cov1[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
        double cov2[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
        for (int i = 0; i < n; i++) {
                vec p = pairs[i].p1 - centroid;
                for (int j = 0; j < 3; j++)
                        for (int k = 0; k < 3; k++)
                                cov1[j][k] += p[j]*p[k];
                p = pairs[i].p2 - centroid;
                for (int j = 0; j < 3; j++)
                        for (int k = 0; k < 3; k++)
                                cov2[j][k] += p[j]*p[k];
        }

        // Compute sqrt of cov
        double csqrt1[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };
        double csqrt2[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };
        double eval[3];
        eigdc(cov1, eval);
        for (int i = 0; i < 3; i++) {
          if (unlikely(eval[i] <= 0)) {
            alignxf = xform();
            return;
            // assert(0);
          }
          eval[i] = sqrt(eval[i] / n);
        }
        for (int i = 0; i < 3; i++)
          eigmult(cov1, eval, csqrt1[i], csqrt1[i]);

        eigdc(cov2, eval);
        for (int i = 0; i < 3; i++) {
          if (unlikely(eval[i] <= 0)) {
            alignxf = xform();
            return;
            // assert(0);
          }
          eval[i] = sqrt(eval[i] / n);
        }
        for (int i = 0; i < 3; i++)
                eigmult(cov2, eval, csqrt2[i], csqrt2[i]);

        if (verbose > 1) {
                fprintf(stderr, "sqrt covariance 1 =");
                for (int j = 0; j < 3; j++) {
                        printf("\n\t");
                        for (int k = 0; k < 3; k++)
                                fprintf(stderr, "%10.3f ", csqrt1[j][k]);
                }
                fprintf(stderr, "\nsqrt covariance 2 =");
                for (int j = 0; j < 3; j++) {
                        printf("\n\t");
                        for (int k = 0; k < 3; k++)
                                fprintf(stderr, "%10.3f ", csqrt2[j][k]);
                }
                printf("\n");
        }


        xform cxf1 = xform(csqrt1[0][0], csqrt1[1][0], csqrt1[2][0], 0,
                           csqrt1[0][1], csqrt1[1][1], csqrt1[2][1], 0,
                           csqrt1[0][2], csqrt1[1][2], csqrt1[2][2], 0,
                           0, 0, 0, 1);
        xform cxf2 = xform(csqrt2[0][0], csqrt2[1][0], csqrt2[2][0], 0,
                           csqrt2[0][1], csqrt2[1][1], csqrt2[2][1], 0,
                           csqrt2[0][2], csqrt2[1][2], csqrt2[2][2], 0,
                           0, 0, 0, 1);
        xform txf = xform::trans(centroid);

        alignxf = txf * cxf1 * inv(cxf2) * inv(txf);
}


// Do one iteration of ICP
static float ICP_iter(ICP_struct &m0, ICP_struct &m1, float &maxdist, int verbose,
                      float &stability, float &max_stability, float &incr,
                      vector<point> *point_vector, vector< pair<int, int> > &ipairs,
                      bool use_indexes, bool do_affine) {
  // Compute pairs
  timestamp t1 = now();
  static vector<PtPair> pairs;        
  pairs.clear();
  ipairs.clear();

  ICP_struct *m[2] = { &m0, &m1 };

  if (use_indexes) {
    select_and_match_idx(m0, m1, DESIRED_PAIRS / 2, maxdist,
                         verbose, pairs, ipairs, false);
    select_and_match_idx(m1, m0, DESIRED_PAIRS, maxdist,
                         verbose, pairs, ipairs, true);
  } else {
    select_and_match_cdf(m0, m1, incr, maxdist, verbose, pairs, ipairs, false);
    select_and_match_cdf(m1, m0, incr, maxdist, verbose, pairs, ipairs, true);
  }

  timestamp t2 = now();
  size_t np = pairs.size();
  if (verbose > 1) {
    fprintf(stderr, "Generated %u pairs in %.2f msec.\n",
            (unsigned int) np, (t2-t1) * 1000.0f);
  }

  // Reject pairs with distance > 2.5 sigma
  float thresh = 13.73818f * median_dist2(pairs) + 0.0001f;
  // float thresh = 3.0f * median_dist2(pairs) + 0.0001f;
  size_t next = 0;
  for (size_t i = 0; i < np; i++) {
    if (dist2(pairs[i].p1, pairs[i].p2) < thresh)
      pairs[next++] = pairs[i];
  }
  pairs.erase(pairs.begin() + next, pairs.end());

  timestamp t3 = now();

  if (verbose > 1) {
    fprintf(stderr, "Rejected %u pairs in %.2f msec.\n",
            (unsigned int) (np - pairs.size()), (t3-t2) * 1000.0f);
  }

  if (pairs.size() < MIN_PAIRS) {
    if (verbose)
      fprintf(stderr, "Too few point pairs.\n");
    return -1.0f;
  }

  // Update incr and maxdist based on what happened here
  incr *= (float) pairs.size() / DESIRED_PAIRS;
  maxdist = sqrt(thresh);

  // Do the minimization
  float evec[6][6], eval[6], b[6], scale, err;
  point centroid;
  xform alignxf;
  compute_ICPmatrix(pairs, evec, eval, b, centroid, scale, err);
  if (verbose > 1) {
    fprintf(stderr, "RMS point-to-plane error = %f\n", err);
    for (int i = 0; i < 5; i++) {
      if (eval[i] < EIG_THRESH * eval[5])
        fprintf(stderr, "Small eigenvalue %f (largest is %f)\n", eval[i], eval[5]);
    }
  }

  // if (check_stability && eval[0] < EIG_THRESH * eval[5])  return -2.0f;
  if (verbose > 1)
    fprintf(stderr, "eigs: %g %g %g %g %g %g\n", eval[0], eval[1], eval[2], eval[3], eval[4], eval[5]);
  stability = eval[0] / eval[5];
  max_stability = eval[0];

  compute_alignxf(evec, eval, b, centroid, scale, alignxf);
  m1.xf = alignxf * m1.xf;

  if (do_affine) {
    for (size_t i = 0; i < pairs.size(); i++)
      pairs[i].p2 = alignxf * pairs[i].p2;
    compute_scale(pairs, alignxf, verbose);
    m1.xf = alignxf * m1.xf;
  }

  timestamp t4 = now();
  if (verbose > 1) {
    fprintf(stderr, "Computed xform in %.2f msec.\n",
            (t4 - t3) * 1000.0f);
  }

  if (point_vector) save_pairs(pairs, *point_vector);

  // Update CDFs, if necessary
  if (!m0.compute_cdf && !m0.compute_indexes &&
      !m1.compute_cdf && !m1.compute_indexes)
    return err;

  float einv[6];
  point centroids[2];
  float evecs[2][6][6], scales[2];
  compute_stable_probs(pairs, evecs[0], m0.evals, centroids[0], scales[0],
                              evecs[1], m1.evals, centroids[1], scales[1]);

  for (int idx = 0; idx < 2; idx++) {
    // we might only need to update one scan or the other
    if (!m[idx]->compute_cdf && !m[idx]->compute_indexes) continue;

    if (m[idx]->compute_cdf) m[idx]->cdf.resize(m[idx]->m->vertices.size());

    for (int i = 0; i < 6; i++)
      einv[i] = 1.0 / max(m[idx]->evals[i], EIG_THRESH * m[idx]->evals[5]);
    float Cinv[6][6];
    for (int i = 0; i < 6; i++) {
      float x[6];
      for (int j = 0; j < 6; j++)
        x[j] = (j == i) ? 1.0 : 0.0;
      eigmult(evecs[idx], einv, x, x);
      for (int j = 0; j < 6; j++)
        Cinv[i][j] = x[j];
    }

    vector< pair<float, int> > *probs = NULL;

    if (m[idx]->compute_indexes)
      probs = new vector< pair<float, int> >[6];

    // xform xfr = norm_xf(m[idx]->xf);
    size_t n = m[idx]->m->vertices.size();
    for (size_t i = 0; i < n; i++) {
      // we need to zero sampcdf[i] here in case
      // the weight is 0 and we continue
      if (m[idx]->compute_cdf) m[idx]->cdf[i] = 0.0;

      if (!m[idx]->weights[i]) continue;
      // if (edge(s1, i)) continue;
      // point p = xf1 * s1->vertices[i];
      point p = m[idx]->m->vertices[i];
      p -= centroids[idx];
      p *= scales[idx];
      // vec n = xf1r * s1->normals[i];
      vec n = m[idx]->m->normals[i];
      assert(!isnan(n[0]));
      assert(!isnan(n[1]));
      assert(!isnan(n[2]));
      vec c = p CROSS n;

      if (m[idx]->compute_indexes) {
        probs[0].push_back(pair<float, int>(c[0], i));
        probs[1].push_back(pair<float, int>(c[1], i));
        probs[2].push_back(pair<float, int>(c[2], i));
        probs[3].push_back(pair<float, int>(n[0], i));
        probs[4].push_back(pair<float, int>(n[1], i));
        probs[5].push_back(pair<float, int>(n[2], i));
      }

      if (m[idx]->compute_cdf) {
        for (int j = 0; j < 6; j++) {
          float tmp = Cinv[j][0] * c[0] + Cinv[j][1] * c[1] +
            Cinv[j][2] * c[2] + Cinv[j][3] * n[0] +
            Cinv[j][4] * n[1] + Cinv[j][5] * n[2];
          if (j < 3)
            m[idx]->cdf[i] += tmp * c[j];
          else
            m[idx]->cdf[i] += tmp * n[j-3];
        }
        m[idx]->cdf[i] *= m[idx]->weights[i];
      }
    }

    if (m[idx]->compute_indexes) {
      for (int i = 0; i < 6; i++) {
        m[idx]->indexes[i].clear();
        sort(probs[i].begin(), probs[i].end());

        // reverse the probabilities
        for (int j = (int) probs[i].size() - 1; j >= 0; j--)
          m[idx]->indexes[i].push_back(probs[i][j].second);
      }
      delete[] probs;
    }

    if (m[idx]->save_pdf) {
      m[idx]->pdf.resize(n);
      for (size_t i = 0; i < n; i++)
        m[idx]->pdf[i] = m[idx]->cdf[i];
    }

    for (size_t i = 1; i < n; i++)
      m[idx]->cdf[i] += m[idx]->cdf[i-1];
    if (!m[idx]->cdf[n - 1]) {
      if (verbose)
        fprintf(stderr, "No overlap.\n");
      return -1.0f;
    }

    float cscale = 1.0f / m[idx]->cdf[n - 1];
    for (size_t i = 0; i < n - 1; i++)
      m[idx]->cdf[i] *= cscale;
    m[idx]->cdf[n - 1] = 1.0f;
  }

  timestamp t5 = now();
  if (verbose > 1) {
    fprintf(stderr, "Updated CDFs in %.2f msec.\n",
            (t5-t4) * 1000.0f);
  }

  return err;
}


// Do ICP.  Aligns mesh s2 to s1, updating xf2 with the new transform.
// Returns alignment error, or -1 on failure
extern float ICP(ICP_struct &m0, ICP_struct &m1, float &stability, float &max_stability,
                 float maxdist, int verbose, vector<point> *point_pairs,
                 vector< std::pair<int, int> > &ipairs, bool use_indexes, bool do_affine) {
  ICP_struct *m[2] = { &m0, &m1 };
  size_t nv[2];

  // Make sure we have everything precomputed
  for (int i = 0; i < 2; i++) {
    nv[i] = m[i]->m->vertices.size();
    m[i]->m->need_normals();
    m[i]->m->need_curvatures();
    assert(m[i]->m->isedge.size == (int) nv[i]);
    // m[i]->m->need_edges();
    // m[i]->m->need_neighbors();
    // m[i]->m->need_adjacentfaces();
  }

  // Make sure weights, maxdist are filled in
  bool need_overlaps = (m0.weights.size() != nv[0]) || (m1.weights.size() != nv[1]);
  if (need_overlaps)
    compute_overlaps(m0, m1, maxdist, verbose);

  timestamp t = now();

  // Compute initial CDFs
  for (int idx = 0; idx < 2; idx++) {
    m[idx]->cdf.resize(nv[idx]);
    m[idx]->cdf[0] = m[idx]->weights[0];
    for (size_t i = 1; i < nv[idx]; i++)
      m[idx]->cdf[i] = m[idx]->weights[i] + m[idx]->cdf[i-1];
    if (!m[idx]->cdf[nv[idx] - 1]) {
      if (verbose)
        fprintf(stderr, "No overlap.\n");
      return -1.0f;
    }
    float cscale = 1.0f / m[idx]->cdf[nv[idx] - 1];
    for (size_t i = 0; i < nv[idx] - 1; i++)
      m[idx]->cdf[i] *= cscale;
    m[idx]->cdf[nv[idx] - 1] = 1.0f;
  }

  // Do an iteration, updating CDFs
  float incr = 4.0f / DESIRED_PAIRS;
  m0.compute_cdf = m1.compute_cdf = true;
  float err = ICP_iter(m0, m1, maxdist, verbose, stability, max_stability,
                       incr, point_pairs, ipairs, false, false);
  m0.compute_cdf = m1.compute_cdf = false;
  m0.compute_indexes = m1.compute_indexes = false;

  if (verbose > 1) {
    timestamp tnow = now();
    fprintf(stderr, "Time for this iteration: %.2f msec.\n\n",
            (tnow-t) * 1000.0f);
    t = tnow;
  }
  if (err < 0.0f)
    return err;

  bool do_affine_now = false;
  int iters = 1;
  vector<int> err_delta_history(TERM_HIST);
  while (iters++ < MAX_ITERS) {
    float lasterr = err;
    if (verbose > 1)
      fprintf(stderr, "Using incr = %f\n", incr);

    err = ICP_iter(m0, m1, maxdist, verbose, stability, max_stability,
                   incr, NULL, ipairs, use_indexes, do_affine_now);

    if (verbose > 1) {
      timestamp tnow = now();
      fprintf(stderr, "Time for this iteration: %.2f msec.\n\n",
              (tnow-t) * 1000.0f);
      t = tnow;
    }
    if (err < 0.0f)
      return err;

    // Check whether the error's been going up or down lately.
    // Specifically, we break out if error has gone up in
    // TERM_THRESH out of the last TERM_HIST iterations.
    for (int i = 0; i < TERM_HIST - 1; i++)
      err_delta_history[i] = err_delta_history[i+1];
    err_delta_history[TERM_HIST - 1] = (err >= lasterr);
    int nincreases = 0;
    for (int i = 0; i < TERM_HIST; i++)
      nincreases += err_delta_history[i];
    if (nincreases >= TERM_THRESH) {
      if (!do_affine || do_affine_now)
        break;
      err_delta_history.clear();
      err_delta_history.resize(TERM_HIST);
      do_affine_now = true;
    }
  }
  if (verbose > 1)
    fprintf(stderr, "Did %d iterations\n\n", iters);

#if 1
  // One final iteration at a higher sampling rate...
  if (verbose > 1)
    fprintf(stderr, "Last iteration...\n");
  incr *= 0.25f;  // 4 times as many pairs as last time
  if (verbose > 1)
    fprintf(stderr, "Using incr = %f\n", incr);
  if (need_overlaps) {
    compute_overlaps(m0, m1, maxdist, verbose);
  }

  err = ICP_iter(m0, m1, maxdist, verbose, stability, max_stability, incr,
                 point_pairs, ipairs, use_indexes, do_affine_now);

  if (verbose > 1) {
    timestamp tnow = now();
    fprintf(stderr, "Time for this iteration: %.2f msec.\n\n",
            (tnow-t) * 1000.0f);
    t = tnow;
  }
#endif
  return err;
}


// Minimalistic ICP
float ICP(EdgeMesh *s1, EdgeMesh *s2, const xform &xf1, xform &xf2) {
  KDtree *kd1 = new KDtree(s1->vertices);
  KDtree *kd2 = new KDtree(s2->vertices);
  vector< pair<int, int> > ipairs;
  float stability, max_stability;

  ICP_struct m1(s1, xf1, kd1, false);
  ICP_struct m2(s2, xf2, kd2, false);
  float icperr = ICP(m1, m2, stability, max_stability, 0.0f, 0, NULL, ipairs, true, false);

  delete kd2;
  delete kd1;
  return icperr;
}

#if 1
// compute stable sampling of npts points
//
// need to fix this
void stable_sample(EdgeMesh *m, int npts, vector<int> &points) {
  m->need_normals();
  int nv = m->vertices.size();

  // Sample some points randomly
  vector<PtPair> p;
  for (int i = 0; i < DESIRED_PAIRS; i++) {
    int ind = int(drand48() * nv);
    p.push_back(PtPair(m->vertices[ind], m->vertices[ind],
                       m->normals[ind],  m->normals[ind]));
  }

  // Compute covariance matrix
  float evec[6][6], eval[6], b[6], scale, err;
  point centroid;
  compute_ICPmatrix(p, evec, eval, b, centroid, scale, err);

  // Compute (clamped) inverse of covariance matrix
  float Cinv[6][6], einv[6];
  for (int i = 0; i < 6; i++)
    einv[i] = 1.0 / max(eval[i], EIG_THRESH * eval[5]);
  for (int i = 0; i < 6; i++) {
    float x[6];
    for (int j = 0; j < 6; j++)
      x[j] = (j == i) ? 1.0 : 0.0;
    eigmult(evec, einv, x, x);
    for (int j = 0; j < 6; j++)
      Cinv[i][j] = x[j];
  }

  // Compute sampling PDF
  vector<float> sampcdf(nv);
  for (int i = 0; i < nv; i++) {
    point p = m->vertices[i];
    p -= centroid;
    p *= scale;
    vec n = m->normals[i];
    vec c = p CROSS n;
    for (int j = 0; j < 6; j++) {
      float tmp = Cinv[j][0] * c[0] + Cinv[j][1] * c[1] +
        Cinv[j][2] * c[2] + Cinv[j][3] * n[0] +
        Cinv[j][4] * n[1] + Cinv[j][5] * n[2];
      if (j < 3)
        sampcdf[i] += tmp * c[j];
      else
        sampcdf[i] += tmp * n[j-3];
    }
  }

  // Accumulate and normalize to get CDF
  for (int i = 1; i < nv; i++)
    sampcdf[i] += sampcdf[i-1];
  float cscale = 1.0f / sampcdf[nv-1];
  for (int i = 0; i < nv-1; i++)
    sampcdf[i] *= cscale;
  sampcdf[nv-1] = 1.0f;

  // Sample according to CDF
  points.clear();
  points.reserve(npts);
  float incr = 2.0f / npts;
  float cval = 0.0f;
  int i = 0;
  while (1) {
    cval += incr * drand48();
    if (cval >= 1.0f)
      break;
    while (sampcdf[i] <= cval)
      i++;
    cval = sampcdf[i];
    points.push_back(i);
  }
}
#endif
