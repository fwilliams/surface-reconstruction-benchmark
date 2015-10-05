#ifndef ICP_H
#define ICP_H
/*
Copyright 2007 Szymon Rusinkiewicz and Benedict Brown
Princeton University

ICP.h
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

#include "EdgeMesh.h"
#include "XForm.h"
#include "KDtree.h"

// groups all the data and settings for ICP
// which are associated with a single mesh
class ICP_struct {
 public:
  EdgeMesh *m;
  xform xf;
  const KDtree *kd;
  vector<float> weights;
  vector<float> pdf;
  vector<float> cdf;
  vector<int> indexes[6];
  float evals[6];

  bool save_pdf;
  bool compute_cdf;
  bool compute_indexes;

  ICP_struct(EdgeMesh *m_, const xform &xf_, const KDtree *kd_, bool save_pdf_) :
    m(m_), xf(xf_), kd(kd_), save_pdf(save_pdf_), compute_cdf(true), compute_indexes(true) {}

 private:
  ICP_struct(void);
};

// Determine which points on m1 and m2 overlap, filling in weights for each one
// Also fills in maxdist, if it is <= 0 on input
extern void compute_overlaps(ICP_struct &m0, ICP_struct &m1, float &maxdist, int verbose);

// Do ICP.  Aligns mesh s2 to s1, updating xf2 with the new transform.
// Returns alignment error, or -1 on failure.
// Pass in 0 for maxdist to figure it out...
// Pass in empty weights to figure it out...
extern float ICP(ICP_struct &m0, ICP_struct &m1, float &stability, float &max_stability,
                 float maxdist, int verbose, vector<point> *point_pairs,
                 vector< std::pair<int, int> > &ipairs, bool use_indexes, bool do_affine);
extern float ICP(EdgeMesh *s1, EdgeMesh *s2, const xform &xf1, xform &xf2);

extern void stable_sample(EdgeMesh *m, int npts, vector<int> &points);

#endif
