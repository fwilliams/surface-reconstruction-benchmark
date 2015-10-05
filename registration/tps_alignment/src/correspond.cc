/*
Copyright 2007 Benedict Brown
Princeton University

correspond.cc
Compute Locally Weighted ICP feature correspondences of all scans to a
specified one.

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <vector>

#include "timestamp.h"
#include "EdgeMesh.h"
#include "TriMesh_algo.h"
#include "ICP.h"

using namespace std;

// this groups all the command line options
// together, with the code to read them
// and a destructor to close any open files
struct opts_t {
  char *tgt_name;
  FILE *src_list;
  FILE *bbox_list;
  FILE *output;

  int   num_overlaps;
  float sample_rate;
  float bbox_dist;
  float falloff;
  float kd;
  bool  use_indexes;
  bool  nolocal;
  bool  affine;
  bool  read_xf;

  static void usage(char *progname = "correspond") {
    fprintf(stderr, "Usage:   %s tgt_mesh mesh_list sample_rate {options}\n"
            "Options: -num_overlaps %%d (minimum number of overlapping points to align two scans)\n"
            "         -kd %%f (maximum nearest neighbor distance for KD-tree)\n"
            "         -falloff %%f (std dev of Gaussian falloff for weighted ICP)\n"
            "         -bbox %%s (name of bbox cache file)\n"
            "         -max_bbox_dist %%f (maximum distance between two bboxes to consider scans overlapping)\n"
            "         -output %%s (file to redirect correspondences to; default is stdout)\n"
            "         -sort_stability (use the most stable points rather than sampling the stability pdf)\n"
            "         -nolocal (do not perform locally weighted ICP)\n"
            "         -affine  (do affine ICP alignment)\n"
            "         -read_xf (read a .xf file for each mesh\n", progname);
  }

  opts_t(int *argc, char **argv) :
    src_list(NULL), bbox_list(NULL), output(stdout), num_overlaps(500),
    sample_rate(-1), bbox_dist(5), falloff(20), kd(20),
    use_indexes(false), nolocal(false), affine(false), read_xf(false)
  {
    bool have_tgt = false, have_src = false, have_sample = false;

    for (int i = 1; i < *argc; i++) {
      if (argv[i][0] != '-') {
        // not an option, so it must be one of tgt_mesh, src_list or sample_rate
        if (!have_tgt) {
          have_tgt = true;
          tgt_name = argv[i];
          fprintf(stderr, "tgt_name: %s\n", tgt_name);
        } else if (!have_src) {
          have_src = true;
          src_list = fopen(argv[i], "r");
          if (!src_list) {
            fprintf(stderr, "Could not open mesh list %s\n", argv[i]);
            exit(-1);
          }
        } else if (!have_sample) {
          have_sample = true;
          sample_rate = atof(argv[i]);
          fprintf(stderr, "sample_rate: %f (%s)\n", sample_rate, argv[i]);
          if (sample_rate <= 0) {
            fprintf(stderr, "Sample rate (%f) must be greater than 0.\n", sample_rate);
            exit(-1);
          }

          if (sample_rate > 1) {
            fprintf(stderr, "Sample rate (%f) must be at most 1.\n", sample_rate);
            exit(-1);
          }
        } else {
          fprintf(stderr, "Unexpected argument %d: %s\n", i, argv[i]);
          exit(-1);
        }
      } else if (!strcmp(argv[i], "-num_overlaps")) {
        num_overlaps = atoi(argv[++i]);
      } else if (!strcmp(argv[i], "-falloff")) {
        falloff = atof(argv[++i]);
      } else if (!strcmp(argv[i], "-kd")) {
        kd = atof(argv[++i]);
      } else if (!strcmp(argv[i], "-bbox")) {
        bbox_list = fopen(argv[++i], "r");
        if (!bbox_list) {
          fprintf(stderr, "Could not open bounding box cache file %s.\n", argv[i]);
          exit(-1);
        }
      } else if (!strcmp(argv[i], "-max_bbox_dist")) {
        bbox_dist = atof(argv[++i]);
      } else if (!strcmp(argv[i], "-output")) {
        output = fopen(argv[++i], "w");
        setbuf(output, NULL); // it's handy not to buffer this stuff
        if (!output) {
          fprintf(stderr, "Could not open output file %s.\n", argv[i]);
          exit(-1);
        }
      } else if (!strcmp(argv[i], "-sort_stability")) {
        use_indexes = true;
      } else if (!strcmp(argv[i], "-nolocal")) {
        nolocal = true;
      } else if (!strcmp(argv[i], "-affine")) {
        affine = true;
      } else if (!strcmp(argv[i], "-read_xf")) {
        read_xf = true;
      } else {
        fprintf(stderr, "Unknown option %s\n", argv[i]);
        exit(-1);
      }
    }
  }

  ~opts_t() {
    if (src_list  != NULL) { fclose(src_list);  src_list  = NULL; }
    if (bbox_list != NULL) { fclose(bbox_list); bbox_list = NULL; }
    if (output  != stdout) { fclose(output);    output  = stdout; }
  }


  // fallof default is good for David, should be 5 for penguin
  // use_indexes should be true for FUR, false for most other things
private:
  opts_t(void) {}
};


struct corr {
  int tgt; // the number of the mesh to which we match
  point p; // the point we match to
  float stability, max_stability, err;
};

// #define SQ(x) ((x) * (x))


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
point closest_on_face(const EdgeMesh *mesh, int i, const point &p) {
  const EdgeMesh::Face &f = mesh->faces[i];
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
int closest_pt(const EdgeMesh *mesh, const KDtree *kd, float mdist, const point &p, point &out) {
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

// get the next file name from the list
// we assume that filenames are whitespace deliminated
// and themselves contain no whitespaces
//
// one would hope that fscanf(flist, " %s ", fname); would
// do exactly this but I guess it wasn't working
static void read_fname(char *fname, FILE *flist) {
  int cpos = 0, c = 0;

  while (!feof(flist) && isspace(c = getc(flist)));
  if (!feof(flist)) fname[cpos++] = c;
  while (!feof(flist) && !isspace(c = getc(flist))) fname[cpos++] = (char) c;
  fname[cpos] = '\0';
}

// read in a mesh that may be a .pre or a .ply
// if it's a .ply, do some normal diffusion
// note that it's much better to preprocess the meshes--
// it's faster, and you can specify the amounts of diffusion
// to geometry, normals and curvature
static EdgeMesh *read_mesh(char *name, bool read_xf) {
  EdgeMesh *mesh = EdgeMesh::read_preprocessed(name);
  if (!mesh) {
    assert(mesh = EdgeMesh::read(name)); putc('\n', stderr);
    if (mesh->normals.size() != mesh->vertices.size()) {
      mesh->normals.clear();
      mesh->need_normals();
      // diffuse_normals(mesh, .5 * mesh->feature_size());
    }
  }

  assert(mesh);

  if (read_xf) {
    char xfname[1024];
    strcpy(xfname, name);
    strcpy(xfname + strlen(xfname) - 3, "xf");
    xform xf;
    xf.read(xfname);
    fprintf(stderr, "Read input transform:\n");
    std::cerr << xf;
    xform xfr = rot_only(xf);
    for (size_t i = 0; i < mesh->vertices.size(); i++) {
      mesh->vertices[i] = xf * mesh->vertices[i];
      mesh->normals[i] = xfr * mesh->normals[i];
    }
  }

  mesh->need_bbox();
  mesh->need_edges();
  mesh->need_curvatures();

  return mesh;
}

// helper functions and definitions for writing out the correspondences
typedef enum {
  FEATURE_ICP_UNSTABLE = -1,
  FEATURE_NO_OVERLAP = -3,
  MESH_NO_OVERLAP = -4,
  BBOX_PRUNED = -5,
  MESH_ICP_UNSTABLE = -6,
} match_t;

int main(int argc, char *argv[]) {
  // assert(argc >= 8);

  // seed the random number generator
  srand48(time(NULL));
  float reading_time = 0, kd_time = 0, icp_time = 0;

  opts_t opts(&argc, argv);

  timestamp tmp_time = now();

  // Read in the target mesh
  EdgeMesh *tgt_mesh = read_mesh(opts.tgt_name, opts.read_xf);
  EdgeMesh::set_verbose(0);

  // select sample vertices on source mesh, disallowing boundary vertices

  reading_time += now() - tmp_time;

  int num_samples = (int) (tgt_mesh->vertices.size() * opts.sample_rate + 0.5f);
  if (num_samples == 0) num_samples++;
  vector<int> samples;
  vector<int> mean_samples;
  vector<point> mean_pos;

  vector<size_t> non_edge_samples;
  for (size_t i = 0; i < tgt_mesh->vertices.size(); i++)
    if (!tgt_mesh->isedge(i)) non_edge_samples.push_back(i);

  if (non_edge_samples.size() > 0) {
    fprintf(stderr, "Selecting %d samples on mesh with %d vertices.\n", num_samples, (int) tgt_mesh->vertices.size());

    stable_sample(tgt_mesh, num_samples / 2, samples);

    fprintf(stderr, "Selected %d stable samples, now selecting uniform samples.\n", (int) samples.size());

    while ((int) samples.size() < num_samples) {
      size_t vnum = (size_t) (non_edge_samples.size() * drand48() + 0.5f);
      if (vnum >= non_edge_samples.size()) continue;
      vnum = non_edge_samples[vnum];
      samples.push_back(vnum);
      // if (vnum == tgt_mesh->vertices.size()) continue;
      // if (!tgt_mesh->isedge(vnum)) samples.push_back(vnum);
    }

    fprintf(stderr, "Selected %d samples.\n", num_samples);
  }

  if (samples.size() == 0) {
    fprintf(stderr, "There are no samples, aborting.\n");
    fprintf(opts.output, "0 %s\n", opts.tgt_name);
    exit(0);
  }
  non_edge_samples.clear();

  vector< vector<corr> > corrs(num_samples);

  tmp_time = now();
  KDtree *tgt_kd = new KDtree(tgt_mesh->vertices);
  kd_time += now() - tmp_time;

  ICP_struct tgt_icp(tgt_mesh, xform::identity(), tgt_kd, false);

  // read in the list of meshes to align to
  int tgt_num = -1;
  vector<char *> mesh_names;
  while (!feof(opts.src_list)) {
    char fname[1024];
    read_fname(fname, opts.src_list);
    if (fname[0] == '\0') break;

    if (!strcmp(opts.tgt_name, fname))
      tgt_num = (int) mesh_names.size();
    mesh_names.push_back(strdup(fname));
  }

  // now we have features selected on the target mesh, so load
  // each source mesh in turn and find the correspondences on it
  for (int src_num = 0; src_num < (int) mesh_names.size(); src_num++) {
    char *fname = mesh_names[src_num];

    bool skip_mesh = false;
    // unsigned int nnormals = 0;
    float bbox[6];

    if (opts.bbox_list) fread(bbox, sizeof(float), 6, opts.bbox_list); // read in bounding box

    // if target and source are the same mesh
    if (src_num == tgt_num) {
      // this is a sanity check that our bbox cache
      // matches what we actually computed
      if (opts.bbox_list) {
        for (int i = 0; i < 3; i++) {
          if (tgt_mesh->bbox.min[i] != bbox[i])
            fprintf(stderr, "Difference in bbox.min[%d]: %g -> %g\n",
                    i, tgt_mesh->bbox.min[i], bbox[i]);
          if (tgt_mesh->bbox.max[i] != bbox[i + 3])
            fprintf(stderr, "Difference in bbox.max[%d]: %g -> %g\n",
                    i, tgt_mesh->bbox.max[i], bbox[i + 3]);
        }
        assert(tgt_mesh->bbox.min[0] == bbox[0]);
        assert(tgt_mesh->bbox.min[1] == bbox[1]);
        assert(tgt_mesh->bbox.min[2] == bbox[2]);
        assert(tgt_mesh->bbox.max[0] == bbox[3]);
        assert(tgt_mesh->bbox.max[1] == bbox[4]);
        assert(tgt_mesh->bbox.max[2] == bbox[5]);
      }

      // since the two meshes are the same, each feature
      // necessarily corresponds to itself; we can just write
      // these out
      for (int i = 0; i < num_samples; i++) {
        corr c = { tgt_num, tgt_mesh->vertices[samples[i]], 1, 10000, 0 };
        corrs[i].push_back(c);
      }
      continue;
    }

    fprintf(stderr, "%s\n", fname);

    if (opts.bbox_list) {
      // check bounding boxes for overlap info
      if (tgt_mesh->bbox.max[0] < bbox[0] - opts.bbox_dist) skip_mesh = true;
      if (tgt_mesh->bbox.max[1] < bbox[1] - opts.bbox_dist) skip_mesh = true;
      if (tgt_mesh->bbox.max[2] < bbox[2] - opts.bbox_dist) skip_mesh = true;
      if (tgt_mesh->bbox.min[0] > bbox[3] + opts.bbox_dist) skip_mesh = true;
      if (tgt_mesh->bbox.min[1] > bbox[4] + opts.bbox_dist) skip_mesh = true;
      if (tgt_mesh->bbox.min[2] > bbox[5] + opts.bbox_dist) skip_mesh = true;
    }

    if (skip_mesh) {
      continue;
    }

    // load in the next source mesh
    tmp_time = now();
    EdgeMesh *src_mesh = read_mesh(fname, opts.read_xf);
    reading_time += now() - tmp_time;

    tmp_time = now();
    KDtree *src_kd = new KDtree(src_mesh->vertices);
    kd_time += now() - tmp_time;

    ICP_struct src_icp(src_mesh, xform::identity(), src_kd, false);

    float maxdist = 15.0f;
    int verbose = 0;

    // Unweighted ICP to initialize things.  Fills in weights based
    // on overlaps.
    tmp_time = now();
    vector< pair<int, int> > ipairs;
    float stability, max_stability;

    vector<bool> src_overlaps(src_mesh->vertices.size());
    vector<bool> tgt_overlaps(tgt_mesh->vertices.size());

    tgt_icp.weights.clear();
    float icpok = ICP(tgt_icp, src_icp, stability, max_stability, 0,
                      verbose, NULL, ipairs, opts.use_indexes, opts.affine);
    xform original_xf = src_icp.xf;

    ipairs.clear();
    if (icpok < 0) {
      delete src_mesh;
      delete src_kd;
      continue;
    }

    // set a bitmap of which points overlap in the two meshes, and count
    // the number of overlaps
    int num_overlaps_tgt = 0, num_overlaps_src = 0;
    for (size_t i = 0; i < tgt_overlaps.size(); i++)
      if (tgt_icp.weights[i]) { tgt_overlaps[i] = true; num_overlaps_tgt++; }
    for (size_t i = 0; i < src_overlaps.size(); i++)
      if (src_icp.weights[i]) { src_overlaps[i] = true; num_overlaps_src++; }
    fprintf(stderr, "num_overlaps_tgt = %d, num_overlaps_src = %d\n", num_overlaps_tgt, num_overlaps_src);

    // don't try and find correspondences if the overlap is too small
    if (max(num_overlaps_src, num_overlaps_tgt) < opts.num_overlaps) {
      delete src_mesh;
      delete src_kd;
      continue;
    }

    // otherwise do try and find correspondences for each feature
    for (int i = 0; i < num_samples; i++) {
      // this feature doesn't overlap the other mesh, so skip it
      if (!tgt_overlaps[samples[i]]) {
        // write_bad_correspondence(opts.output, FEATURE_NO_OVERLAP);
        continue;
      }

      // compute weighted pdf for this feature on target mesh
      // (falloff by squared distance from feature)
      const point &p = tgt_mesh->vertices[samples[i]];
      for (size_t j = 0; j < tgt_mesh->vertices.size(); j++) {
        if (tgt_overlaps[j]) {
          float sqdist = dist2(p, tgt_mesh->vertices[j]);
          tgt_icp.weights[j] = 1.0f / (sqr(opts.falloff) + sqdist);
        } else tgt_icp.weights[j] = 0;
      }

      // compute weighted pdf for this feature on source mesh
      // (squared distance from feature after rigid-body transform applied)
      point q = inv(src_icp.xf) * p;
      for (size_t j = 0; j < src_mesh->vertices.size(); j++) {
        if (src_overlaps[j]) {
          float sqdist = dist2(q, src_mesh->vertices[j]);
          src_icp.weights[j] = 1.0f / (sqr(opts.falloff) + sqdist);
        } else src_icp.weights[j] = 0;
      }

      // Do weighted ICP
      src_icp.xf = original_xf;
      float icp_err = icpok;
      if (!opts.nolocal) {
        vector< pair<int, int> > ipairs;
        // here use_indexes should always be false
        icp_err = ICP(tgt_icp, src_icp, stability, max_stability,
                      maxdist, verbose, NULL, ipairs, false, opts.affine);
        // bool icp_stable = (icp_err >= 0.0f);
        ipairs.clear();
      }

      q = inv(src_icp.xf) * p;

      point r;
      int match_num = closest_pt(src_mesh, src_kd, sqr(opts.kd), q, r);

      if (match_num > 0) {
        // we found a match
        corr c = { src_num, r, stability, max_stability, icp_err };
        corrs[i].push_back(c);
      } else {
        // we didn't find a match, but we might as well record the stability, etc.
        // the factor that match_num < 0 tells the global_reg that this correpsondence is bad
      }
    }
    icp_time += now() - tmp_time;

    // putc('\n', opts.output); // only useful for old style correspond output

    tmp_time = now();
    delete src_kd;
    delete src_mesh;
    kd_time += now() - tmp_time;
  }

  // print out correspondences to opts.output
  int real_num_points = 0;
  for (int i = 0; i < num_samples; i++)
    if (corrs[i].size() > 1) real_num_points++;

  // we often work with .pre file, and need to
  // replace the extension with .ply for output
  char tgt_name[1024];
  strcpy(tgt_name, opts.tgt_name);
  strcpy(tgt_name + strlen(tgt_name) - 3, "ply");
  fprintf(opts.output, "%d %s\n", real_num_points, tgt_name);
  for (int i = 0; i < num_samples; i++) {
    if (corrs[i].size() <= 1) continue;
    point p = tgt_mesh->vertices[samples[i]];
    fprintf(opts.output, "\n%f %f %f %d %d\n", p[0], p[1], p[2], samples[i], (int) corrs[i].size());
    for (size_t j = 0; j < corrs[i].size(); j++) {
      corr c = corrs[i][j];
      fprintf(opts.output, "%d %f %f %f %f %f %f\n", c.tgt, c.p[0], c.p[1], c.p[2], c.stability, c.max_stability, c.err);
    }
  }

  delete tgt_kd;
  delete tgt_mesh;

  fprintf(stderr, "Timing Info: %f reading, %f kd_tree, %f icp\n",
          reading_time, kd_time, icp_time);
}
