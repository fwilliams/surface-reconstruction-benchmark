#ifndef GLOBAL_REG_H
#define GLOBAL_REG_H

/*
Copyright 2007 Benedict Brown
Princeton University

global_reg.h
Prototypes, etc. for global range scan alignment

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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>

#include <vector>

#include "XForm.h"
#include "EdgeMesh.h"
#include "tps.h"

using namespace std;

struct corr_t {
  int vnum, tgt, stable;
  point pnt;

  corr_t(void) : vnum(-1), tgt(-1) {}
  corr_t(int v, int t, point &p, int s) :
    vnum(v), tgt(t), stable(s), pnt(p) {}
};
typedef vector<corr_t> corr_vector;

struct spring_t{
  int tgt;
  float len;
  float wgt;
  bool stable;

  spring_t(void) : tgt(-1), len(0), wgt(0), stable(false) {}
  spring_t(int t, float l, float w, bool s) :
    tgt(t), len(l), wgt(w), stable(s) {}
};
typedef vector<spring_t> spring_vector;


const char *get_mesh_name(const char *plyfilename);
float point_error(int pnum, const spring_vector &springs, bool skip_unstable,
                  const vector<point> &targets, const vector<bool> &use_points);
void set_errors(vector<float> &errors, const vector<spring_vector> &springs,
                const vector<point> &targets, const vector<bool> &use_points);
void compute_rigid_alignment(const farr &x, const farr &y, farr &R, point &T);
void write_cube_vertex(point *v, point &o, float w);
void write_cube_face(vector<TriMesh::Face> &faces, int f, int v);

// process command-line options for global_reg
// the prefix for mesh names is specified stand-alone
// and is required.  The global points file is specified on stdin.
struct opts_t {
  // bool rigid; // true for rigid-body, false for TPS alignment
  bool write_individual_points; // should meshes be written for individual points for debugging

  char *rigid_prefix; // what to prepend to mesh names for rigid alignment
  char *nonrigid_prefix; // what to prepend to mesh names for rigid alignment
  char *affine_prefix; // what to prepend to mesh names for rigid alignment

  float max_allowed_err, min_allowed_stability, min_allowed_max_stability, max_allowed_divergence;
  float min_target_dist2;

  int nthreads; // number of alignment threads to use

  FILE *disable_pts_file; // points to not use
  FILE *pts_file;         // which points we ultimately used

  // should an accompanying xf be read for each mesh, and should rigid/affine alignments
  // write an xf instead of a ply?
  bool read_xf, write_xf;


//???????????????????????????????
  opts_t(int argc = 0, char **argv = NULL) : write_individual_points(false),
       rigid_prefix(NULL), nonrigid_prefix(NULL), affine_prefix(NULL),
       max_allowed_err(1000), min_allowed_stability(0.001f), min_allowed_max_stability(1),
       max_allowed_divergence(50), min_target_dist2(64), nthreads(1),
       disable_pts_file(NULL), pts_file(NULL), read_xf(false), write_xf(false) {
    // we rely on default initialization of everything to false/NULL
    for (int i = 1; i < argc; i++) {
      if (!strcmp(argv[i], "-r")) {
        // Rigid alignment; the argument is the prefix to use
        rigid_prefix = argv[++i];
        fprintf(stderr, "Performing rigid alignments, with prefix \"%s\"\n", rigid_prefix);
      } else if (!strcmp(argv[i], "-nr")) {
        // Nonrigid alignment
        nonrigid_prefix = argv[++i];
        fprintf(stderr, "Performing nonrigid alignments, with prefix \"%s\"\n", nonrigid_prefix);
      } else if (!strcmp(argv[i], "-aff")) {
        // Affine alignment; implemented using only the affine
        // portion of the thin-plate spline
        affine_prefix = argv[++i];
        fprintf(stderr, "Performing affine alignments, with prefix \"%s\"\n", affine_prefix);
      } else if (!strcmp(argv[i], "-wp")) {
        // Specifies that meshes should be written
        // for each individual point; this is useful
        // for debugging which points went bad, but
        // it takes a long time to write these out if there
        // are many points.
        write_individual_points = true;
      } else if (!strcmp(argv[i], "-pf")) {
        // The list of points should be written
        // to the specified file; this can be
        // used to specifically disable points later
        assert(i < (argc - 1));
        i++;
        assert(pts_file = fopen(argv[i], "w"));
      } else if (!strcmp(argv[i], "-dpf")) {
        // argument is the list of points
        // which should be disabled for alignment
        assert(i < (argc - 1));
        i++;
        assert(disable_pts_file = fopen(argv[i], "r"));
      } else if (!strcmp(argv[i], "-min_target_dist")) {
        // minimum allowed distance between global features
        assert(i < (argc - 1));
        i++;
        min_target_dist2 = sqr(atof(argv[i]));
      } else if (!strcmp(argv[i], "-nthreads")) {
        // multithreaded warping/alignment
        // argument is number of threads to use
        assert(i < (argc - 1));
        nthreads = atoi(argv[++i]);
        fprintf(stderr, "Using %d threads.\n", nthreads);
        assert(nthreads >= 1);
      } else if (!strcmp(argv[i], "-read_xf")) {
        // an accompanying .xf should be read for
        // each mesh if it exists
        read_xf = true;
      } else if (!strcmp(argv[i], "-write_xf")) {
        // rigid and affine alignments should
        // write an xf rather than a ply
        write_xf = true;
      } else {
        // There should be no naked parameters
        assert(0);
#if 0
        assert(!have_prefix);
        have_prefix = true;
        prefix = argv[i];
        fprintf(stderr, "Assigned prefix: %s\n", prefix);
#endif
      }
    }

    // assert(have_prefix);
    if (!(rigid_prefix || nonrigid_prefix || affine_prefix || write_individual_points || pts_file)) {
      fprintf(stderr, "No options which produce persistent output have been specified.  Exiting.\n");
      exit(0);
    }
  }

  ~opts_t() {
    if (disable_pts_file) fclose (disable_pts_file);
    if (pts_file) fclose(pts_file);
  }
};

// This class hold the state necessary to run multiple warps in
// parallel.  It spawns multiple threads, which each query the shared
// variable counter (protected by counter_mutex) to get the index of
// the next mesh to warp.  The only other shared state that is not
// read-only is points_mesh.  All accesses to this should be protected
// by locking mesh_mutex.
//
// Note that the constructor does all the work, so threaded_alignments
// should be created an destroyed only.
void align_scan(const opts_t &opts, const char *mesh_name, const corr_vector &corrs,
                const vector<point> &targets, const vector<bool> &use_points, const vector<float> &confidence,
                TriMesh *points_mesh, const vector<TriMesh::Face> &all_faces, pthread_mutex_t *mutex);
class threaded_alignment {
  mutable pthread_mutex_t counter_mutex, mesh_mutex;
  volatile size_t counter;
  vector<pthread_t> threads;

  // This stuff really doesn't need to be in a struct, it
  // could just be member data.  But this leads to a marginally
  // less verbose initialization using memcpy.
  struct align_scan_struct {
    const opts_t opts;
    const vector<char *> mesh_names;
    const vector<corr_vector> corrs;
    const vector<point> targets;
    const vector<bool> use_points;
    const vector<float> confidence;
    TriMesh *points_mesh;
    const vector<TriMesh::Face> all_faces;
    threaded_alignment *ta;

    // align_scan_struct(const opts_t &opts, const vector<char *> &mesh_names, const vector<corr_vector> &corrs, const vector<point> &targets, const vector<bool> &use_points, const vector<float> &confidence, const vector<TriMesh::Face> &all_faces)
    // : opts(opts), mesh_names(mesh_names), corrs(corrs), targets(targets), use_points(use_points), confidence(confidence), all_faces(all_faces)
    // { }
  } data;

  static void *alignment_thread(void *data) {
    size_t i;
    align_scan_struct *d = (align_scan_struct *) data;
    while ((i = d->ta->count()) < d->mesh_names.size()) {
      align_scan(d->opts, d->mesh_names[i], d->corrs[i],
                  d->targets, d->use_points, d->confidence,
                  d->points_mesh, d->all_faces, &d->ta->mesh_mutex);
    }
    return NULL;
  }

  // lock a mutex safely
  static bool lock(pthread_mutex_t *m) {
    bool  success;
    success = !pthread_mutex_lock(m);
    assert(success);
    return success;
  }

  // unlock a mutex
  static void unlock(pthread_mutex_t *m) { pthread_mutex_unlock(m); }

  int count() { lock(&counter_mutex); size_t rc = counter++; unlock(&counter_mutex); return rc; }

 public:
  threaded_alignment(const opts_t &opts, const vector<char *> &mesh_names,
                     const vector<corr_vector> &corrs, const vector<point> &targets,
                     const vector<bool> &use_points, const vector<float> &confidence,
                     TriMesh *points_mesh, const vector<TriMesh::Face> &all_faces) : counter(0) {
    pthread_mutex_init(&counter_mutex, NULL);
    pthread_mutex_init(&mesh_mutex, NULL);

    align_scan_struct d = { opts, mesh_names, corrs, targets, use_points,
             confidence, points_mesh, all_faces, this };
    memcpy(&data, &d, sizeof(align_scan_struct));

    threads.resize(opts.nthreads);
    for (int i = 0; i < opts.nthreads; i++) {
      pthread_create(&threads[i], NULL, alignment_thread, (void *) &data);
    }
  }

  ~threaded_alignment() {
    for (size_t i = 0; i < threads.size(); i++)
      pthread_join(threads[i], NULL);
    pthread_mutex_destroy(&counter_mutex);
    pthread_mutex_destroy(&mesh_mutex);
  }
};

#endif /* GLOBAL_REG_H */
