/*
Copyright 2007 Benedict Brown
Princeton University

preprocess.cc

Preprocess a .ply file, writing out different information (vertices,
smoothed normals, etc.)depending on the type into a .pre file.  .pre
file are a very simple endian-dependent dump of the data.

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
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <GL/glut.h>
#include <assert.h>

#include <vector>

#include <string.h>
#include <algorithm>
#include <utility>
#include <iostream>
#include <fstream>

#include <tnt_array2d.h>

#include "EdgeMesh.h"
#include "TriMesh_algo.h"
#include "ICP.h"

using namespace std;

// process command-line options for preprocess
struct opts_t {
  char *bbox_file;
  enum { BBOX, FULL, FACE, CURVE, COLOR } type;
  float dn, dc; // diffuse_normals and diffuse_curvature


  opts_t(int argc, char **argv) : bbox_file(NULL), type(BBOX), dn(2), dc(1) {
    bool have_type = false;

    fprintf(stderr, "argc == %d\n", argc);

    // we rely on default initialization of everything to false/NULL
    for (int i = 1; i < argc; i++) {
      fprintf(stderr, "matching option %d: %s\n", i, argv[i]);
      if (!strcmp(argv[i], "-bbox")) {
        have_type = true;
        assert(i < argc - 1);
        type = BBOX;
        bbox_file = argv[++i];
        fprintf(stderr, "bbox type, file = %s\n", bbox_file);
      } else if (!strcmp(argv[i], "-full")) {
        have_type = true;
        assert(i < argc - 1);
        type = FULL;
        bbox_file = argv[++i];
        fprintf(stderr, "full type, file = %s\n", bbox_file);
      } else if (!strcmp(argv[i], "-face")) {
        have_type = true;
        assert(i < argc - 1);
        type = FACE;
        bbox_file = argv[++i];
        fprintf(stderr, "face type, file = %s\n", bbox_file);
      } else if (!strcmp(argv[i], "-curve")) {
        have_type = true;
        assert(i < argc - 1);
        type = CURVE;
        bbox_file = argv[++i];
        fprintf(stderr, "curve type, file = %s\n", bbox_file);
      } else if (!strcmp(argv[i], "-color")) {
        have_type = true;
        assert(i < argc - 1);
        type = COLOR;
        bbox_file = argv[++i];
        fprintf(stderr, "color type, file = %s\n", bbox_file);
      } else if (!strcmp(argv[i], "-dn")) {
        assert(i < argc - 1);
        dn = atof(argv[++i]);
        fprintf(stderr, "normal smoothing: %f\n", dn);
      } else if (!strcmp(argv[i], "-dc")) {
        assert(i < argc - 1);
        dc = atof(argv[++i]);
        fprintf(stderr, "normal smoothing: %f\n", dc);
      }
    }

    assert(have_type);
  }
};

// Determine whether point i is an edge point
// For point clouds without faces, assume nothing is an edge
static inline bool edge(const EdgeMesh *s, int i) {
  return s->faces.size() && (s->adjacentfaces[i].size() != s->neighbors[i].size());
}


FILE *outfile = stdout;

void read_fname(char *fname, FILE *flist) {
  int cpos = 0, c = 0;

  while (!feof(flist) && isspace(c = getc(flist)));
  if (!feof(flist)) fname[cpos++] = c;
  while (!feof(flist) && !isspace(c = getc(flist))) fname[cpos++] = (char) c;
  fname[cpos] = '\0';
}

int main(int argc, char *argv[])
{
  assert(argc >= 1);

  srand48(0);

  opts_t opts(argc, argv);

  FILE *outfile = fopen(opts.bbox_file, "w");

  // this needs to be supported in opts_t
  FILE *fstable = NULL;

  if (opts.type == opts_t::COLOR) {
    for (int i = 0; i < 4; i++) putc('F', outfile);
    putc('C', outfile); putc('O', outfile); putc('L', outfile); putc('R', outfile);
  }

  while (!feof(stdin)) {
    char fname[1024];
    read_fname(fname, stdin);
    if (fname[0] == '\0') continue;

    fprintf(stderr, "%s\n", fname);

    srand48(0);
    EdgeMesh *mesh = EdgeMesh::read(fname); putc('\n', stderr);

    mesh->need_bbox();

    float min_color = 0, max_color = 1;
    if ((opts.type == opts_t::COLOR) && !mesh->colors.empty()) {
      vector<Color> &c = mesh->colors;
      min_color = min(min(c[0][0], c[0][1]), c[0][2]);
      max_color = max(max(c[0][0], c[0][1]), c[0][2]);
      for (size_t i = 1; i < c.size(); i++) {
        min_color = min(min(min(c[i][0], c[i][1]), c[i][2]), min_color);
        max_color = max(max(max(c[i][0], c[i][1]), c[i][2]), max_color);
      }
    }

    // example diffusion settings:
    // david: normals = 2, curve = 1
    // fur:   normals = 1, curve = 1

    if (fstable || opts.type == opts_t::FULL || opts.type == opts_t::FACE || opts.type == opts_t::CURVE) {
      mesh->normals.clear();
      mesh->need_normals();
      mesh->need_neighbors();
      mesh->need_adjacentfaces();
      if (opts.dn > 0)
        diffuse_normals(mesh, opts.dn * mesh->feature_size());
    }

    if (opts.type == opts_t::CURVE) {
      mesh->need_curvatures();
      if (opts.dc > 0)
        diffuse_curv(mesh, opts.dc * mesh->feature_size());
    }

    if (fstable) {
      KDtree *kd = new KDtree(mesh->vertices);
      vector< pair<int, int> > ipairs;
      float stability, max_stability;

      ICP_struct m0(mesh, xform::identity(), kd, false);
      ICP_struct m1(mesh, xform::identity(), kd, false);
      float icperr = ICP(m0, m1, stability, max_stability, 100.0f, 0, NULL, ipairs, false, false);
      delete kd;

      if (icperr < 0 || stability < 0.1f) {
        fprintf(stderr, "Unstable %f\n", icperr);
        delete mesh;
        read_fname(fname, stdin);
        continue;
      }

      fprintf(fstable, "%s\n", fname);
    }

    if (opts.type == opts_t::FULL || opts.type == opts_t::FACE || opts.type == opts_t::CURVE) {
      char vert_header[] = { 'F', 'F', 'F', 'F', 'V', 'E', 'R', 'T' };
      char face_header[] = { 'F', 'F', 'F', 'F', 'F', 'A', 'C', 'E' };
      char curv_header[] = { 'F', 'F', 'F', 'F', 'C', 'U', 'R', 'V' };
      read_fname(fname, stdin);
      FILE *premesh = fopen(fname, "w");
      assert(premesh);
      unsigned int vsize = mesh->vertices.size();
      unsigned int fsize = mesh->faces.size();
      assert(mesh->normals.size() == vsize);
      fwrite(opts.type == opts_t::FACE ? face_header : (opts.type == opts_t::CURVE ? curv_header : vert_header),
             sizeof(char), 8, premesh);
      fwrite(&vsize, sizeof(int), 1, premesh);
      if (opts.type == opts_t::FACE) fwrite(&fsize, sizeof(int), 1, premesh);
      fwrite(&mesh->vertices[0], sizeof(float), 3 * vsize, premesh);
      fwrite(&mesh->normals[0],  sizeof(float), 3 * vsize, premesh);

      if (opts.type == opts_t::CURVE) {
        // just write mean curvature, not principal curvatures
        for (unsigned int i = 0; i < vsize; i++) {
          float c = 0.5f * (mesh->curv1[i] + mesh->curv2[i]);
          fwrite(&c, sizeof(float), 1, premesh);
        }
      }

      int isedge[(vsize + 31) / 32];
      for (unsigned int i = 0; i < (vsize + 31) / 32; isedge[i++] = 0);
      for (unsigned int i = 0; i < vsize; i++) {
        int bit = (mesh->neighbors[i].size() != mesh->adjacentfaces[i].size());
        isedge[i / 32] |= (bit << (i % 32));
      }
      fwrite(isedge, 4, (vsize + 31) / 32, premesh);
      if (opts.type == opts_t::FACE) fwrite(&mesh->faces[0], sizeof(int), 3 * fsize, premesh);
      fclose(premesh);
    }

    float bbox[6];
    bbox[0] = mesh->bbox.min[0];
    bbox[1] = mesh->bbox.min[1];
    bbox[2] = mesh->bbox.min[2];
    bbox[3] = mesh->bbox.max[0];
    bbox[4] = mesh->bbox.max[1];
    bbox[5] = mesh->bbox.max[2];
    fwrite((const void *) bbox, 6, sizeof(float), outfile);
    if (opts.type == opts_t::COLOR) {
      fwrite((const void *) &min_color, 1, sizeof(float), outfile);
      fwrite((const void *) &max_color, 1, sizeof(float), outfile);
    }

    delete mesh;
  }

  fclose(outfile);
  if (fstable) fclose(fstable);
}
