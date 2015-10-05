/*
Copyright 2007 Benedict Brown
Princeton University

EdgeMesh.cc
A subclass of TriMesh.cc which explicitly represents mesh boundary information.

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
#include "TriMesh_algo.h"

void EdgeMesh::need_edges(int bsize) {
  if (isedge.size != 0) return;

  need_neighbors();
  need_adjacentfaces();

  isedge.resize(vertices.size());
  assert(isedge.size <= (isedge.words * 32));
  memset(isedge.isedge, 0, isedge.words * sizeof(int));
  for (unsigned int i = 0; i < vertices.size(); i++) {
    int bit = (neighbors[i].size() != adjacentfaces[i].size());
    bit |= (neighbors[i].size() == 0);
    isedge.isedge[i / 32] |= (bit << (i % 32));
  }
 
  // if this is a point cloud without faces, assume nothing is an edge
  if (!faces.size()) return;
 
  if (bsize > 0) {
    vector<int> frontier(vertices.size()), next(vertices.size());

    for (unsigned int i = 0; i < vertices.size(); i++)
      frontier[i] = isedge(i);

    for (int i = 0; i < bsize; i++) {
      for (unsigned int j = 0; j < vertices.size(); j++) {
        if (!frontier[j]) continue;
        for (unsigned int k = 0; k < neighbors[j].size(); k++) {
          int l = neighbors[j][k];
          next[l] = !(frontier[l] | isedge(l));
        }
      }

      for (unsigned int j = 0; j < vertices.size(); j++) {
        assert((int) (j / 32) < isedge.words);
        int bit = frontier[j];
        isedge.isedge[j / 32] |= (bit << (j % 32));
        frontier[j] = next[j];
        next[j] = false;
      }
    }
  }
}

EdgeMesh *EdgeMesh::read_preprocessed(const char *filename) {
  int numvertices, numfaces;
  bool read_faces = false, read_curve = false;
  char format[8];
  FILE *f = !strcmp(filename, "-") ? stdin : fopen(filename, "r");
  if (!f) return NULL;

  fread(format, sizeof(char), 8, f);
  if (strncmp(format, "FFFF", 4)) return NULL;
  if (!strncmp(format + 4, "FACE", 4)) read_faces = true;
  else if (!strncmp(format + 4, "CURV", 4)) read_curve = true;
  else if (!strncmp(format + 4, "VERT", 4));
  else assert(0);
  fread(&numvertices, sizeof(int), 1, f);
  if (read_faces) fread(&numfaces, sizeof(int), 1, f);

#if 0
  fprintf(stderr, "Reading %d vertices...\nReading %d faces...\n",
          numvertices, read_faces ? numfaces : 0);
#endif

  EdgeMesh *m = new EdgeMesh();
  m->vertices.resize(numvertices);
  m->normals.resize(numvertices);
  m->isedge.resize(numvertices);
  fread(&m->vertices[0], sizeof(float), 3 * numvertices, f);
  fread(&m->normals[0],  sizeof(float), 3 * numvertices, f);
  if (read_curve) {
    // fprintf(stderr, "Reading curvatures...\n");
    m->curv1.resize(numvertices);
    fread(&m->curv1[0], sizeof(float), numvertices, f);
    m->curv2 = m->curv1;
  }
  fread(m->isedge.isedge, sizeof(int), (numvertices + 31) / 32, f);

  if (read_faces) {
    m->faces.resize(numfaces);
    fread(&m->faces[0], sizeof(int), 3 * numfaces, f);
    m->need_adjacentfaces();
  }

  if (f != stdin) fclose(f);
  return m;
}

EdgeMesh *EdgeMesh::read(const char *filename) {
  EdgeMesh *mesh = new EdgeMesh();

  if (read_helper(filename, mesh))
    return mesh;

  delete mesh;
  return NULL;
}
