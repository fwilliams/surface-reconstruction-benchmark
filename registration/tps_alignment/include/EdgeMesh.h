#ifndef EDGE_MESH_H
#define EDGE_MESH_H

/*
Copyright 2007 Benedict Brown
Princeton University

EdgeMesh.h
A subclass of TriMesh that caches edge boundaries.  This is necessary
because .pre file generally store boundary information but not face
information.

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
#include <string.h>
#include "TriMesh.h"

// this is an edge mesh with a cached edge structure
class EdgeMesh : public TriMesh {
 public:
  // We need a bitvector to mark whether each vertex is an edge, but
  // we'd like to stream this efficiently to and from disk.  That
  // doesn't seem to be possible with vector<bool>, so we do this
  // instead.
  struct EdgeList {
    int *isedge;
    int size, words;
    void free()  { clear(); if (isedge) delete[] isedge; isedge = NULL; size = words = 0; }
    void clear() { size = 0; }
    void resize(int newsize) {
      if (newsize <= (32 * words)) size = newsize;
      else {
        int newwords = (newsize + 31) / 32;
        assert(newwords > words);
        int *tmp = new int[newwords];
        if (words) memcpy(tmp, isedge, words * sizeof(int));
        delete[] isedge;
        isedge = tmp;
        size = newsize;
        words = newwords;
      }
    }

    bool operator()(int idx) const { return !!(isedge[idx / 32] & (1 << (idx % 32))); }
    void set(int idx, bool val) {
      if (val) isedge[idx / 32] |= (1 << (idx % 32));
      else isedge[idx / 32] &= ~(1 << (idx % 32));
    }

    EdgeList() : isedge(NULL), size(0), words(0) {};
    ~EdgeList() { free(); }
  } isedge;

  void need_edges(int bsize = 0);

  static EdgeMesh *read(const char *filename);
  static EdgeMesh *read_preprocessed(const char *filename);
};

#endif // EDGE_MESH_H
