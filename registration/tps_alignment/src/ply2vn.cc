/*
Copyright 2007 Benedict Brown
Princeton University

ply2vn.cc
Take a mesh and spit out the list of vertices and normals

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

#include "TriMesh.h"

using namespace std;

int main(int argc, char *argv[])
{
  assert(argc >= 1);

  srand48(0);

  TriMesh::set_verbose(0);

  for (int f = 1; f < argc; f++) {
    TriMesh *mesh = TriMesh::read(argv[f]);
    mesh->need_normals();
    for (int i = 0; i < (int) mesh->vertices.size(); i++) {
#if 0
      // ascii
      printf("%f %f %f %f %f %f\n", mesh->vertices[i][0],
             mesh->vertices[i][1],  mesh->vertices[i][2],
             mesh->normals[i][0], mesh->normals[i][1], mesh->normals[i][2]);
#else
      // binary
      fwrite(&mesh->vertices[i][0], sizeof(float), 3, stdout);
      fwrite(&mesh->normals[i][0],  sizeof(float), 3, stdout);
#endif
    }
    delete mesh;
  }
}
