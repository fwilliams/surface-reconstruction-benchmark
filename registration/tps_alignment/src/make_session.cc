/*
Copyright 2007 Benedict Brown
Princeton University

Make a .session file for scanalyze for a list of meshes

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
#include "TriMesh.h"

int main(int argc, char**argv) {
  for (int i = 1; i < argc; i++) {
    TriMesh *mesh = TriMesh::read(argv[i]);
    mesh->need_bbox();
    printf("%s: [%f %f %f] * [%f %f %f]\n", argv[i],
           mesh->bbox.min[0], mesh->bbox.min[1], mesh->bbox.min[2],
           mesh->bbox.max[0], mesh->bbox.max[1], mesh->bbox.max[2]);
    delete mesh;
  }

  return 0;
}
