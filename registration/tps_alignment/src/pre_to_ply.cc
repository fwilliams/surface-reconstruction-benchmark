/*
Copyright 2007 Benedict Brown
Princeton University

pre_to_ply.cc
Convert a .pre file back to a .ply

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

using namespace std;

FILE *outfile = stdout;

// stdin is list of input/output mesh pairs for full preprocess, input only for bbox
// argv[1] indicates full preprocess or bbox only
// argv[2] is output data file
// argv[3] if present is stable file list (only writes out these files/bboxes
int main(int argc, char *argv[])
{
  assert(argc >= 1);

  srand48(0);

  const char *fname = (argc > 1) ? argv[1] : "-";

  EdgeMesh *mesh = EdgeMesh::read_preprocessed(argv[1]);
  mesh->write(fname);
}
