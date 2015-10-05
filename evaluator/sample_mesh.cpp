/*
Copyright (c) 2011, Joshua A. Levine
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
- Neither the name of the University of Utah nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "mesh_sampler.h"
#include "../sampler/OrientedPointCloud.h"

#include <fstream>
#include <sstream>
#include <cstdio>

using namespace std;

int main(int argc, char **argv) {
	if(argc < 4)  {
		std::cerr << "usage: " << argv[0] << " mesh num_pts output_file [num_iters] [min_iters] [energy_eps]" << std::endl;
		return 1;
	}

	int arg_num = 1;

	// --- input --- //
	BasicTriMesh mesh;
   OpenMesh::IO::Options opt;
	if(!OpenMesh::IO::read_mesh(mesh, argv[arg_num++], opt))  {
		std::cerr << "unable to read mesh from file " << argv[arg_num] << std::endl;
		return 1;
	}

   int num_samples = atoi(argv[arg_num++]);
   string pc_filename = argv[arg_num++];

	int num_iters = 15, min_iters = 10;
	double energy_eps = 0.1;

	if(argc > 4)
		num_iters = atoi(argv[arg_num++]);
   if(num_iters < 15)
      num_iters = 15;

	if(argc > 5)
		min_iters = atoi(argv[arg_num++]);
   if(min_iters < 10)
      min_iters = 10;

	if(argc > 6)
		energy_eps = atof(argv[arg_num++]);
	if(energy_eps > 0.1)
		energy_eps = 0.1;

   Vector3* points = new Vector3[num_samples];
   Vector3* normals = new Vector3[num_samples];

   MeshSampler ms(&mesh, num_samples, num_iters, min_iters, energy_eps);

   ms.sample(points, normals);

	OrientedPointCloud pc;

   for (int i=0; i<num_samples; i++) {
      pc.addOrientedPoint(points[i], normals[i]);
   }

   pc.write_to_file(pc_filename);


   delete[] points;
   delete[] normals;

   return 0;

}
