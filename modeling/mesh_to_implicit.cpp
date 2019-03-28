/*
Copyright (c) 2010, Matt Berger
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list
  of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list
  of conditions and the following disclaimer in the documentation and/or other
  materials provided with the distribution.
- Neither the name of the University of Utah nor the names of its contributors may be used
  to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <OpenMesh/Core/IO/MeshIO.hh>
#include "shape_loader.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
using namespace std;

#define MAX_DIM_MM 75

int main(int argc, char** argv)  {
	if(argc != 6)  {
		cerr << "usage: " << argv[0] << " mesh_in implicit_out min_samples fit_eps covering" << endl;
		return 0;
	}

	int arg_num=1;
	string mesh_in = argv[arg_num++];
	string implicit_out = argv[arg_num++];

	BasicTriMesh surface_mesh;
	OpenMesh::IO::Options opt;
	OpenMesh::IO::read_mesh(surface_mesh, mesh_in, opt);

	// resize mesh to largest dimension being 75mm
	BasicTriMesh::Point ll, ur;
	for(int i = 0; i < 3; i++)  {
		ll[i] = 1e10;
		ur[i] = -1e10;
	}

	for(int p = 0; p < surface_mesh.n_vertices(); p++)  {
		BasicTriMesh::Point pt = surface_mesh.point(BasicTriMesh::VertexHandle(p));
		for(int i = 0; i < 3; i++)  {
			ll[i] = pt[i] < ll[i] ? pt[i] : ll[i];
			ur[i] = pt[i] > ur[i] ? pt[i] : ur[i];
		}
	}

	BasicTriMesh::Point bounds = ur-ll;
	double max_bound = -1e10;
	int max_dim = -1;
	for(int i = 0; i < 3; i++)  {
		if(bounds[i] > max_bound)  {
			max_dim = i;
			max_bound = bounds[i];
		}
	}

	BasicTriMesh::Point rescaled_bounds;
	for(int i = 0; i < 3; i++)
		rescaled_bounds[i] = (MAX_DIM_MM/max_bound)*bounds[i];

	BasicTriMesh::Point center = (ur+ll)*0.5;
	for(int p = 0; p < surface_mesh.n_vertices(); p++)  {
		BasicTriMesh::Point pt = surface_mesh.point(BasicTriMesh::VertexHandle(p));
		BasicTriMesh::Point normalized_pt = pt-center;
		for(int i = 0; i < 3; i++)
			normalized_pt[i] = normalized_pt[i] / (0.5*bounds[i]);
		BasicTriMesh::Point new_pt;
		for(int i = 0; i < 3; i++)
			new_pt[i] = normalized_pt[i] * (0.5*rescaled_bounds[i]);
		surface_mesh.set_point(BasicTriMesh::VertexHandle(p), new_pt);
	}

	BasicTriMesh::Point new_ll, new_ur;
	for(int i = 0; i < 3; i++)  {
		new_ll[i] = 1e10;
		new_ur[i] = -1e10;
	}

	for(int p = 0; p < surface_mesh.n_vertices(); p++)  {
		BasicTriMesh::Point pt = surface_mesh.point(BasicTriMesh::VertexHandle(p));
		for(int i = 0; i < 3; i++)  {
			new_ll[i] = pt[i] < new_ll[i] ? pt[i] : new_ll[i];
			new_ur[i] = pt[i] > new_ur[i] ? pt[i] : new_ur[i];
		}
	}

	cout << "old bounds: " << ll << " : " << ur << endl;
	cout << "new bounds: " << new_ll << " : " << new_ur << endl;

	int min_samples = atoi(argv[arg_num++]);
	double fit_eps = atof(argv[arg_num++]);
	double covering = atof(argv[arg_num++]);

	cout << "mesh to mpu..." << endl;
	MPUSoup* mpu_surface = new MPUSoup(&surface_mesh, min_samples, fit_eps, covering);
	mpu_surface->write_to_file(implicit_out);
	cout << "... done with mesh to mpu" << endl;
}
