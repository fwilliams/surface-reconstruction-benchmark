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

#include "mesh_types.h"
#include "kd_tree.h"

int main(int argc, char **argv) {
	if(argc != 3)  {
		std::cerr << "usage: " << argv[0] << " input_mesh output_mesh" << std::endl;
		return 1;
	}

	int arg_num = 1;
	string input_filename = argv[arg_num++];
	string output_filename = argv[arg_num++];

	BasicTriMesh input_mesh;
   OpenMesh::IO::Options opt;
	if(!OpenMesh::IO::read_mesh(input_mesh, input_filename, opt))  {
		std::cerr << "unable to read mesh from file " << argv[arg_num] << std::endl;
		return 1;
	}

	// setup kd-tree of vertices
	SimpleKdTree kd_tree(input_mesh.n_vertices());
	Vector3 min_bound(1e10,1e10,1e10), max_bound(-1e10,-1e10,-1e10);
	for(unsigned i = 0; i < input_mesh.n_vertices(); i++)  {
		BasicTriMesh::Point mesh_pt = input_mesh.point(OpenMesh::VertexHandle(i));
		Vector3 next_pt(mesh_pt[0], mesh_pt[1], mesh_pt[2]);
		min_bound.expand_min_bound(next_pt);
		max_bound.expand_max_bound(next_pt);
		kd_tree.add_point(i, next_pt);
	}
	kd_tree.finalize_tree();
	double tree_eps = 1e-12*(max_bound-min_bound).length();

	Vector3 answer_pts[25];
	int answer_ids[25];
	int vertex_remapped = 0;
	int* vertex_mapping = new int[input_mesh.n_vertices()];
	bool* is_unique_vertex = new bool[input_mesh.n_vertices()];
	for(int v = 0; v < input_mesh.n_vertices(); v++)  {
		vertex_mapping[v] = -1;
		is_unique_vertex[v] = false;
	}

	for(int v = 0; v < input_mesh.n_vertices(); v++)  {
		BasicTriMesh::Point mesh_pt = input_mesh.point(OpenMesh::VertexHandle(v));
		Vector3 next_pt(mesh_pt[0], mesh_pt[1], mesh_pt[2]);
		int num_found;
		kd_tree.closest_point_query(next_pt, 25, tree_eps, answer_pts, answer_ids, num_found);

		if(num_found == 0)
			cout << "NOTHING HERE?" << endl;
		if(num_found == 1)  {
			vertex_mapping[v] = vertex_remapped++;
			is_unique_vertex[v] = true;
			cout << "vertex is isolated?" << endl;
		}
		else  {
			// go through ids, and find the one which has been mapped
			bool found_vertex = false;
			for(int i = 0; i < num_found; i++)  {
				int found_vert = answer_ids[i];
				if(vertex_mapping[found_vert] != -1)  {
					found_vertex = true;
					vertex_mapping[v] = vertex_mapping[found_vert];
					break;
				}
			}
			if(!found_vertex)  {
				vertex_mapping[v] = vertex_remapped++;
				is_unique_vertex[v] = true;
			}
		}
	}

	BasicTriMesh output_mesh;
	for(int v = 0; v < input_mesh.n_vertices(); v++)  {
		if(!is_unique_vertex[v])
			continue;
		output_mesh.add_vertex(input_mesh.point(OpenMesh::VertexHandle(v)));
	}

	int num_degenerate_edges = 0;
	for(int f = 0; f < input_mesh.n_faces(); f++)  {
		BasicTriMesh::FaceVertexIter fv_it = input_mesh.fv_iter(OpenMesh::FaceHandle(f));
		int v0 = fv_it.handle().idx();
		int v1 = (++fv_it).handle().idx();
		int v2 = (++fv_it).handle().idx();
		int proper_v0 = vertex_mapping[v0], proper_v1 = vertex_mapping[v1], proper_v2 = vertex_mapping[v2];

		if(proper_v0 == proper_v1 || proper_v0 == proper_v2 || proper_v1 == proper_v2)  {
			num_degenerate_edges++;
			continue;
		}
		OpenMesh::FaceHandle new_face =
			output_mesh.add_face(OpenMesh::VertexHandle(proper_v0), OpenMesh::VertexHandle(proper_v1), OpenMesh::VertexHandle(proper_v2));
	}
	cout << "num degenerate edges: " << num_degenerate_edges << endl;

	delete [] vertex_mapping;
	delete [] is_unique_vertex;

	OpenMesh::IO::write_mesh(output_mesh, output_filename, opt);
}
