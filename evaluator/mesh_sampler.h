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

#ifndef MESHSAMPLER_H
#define MESHSAMPLER_H

#include "UniformGrid.h"
#include "../modeling/Vector3.h"
#include "mesh_types.h"
#include "ComputationTimer.h"
#include <algorithm>
using namespace std;

struct TriIndex  {
	TriIndex() {}
	TriIndex(double _area, int _tri) : area(_area) , tri(_tri) {}
	inline bool operator<(const TriIndex& _otherTri) const {
		return area < _otherTri.area;
	}
	double area;
	int tri;
};

struct SamplerTri  {
	SamplerTri()  {}
	SamplerTri(Vector3 _v1, Vector3 _v2, Vector3 _v3)  {
		v1 = _v1;
		v2 = _v2;
		v3 = _v3;
		normal = (v2-v1).cross(v3-v1);
		area = 0.5*normal.length();
		normal.normalize();
	}
	Vector3 v1, v2, v3, normal;
	double area;
};

struct SamplerPoint  {
	SamplerPoint()  {}
	SamplerPoint(Vector3 _pt, int _tri, int _s) : pt(_pt) , tri(_tri) , s(_s)  {}
	Vector3 pt;
	int tri, s;
};

class MeshSampler  {
	public:
		MeshSampler(BasicTriMesh* _surfaceMesh, int _numPoints);
		MeshSampler(BasicTriMesh* _surfaceMesh, int _numPoints, int _numIterations);
		MeshSampler(BasicTriMesh* _surfaceMesh, int _numPoints, int _numIterations, int _minIterations);
		MeshSampler(BasicTriMesh* _surfaceMesh, int _numPoints, int _numIterations, int _minIterations, double _energyEps);
		~MeshSampler();

		int sampling_size()  { return num_points; }
		void sample(Vector3* points, int* tris);
		void sample(Vector3* points, Vector3* normals);

	private:
		BasicTriMesh* surface_mesh;
		int num_points, num_iterations, min_iterations;
		double energy_eps;
		double surface_area;

		Vector3 ll_bound, ur_bound;
		vector<SamplerTri> sampler_tris;

		int grid_res_x, grid_res_y, grid_res_z;
		double global_radius, pi_over_2;

		SamplerPoint* sampled_pts;

		void init();

		int tri_search(vector<TriIndex>* _tris, double _query, int _start, int _end);

		void gather_neighborhood(SamplerPoint _pt, UniformGrid* _grid, vector<SamplerPoint>* neighborhood);
		Vector3 construct_force(SamplerPoint _point, vector<SamplerPoint>* _neighborhood);
		double compute_energy(SamplerPoint _point, vector<SamplerPoint>* _neighborhood);

		SamplerPoint projection(Vector3 _pt, SamplerPoint _samplerPt);
		bool point_triangle_projection(Vector3 _pt, SamplerTri _tri, Vector3& newPt);

		void update_samples(SamplerPoint _newPt, UniformGrid* _grid);

		double signed_dist(Vector3 _p, Vector3 _x, Vector3 _normal)  {
			return ((_p-_x).dotProduct(_normal)) / _normal.dotProduct(_normal);
		}
		Vector3 edge_projection(Vector3 _v1, Vector3 _v2, Vector3 _normal, Vector3 _point)  {
			Vector3 edge_normal = (_v2-_v1).cross(_normal);
			edge_normal.normalize();
			return _point + edge_normal*(_v1-_point).dotProduct(edge_normal);
		}

		void barycentric_coordinate(SamplerPoint _pt, double& b1, double& b2, double& b3)  {
			SamplerTri tri = sampler_tris[_pt.tri];
			Vector3 pt = _pt.pt;
			double area1 = 0.5*((pt-tri.v2).cross(pt-tri.v3)).length();
			double area2 = 0.5*((pt-tri.v3).cross(pt-tri.v1)).length();
			double area3 = 0.5*((pt-tri.v1).cross(pt-tri.v2)).length();
			double inv_area = 1.0 / tri.area;
			b1 = area1*inv_area;
			b2 = area2*inv_area;
			b3 = area3*inv_area;
		}
};

#endif
