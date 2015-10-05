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

#ifndef IMPLICITSAMPLER_H
#define IMPLICITSAMPLER_H

#include "UniformGrid.h"
#include "mesh_sampler.h"
#include "../modeling/ImplicitFunction.h"
#include "mesh_types.h"
#include "ComputationTimer.h"
#include <algorithm>
using namespace std;

class ImplicitSampler  {
	public:
		ImplicitSampler(ImplicitFunction* _surface, int _numPoints);
		ImplicitSampler(ImplicitFunction* _surface, int _numPoints, int _numIterations);
		~ImplicitSampler();

		int sampling_size()  { return num_points; }
		void sample(Vector3* points, int* tris);
		void sample(Vector3* points, Vector3* normals);

	private:
		ImplicitFunction* surface;
		int num_points, num_iterations, min_iterations;
		double surface_area;

		Vector3 ll_bound, ur_bound;
		vector<SamplerTri> sampler_tris;

		int grid_res_x, grid_res_y, grid_res_z;
		double global_radius, pi_over_2;

		GridPoint* sampled_pts;

		void init();

		int tri_search(vector<TriIndex>* _tris, double _query, int _start, int _end);

		void gather_neighborhood(GridPoint _pt, UniformGrid* _grid, vector<GridPoint>* neighborhood);
		Vector3 construct_force(GridPoint _point, vector<GridPoint>* _neighborhood);
		double compute_energy(GridPoint _point, vector<GridPoint>* _neighborhood);

		GridPoint projection(Vector3 _pt, GridPoint _samplerPt);

		void update_samples(GridPoint _newPt, UniformGrid* _grid);
};

#endif
