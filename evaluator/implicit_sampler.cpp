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

#include "implicit_sampler.h"

ImplicitSampler::ImplicitSampler(ImplicitFunction* _surface, int _numPoints)  {
	surface = _surface;
	num_points = _numPoints;
	num_iterations = 50;
	min_iterations = 50;
	this->init();
}

ImplicitSampler::ImplicitSampler(ImplicitFunction* _surface, int _numPoints, int _numIterations)  {
	surface = _surface;
	num_points = _numPoints;
	num_iterations = _numIterations;
	min_iterations = 10;
	this->init();
}

ImplicitSampler::~ImplicitSampler()  {
	delete [] sampled_pts;
}

void ImplicitSampler::sample(Vector3* points, int* tris)  {
	UniformGrid* uniform_grid = new UniformGrid(grid_res_x, grid_res_y, grid_res_z, ll_bound, ur_bound);

	int iteration = 0;
	bool sampling_satisfied = iteration == num_iterations;

	double start_radius = global_radius, end_radius = 2.0*global_radius;

	double prior_energy = 0;
	while(!sampling_satisfied)  {
		// populate grid with new points
		for(int s = 0; s < num_points; s++)  {
			if(iteration > 0)
				break;
			uniform_grid->add_point(sampled_pts[s]);
		}

		int periodicity = num_points/10;

		vector<GridPoint> new_samples;
		ComputationTimer timer("relaxation");
		timer.start();
		double total_energy = 0;
		// and perform repulsion
		for(int s = 0; s < num_points; s++)  {
			if((s % periodicity) == 0)
				cout << "point["<<s<<"]" << endl;
			GridPoint sampler_pt = sampled_pts[s];
			Vector3 pt = sampler_pt.pt;

			// gather neighborhood
			vector<GridPoint> neighborhood;
			this->gather_neighborhood(sampler_pt, uniform_grid, &neighborhood);
			if(neighborhood.size() == 1)
				continue;

			// compute energy & repulsed point
			total_energy += compute_energy(sampler_pt, &neighborhood);
			Vector3 force_vector = construct_force(sampler_pt, &neighborhood);
			Vector3 repulsed_pt = pt + force_vector;

			// project point onto nearest triangle
			GridPoint new_pt = this->projection(repulsed_pt, sampler_pt);

			// assign new point, and update grid
			new_samples.push_back(new_pt);
		}

		for(unsigned s = 0; s < new_samples.size(); s++)
			update_samples(new_samples[s], uniform_grid);

		// total_energy /= num_points;
		timer.end();
		cout << timer.getComputation() << " : " << timer.getElapsedTime() << "s energy: " << total_energy << endl;

		bool energy_achieved = false;
		if(iteration > 0)  {
			double energy_fraction = (prior_energy-total_energy)/total_energy;
			energy_achieved = iteration >= (min_iterations-1) && energy_fraction < 0.1;
		}

		iteration++;
		sampling_satisfied = energy_achieved || (iteration == num_iterations);
		prior_energy = total_energy;
	}

	delete uniform_grid;

	for(int s = 0; s < num_points; s++)
		points[s] = sampled_pts[s].pt;
}

void ImplicitSampler::sample(Vector3* points, Vector3* normals)  {
	int* tris = new int[num_points];
	this->sample(points, tris);
	for(int s = 0; s < num_points; s++)
		normals[s] = surface->normal(points[s]);
	delete [] tris;
}

void ImplicitSampler::init()  {
	surface_area = 0;
	pi_over_2 = 0.5*acos(-1);
	vector<TriIndex> all_tris;
	BasicTriMesh* mc_mesh = surface->isosurface(64,256,0.01);
	ll_bound = Vector3(1e10,1e10,1e10);
	ur_bound = Vector3(-1e10,-1e10,-1e10);

	for(int t = 0; t < mc_mesh->n_faces(); t++)  {
		BasicTriMesh::FaceVertexIter fv_it = mc_mesh->fv_iter(OpenMesh::FaceHandle(t));
		BasicTriMesh::VertexHandle v0 = fv_it.handle();
		BasicTriMesh::VertexHandle v1 = (++fv_it).handle();
		BasicTriMesh::VertexHandle v2 = (++fv_it).handle();

		BasicTriMesh::Point p0 = mc_mesh->point(v0), p1 = mc_mesh->point(v1), p2 = mc_mesh->point(v2);
		BasicTriMesh::Point raw_normal = (p1-p0) % (p2-p0);
		double tri_area = 0.5*raw_normal.length();

		Vector3 vert0(p0[0],p0[1],p0[2]);
		Vector3 vert1(p1[0],p1[1],p1[2]);
		Vector3 vert2(p2[0],p2[1],p2[2]);
		sampler_tris.push_back(SamplerTri(vert0,vert1,vert2));
		ll_bound.expand_min_bound(vert0); ur_bound.expand_max_bound(vert0);
		ll_bound.expand_min_bound(vert1); ur_bound.expand_max_bound(vert1);
		ll_bound.expand_min_bound(vert2); ur_bound.expand_max_bound(vert2);

		if(p0 == p1 || p0 == p2 || p1 == p2 || tri_area < 1e-9)
			continue;
		all_tris.push_back(TriIndex(tri_area, t));
		surface_area += tri_area;
	}
	std::sort(all_tris.begin(),all_tris.end());
	global_radius = 2.5*sqrt(surface_area / ((double)num_points));

	vector<TriIndex> cum_tris;
	for(unsigned t = 0; t < all_tris.size(); t++)  {
		TriIndex next_tri = all_tris[t];
		if(t == 0)
			cum_tris.push_back(next_tri);
		else  {
			TriIndex prev_tri = cum_tris[t-1];
			double cum_area = next_tri.area + prev_tri.area;
			cum_tris.push_back(TriIndex(cum_area, next_tri.tri));
		}
	}

	time_t cur_time;
	time(&cur_time);
	srand(cur_time);
	for(int i = 0; i < 20; i++)
		rand();

	sampled_pts = new GridPoint[num_points];
	double total_area = cum_tris[cum_tris.size()-1].area;
	for(int s = 0; s < num_points; s++)  {
		double rand_unit = (double)rand() / (double)RAND_MAX;
		double rand_area = rand_unit*total_area;
		int rand_ind = this->tri_search(&cum_tris, rand_area, 0, cum_tris.size()-1);
		TriIndex rand_tri_index = cum_tris[rand_ind];

		SamplerTri rand_tri = sampler_tris[rand_tri_index.tri];
		double r1 = (double)rand() / (double)RAND_MAX;
		double r2 = (double)rand() / (double)RAND_MAX;
		double r1_sqrt = sqrt(r1);
		Vector3 rand_pt = rand_tri.v1*(1.0-r1_sqrt) + rand_tri.v2*(r1_sqrt*(1.0-r2)) + rand_tri.v3*(r1_sqrt*r2);

		sampled_pts[s] = GridPoint(rand_pt, s);
	}

	grid_res_x = (ur_bound.x-ll_bound.x)/global_radius;
	grid_res_y = (ur_bound.y-ll_bound.y)/global_radius;
	grid_res_z = (ur_bound.z-ll_bound.z)/global_radius;
	grid_res_x = grid_res_x > 300 ? 300 : grid_res_x;
	grid_res_y = grid_res_y > 300 ? 300 : grid_res_y;
	grid_res_z = grid_res_z > 300 ? 300 : grid_res_z;
	cout << "grid res: " << grid_res_x << " : " << grid_res_y << " : " << grid_res_z << endl;

	delete mc_mesh;
}

int ImplicitSampler::tri_search(vector<TriIndex>* _tris, double _query, int _start, int _end)  {
	int start = _start, end = _end;
	double query = _query;

	while(start <= end)  {
		int mid = (start+end)/2;
		double mid_value = _tris->at(mid).area;
		double prev_mid_value = mid == 0 ? 0 : _tris->at(mid-1).area;

		if(query >= prev_mid_value && query <= mid_value)  {
			//cout << "got you! " << query << " : " << prev_mid_value << " : " << mid_value << endl;
			return mid;
		}
		else if(query > mid_value)
			start = mid+1;
		else if(query < prev_mid_value)
			end = mid-1;
		else  {
			cerr << "error in tri search!" << endl;
			break;
		}
	}

	cerr << "unable to find " << query << " in anything!" << endl;
	return 0;
}

void ImplicitSampler::gather_neighborhood(GridPoint _pt, UniformGrid* _grid, vector<GridPoint>* neighborhood)  {
	vector<GridPoint> grid_neighborhood;
	_grid->neighborhood_query(_pt.pt, global_radius, &grid_neighborhood);

	Vector3 pt = _pt.pt;
	// TODO: incorporate normals!
	//Vector3 normal;

	for(unsigned g = 0; g < grid_neighborhood.size(); g++)  {
		GridPoint grid_point = grid_neighborhood[g];
		GridPoint sampler_pt = sampled_pts[grid_point.ind];
		/*
		Vector3 next_normal = sampler_tris[sampler_pt.tri].normal;
		if(next_normal.dotProduct(normal) < 0)
			continue;
			*/
		neighborhood->push_back(sampler_pt);
	}
}

GridPoint ImplicitSampler::projection(Vector3 _pt, GridPoint _samplerPt)  {
	int projection_iterations = 10;
	double lambda = 0.75;
	double prev_dist = 1e10;
	Vector3 pt = _pt;

	for(int j = 0; j < projection_iterations; j++)  {
		double pt_function = surface->function(pt);
		Vector3 pt_gradient = surface->gradient(pt);
		Vector3 pt_normal = pt_gradient;
		pt_normal.normalize();
		pt_normal.scale(-(pt_function*lambda)/pt_gradient.length());
		pt = pt + pt_normal;

		double next_dist = pt_function / pt_gradient.length();
		double next_abs = fabs(next_dist);
		if(next_abs < 1e-3)
			break;
		if(next_abs > prev_dist)
			cerr << "error in projection: " << next_abs << " : " << prev_dist << endl;
		prev_dist = next_abs;
	}
	return GridPoint(pt, _samplerPt.ind);
}

void ImplicitSampler::update_samples(GridPoint _newPt, UniformGrid* _grid)  {
	GridPoint pt = sampled_pts[_newPt.ind];

	if(!_grid->remove_point(pt.pt))
		cerr << "error: should have removed point!" << endl;

	_grid->add_point(GridPoint(_newPt.pt, _newPt.ind));
	sampled_pts[_newPt.ind] = _newPt;
}

Vector3 ImplicitSampler::construct_force(GridPoint _point, vector<GridPoint>* _neighborhood)  {
	Vector3 total_force(0,0,0);
	Vector3 pt = _point.pt;
	// TODO: precompute normals!
	Vector3 normal = surface->normal(pt);

	double total_weight = 0;
	for(int n = 0; n < _neighborhood->size(); n++)  {
		Vector3 next_pt = _neighborhood->at(n).pt;
		if(next_pt == pt)
			continue;

		Vector3 tangent_pt = next_pt + normal*((pt-next_pt).dotProduct(normal));
		double dist = pt.length(next_pt);
		double sin_weight_sqd = 1.0 / sin(pi_over_2*(dist/global_radius));
		sin_weight_sqd *= sin_weight_sqd;
		double weight = dist > global_radius ? 0 : (-pi_over_2*(1.0-sin_weight_sqd));
		total_weight += weight;

		Vector3 repulsion_vec = pt-tangent_pt;
		total_force = total_force + repulsion_vec*weight;
	}
	total_force.scale(1.0/total_weight);

	return total_force;
}

double ImplicitSampler::compute_energy(GridPoint _point, vector<GridPoint>* _neighborhood)  {
	Vector3 pt = _point.pt;
	double total_energy = 0;
	for(int n = 0; n < _neighborhood->size(); n++)  {
		Vector3 next_pt = _neighborhood->at(n).pt;
		if(next_pt == pt)
			continue;

		double dist = pt.length(next_pt);
		if(dist > global_radius)
			continue;

		double angle = pi_over_2*(dist/global_radius);
		double cot_weight = cos(angle) / sin(angle);
		total_energy += cot_weight + (dist/global_radius)*pi_over_2 - pi_over_2;
	}
	return total_energy;
}
