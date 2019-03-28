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

#include "mesh_sampler.h"

MeshSampler::MeshSampler(BasicTriMesh* _surfaceMesh, int _numPoints)  {
	surface_mesh = _surfaceMesh;
	num_points = _numPoints;
	num_iterations = 15;
	min_iterations = 10;
	energy_eps = 0.1;
	this->init();
}

MeshSampler::MeshSampler(BasicTriMesh* _surfaceMesh, int _numPoints, int _numIterations)  {
	surface_mesh = _surfaceMesh;
	num_points = _numPoints;
	num_iterations = _numIterations;
	min_iterations = 10;
	energy_eps = 0.1;
	this->init();
}

MeshSampler::MeshSampler(BasicTriMesh* _surfaceMesh, int _numPoints, int _numIterations, int _minIterations)  {
	surface_mesh = _surfaceMesh;
	num_points = _numPoints;
	num_iterations = _numIterations;
	min_iterations = _minIterations;
	energy_eps = 0.1;
	this->init();
}

MeshSampler::MeshSampler(BasicTriMesh* _surfaceMesh, int _numPoints, int _numIterations, int _minIterations, double _energyEps)  {
	surface_mesh = _surfaceMesh;
	num_points = _numPoints;
	num_iterations = _numIterations;
	min_iterations = _minIterations;
	energy_eps = _energyEps;
	this->init();
}

MeshSampler::~MeshSampler()  {
	delete [] sampled_pts;
}

void MeshSampler::sample(Vector3* points, int* tris)  {
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
			uniform_grid->add_point(GridPoint(sampled_pts[s].pt,s));
		}

		int periodicity = std::max(1, num_points/10);

		vector<SamplerPoint> new_samples;
		ComputationTimer timer("relaxation");
		timer.start();
		double total_energy = 0;
		// and perform repulsion
		for(int s = 0; s < num_points; s++)  {
			if((s % periodicity) == 0)
				cout << "point["<<s<<"]" << endl;
			SamplerPoint sampler_pt = sampled_pts[s];
			Vector3 pt = sampler_pt.pt;

			// gather neighborhood
			vector<SamplerPoint> neighborhood;
			this->gather_neighborhood(sampler_pt, uniform_grid, &neighborhood);
			if(neighborhood.size() == 1)
				continue;

			// compute energy & repulsed point
			total_energy += compute_energy(sampler_pt, &neighborhood);
			Vector3 force_vector = construct_force(sampler_pt, &neighborhood);
			Vector3 repulsed_pt = pt + force_vector;

			// project point onto nearest triangle
			SamplerPoint new_pt = this->projection(repulsed_pt, sampler_pt);

			// assign new point, and update grid
			new_samples.push_back(new_pt);
		}

		for(unsigned s = 0; s < new_samples.size(); s++)
			update_samples(new_samples[s], uniform_grid);

		timer.end();
		cout << timer.getComputation() << " : " << timer.getElapsedTime() << "s iteration: " << iteration << " energy: " << total_energy << endl;

		bool energy_achieved = false;
		if(iteration > 0)  {
			double energy_fraction = (prior_energy-total_energy)/total_energy;
			energy_achieved = iteration >= (min_iterations-1) && energy_fraction < energy_eps;
		}

		iteration++;
		sampling_satisfied = energy_achieved || (iteration == num_iterations);
		prior_energy = total_energy;
	}

	delete uniform_grid;

	for(int s = 0; s < num_points; s++)  {
		points[s] = sampled_pts[s].pt;
		int tri_ind = sampled_pts[s].tri;
		tris[s] = tri_ind;
	}
}

void MeshSampler::sample(Vector3* points, Vector3* normals)  {
	int* tris = new int[num_points];
	this->sample(points, tris);
	for(int s = 0; s < num_points; s++)
		normals[s] = sampler_tris[ tris[s] ].normal;
	delete [] tris;
}

void MeshSampler::init()  {
	surface_area = 0;
	pi_over_2 = 0.5*acos(-1);
	vector<TriIndex> all_tris;
	ll_bound = Vector3(1e10,1e10,1e10);
	ur_bound = Vector3(-1e10,-1e10,-1e10);

	for(int t = 0; t < surface_mesh->n_faces(); t++)  {
		BasicTriMesh::FaceVertexIter fv_it = surface_mesh->fv_iter(OpenMesh::FaceHandle(t));
		BasicTriMesh::VertexHandle v0 = fv_it.handle();
		BasicTriMesh::VertexHandle v1 = (++fv_it).handle();
		BasicTriMesh::VertexHandle v2 = (++fv_it).handle();

		BasicTriMesh::Point p0 = surface_mesh->point(v0), p1 = surface_mesh->point(v1), p2 = surface_mesh->point(v2);
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
	global_radius = 2.0*sqrt(surface_area / ((double)num_points));

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

	sampled_pts = new SamplerPoint[num_points];
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

		sampled_pts[s] = SamplerPoint(rand_pt, rand_tri_index.tri, s);
	}

	grid_res_x = (ur_bound.x-ll_bound.x)/global_radius;
	grid_res_y = (ur_bound.y-ll_bound.y)/global_radius;
	grid_res_z = (ur_bound.z-ll_bound.z)/global_radius;
	grid_res_x = grid_res_x > 300 ? 300 : grid_res_x;
	grid_res_y = grid_res_y > 300 ? 300 : grid_res_y;
	grid_res_z = grid_res_z > 300 ? 300 : grid_res_z;
	cout << "grid res: " << grid_res_x << " : " << grid_res_y << " : " << grid_res_z << endl;
}

int MeshSampler::tri_search(vector<TriIndex>* _tris, double _query, int _start, int _end)  {
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

void MeshSampler::gather_neighborhood(SamplerPoint _pt, UniformGrid* _grid, vector<SamplerPoint>* neighborhood)  {
	vector<GridPoint> grid_neighborhood;
	_grid->neighborhood_query(_pt.pt, global_radius, &grid_neighborhood);

	Vector3 pt = _pt.pt;
	Vector3 normal = sampler_tris[_pt.tri].normal;

	for(unsigned g = 0; g < grid_neighborhood.size(); g++)  {
		GridPoint grid_point = grid_neighborhood[g];
		SamplerPoint sampler_pt = sampled_pts[grid_point.ind];
		Vector3 next_normal = sampler_tris[sampler_pt.tri].normal;
		if(next_normal.dotProduct(normal) < 0)
			continue;
		neighborhood->push_back(sampler_pt);
	}
}

SamplerPoint MeshSampler::projection(Vector3 _pt, SamplerPoint _samplerPt)  {
	int tri_ind = _samplerPt.tri;
	SamplerTri tri = sampler_tris[tri_ind];

	Vector3 projected_pt;
	if(this->point_triangle_projection(_pt, tri, projected_pt))
		return SamplerPoint(projected_pt, tri_ind, _samplerPt.s);

	double min_dist = projected_pt.sqdDist(_pt);
	int best_tri = tri_ind;
	/*
	for(BasicTriMesh::FaceFaceIter ff_it = surface_mesh->ff_iter(BasicTriMesh::FaceHandle(tri_ind)); ff_it; ++ff_it)  {
		BasicTriMesh::FaceHandle f_handle = ff_it.handle();
		int neighbor_face = f_handle.idx();
		Vector3 face_projection;
		this->point_triangle_projection(_pt, sampler_tris[neighbor_face], face_projection);

		double sqd_dist = face_projection.sqdDist(_pt);
		if(sqd_dist < min_dist)  {
			min_dist = sqd_dist;
			best_tri = neighbor_face;
			projected_pt = face_projection;
		}
	}
	*/

	BasicTriMesh::FaceVertexIter fv_it = surface_mesh->fv_iter(OpenMesh::FaceHandle(tri_ind));
	BasicTriMesh::VertexHandle vh0 = fv_it.handle();
	BasicTriMesh::VertexHandle vh1 = (++fv_it).handle();
	BasicTriMesh::VertexHandle vh2 = (++fv_it).handle();
	BasicTriMesh::VertexHandle all_vhs[] = { vh0, vh1, vh2 };

	vector<int> face_ring;
	for(int v = 0; v < 3; v++)  {
		BasicTriMesh::VertexHandle next_vh = all_vhs[v];
		for(BasicTriMesh::VertexFaceIter vf_it = surface_mesh->vf_iter(next_vh); vf_it; ++vf_it)  {
			int next_face = vf_it.handle().idx();
			if(next_face == tri_ind)
				continue;
			bool face_found = false;
			for(unsigned f = 0; f < face_ring.size(); f++)  {
				if(next_face == face_ring[f])  {
					face_found = true;
					break;
				}
			}
			if(face_found)
				continue;
			face_ring.push_back(next_face);
		}
	}

	for(unsigned f = 0; f < face_ring.size(); f++)  {
		int neighbor_face = face_ring[f];
		Vector3 face_projection;
		this->point_triangle_projection(_pt, sampler_tris[neighbor_face], face_projection);
		double sqd_dist = face_projection.sqdDist(_pt);
		if(sqd_dist < min_dist)  {
			min_dist = sqd_dist;
			best_tri = neighbor_face;
			projected_pt = face_projection;
		}
	}

	return SamplerPoint(projected_pt, best_tri, _samplerPt.s);
}

bool MeshSampler::point_triangle_projection(Vector3 _pt, SamplerTri _tri, Vector3& newPt)  {
	Vector3 v1 = _tri.v1, v2 = _tri.v2, v3 = _tri.v3;
	Vector3 normal = _tri.normal;
	if(_pt == v1)  {
		newPt = v1;
		return true;
	}
	else if(_pt == v2)  {
		newPt = v2;
		return true;
	}
	else if(_pt == v3)  {
		newPt = v3;
		return true;
	}

	// --- first step: project point onto plane, and subdivide triangle based on projected point --- //
	Vector3 orthog_pt = _pt + normal*(v1-_pt).dotProduct(normal);
	Vector3 area_normal1 = (v2-orthog_pt).cross(v3-orthog_pt);
	Vector3 area_normal2 = (v3-orthog_pt).cross(v1-orthog_pt);
	Vector3 area_normal3 = (v1-orthog_pt).cross(v2-orthog_pt);
	bool is_area1_positive = area_normal1.dotProduct(normal) > 0;
	bool is_area2_positive = area_normal2.dotProduct(normal) > 0;
	bool is_area3_positive = area_normal3.dotProduct(normal) > 0;

	// is our projection in the triangle?
	if(is_area1_positive && is_area2_positive && is_area3_positive)  {
		newPt = orthog_pt;
		return true;
	}

	// otherwise, find the point on the border of the triangle closest to orthog_pt
	Vector3 tri_verts[3];
	tri_verts[0] = v1; tri_verts[1] = v2; tri_verts[2] = v3;

	// first try the edges ...
	bool is_contained = false;
	for(int i = 0; i < 3; i++)  {
		Vector3 q1 = tri_verts[i], q2 = tri_verts[(i+1)%3], q3 = tri_verts[(i+2)%3];
		Vector3 e1 = q2-q1, e2 = q1-q2;
		Vector3 edge_perp = e2.cross(normal);

		is_contained = signed_dist(orthog_pt, q1, edge_perp) > 0 && signed_dist(orthog_pt, q1, e1) > 0 && signed_dist(orthog_pt, q2, e2) > 0;
		if(is_contained)  {
			newPt = this->edge_projection(q1, q2, normal, orthog_pt);
			return false;
		}
	}

	// ... and last, find closest vertex
	int min_ind = 2;
	double sqd_dist1 = v1.sqdDist(orthog_pt), sqd_dist2 = v2.sqdDist(orthog_pt), sqd_dist3 = v3.sqdDist(orthog_pt);
	if(sqd_dist1 < sqd_dist2 && sqd_dist1 < sqd_dist3)
		min_ind = 0;
	else if(sqd_dist2 < sqd_dist3)
		min_ind = 1;

	newPt = tri_verts[min_ind];
	return false;
}

void MeshSampler::update_samples(SamplerPoint _newPt, UniformGrid* _grid)  {
	SamplerPoint pt = sampled_pts[_newPt.s];

	if(!_grid->remove_point(pt.pt))
		cerr << "error: should have removed point!" << endl;

	_grid->add_point(GridPoint(_newPt.pt, _newPt.s));
	sampled_pts[_newPt.s] = _newPt;
}

Vector3 MeshSampler::construct_force(SamplerPoint _point, vector<SamplerPoint>* _neighborhood)  {
	Vector3 total_force(0,0,0);
	Vector3 pt = _point.pt;
	Vector3 normal = sampler_tris[_point.tri].normal;

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

double MeshSampler::compute_energy(SamplerPoint _point, vector<SamplerPoint>* _neighborhood)  {
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
