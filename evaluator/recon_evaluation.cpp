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

#include "recon_evaluation.h"

ReconEvaluation::ReconEvaluation(ImplicitFunction* _shape, OrientedPointCloud* _pointShape, BasicTriMesh* _surfaceMesh)  {
	shape = _shape;
	point_shape = _pointShape;
	input_mesh = _surfaceMesh;
	determine_inversion = true;
	to_invert = false;

	// very first thing: deal with meta-data
	this->pre_process();

	n_shape_pts = point_shape->size();
	n_mesh_vertices = surface_mesh->n_vertices();
	n_mesh_tris = surface_mesh->n_faces();
	n_mesh_points = n_mesh_tris >= MIN_MESH_SAMPLES ? n_mesh_tris : MIN_MESH_SAMPLES;
	n_mesh_points = n_mesh_points <= MAX_MESH_SAMPLES ? n_mesh_points : MAX_MESH_SAMPLES;

	kd_tree = new KdOpenMesh(surface_mesh);

	sampled_shape = new Vector3[n_shape_pts];
	sampled_mesh = new Vector3[n_mesh_points];

	shape_tree = new ImplicitKdTree(n_shape_pts, shape);
	mesh_tree = new MeshKdTree(n_mesh_points, surface_mesh);

	ll_bound.x = ll_bound.y = ll_bound.z = 1e10;
	ur_bound.x = ur_bound.y = ur_bound.z = -1e10;
	for(int p = 0; p < n_shape_pts; p++)  {
		sampled_shape[p] = point_shape->getPoint(p);
		ll_bound.x = sampled_shape[p].x < ll_bound.x ? sampled_shape[p].x : ll_bound.x;
		ll_bound.y = sampled_shape[p].y < ll_bound.y ? sampled_shape[p].y : ll_bound.y;
		ll_bound.z = sampled_shape[p].z < ll_bound.z ? sampled_shape[p].z : ll_bound.z;
		ur_bound.x = sampled_shape[p].x > ur_bound.x ? sampled_shape[p].x : ur_bound.x;
		ur_bound.y = sampled_shape[p].y > ur_bound.y ? sampled_shape[p].y : ur_bound.y;
		ur_bound.z = sampled_shape[p].z > ur_bound.z ? sampled_shape[p].z : ur_bound.z;

		shape_tree->add_point(p, sampled_shape[p]);
	}
	shape_tree->finalize_tree();

	int* sampled_tris = new int[n_mesh_points];
	MeshSampler mesh_sampler(surface_mesh, n_mesh_points, 15);
	mesh_sampler.sample(sampled_mesh, sampled_tris);
	for(int i = 0; i < n_mesh_points; i++)
		mesh_tree->add_tri_point(i, sampled_tris[i], sampled_mesh[i]);
	delete [] sampled_tris;

	mesh_tree->finalize_tree();

	bool use_occupancy_grid = false;
	if(use_occupancy_grid)  {
		Vector3 ll_grid, ur_grid;
		shape->bounds(ll_grid, ur_grid);
		UniformOccupancyGrid* occupancy_grid = new UniformOccupancyGrid(ll_grid, ur_grid, 100);
		for(int p = 0; p < n_shape_pts; p++)
			occupancy_grid->addPoint(sampled_shape[p]);
		occupancy_grid->post_process();
		root_finder = new ImplicitNewton(shape, occupancy_grid, 1.0, 1e-7);
	}
	else
		root_finder = new ImplicitNewton(shape, 1.0, 1e-7);

	shape_to_mesh = new ShortestDistanceMap();
	mesh_to_shape = new ShortestDistanceMap();
}

ReconEvaluation::~ReconEvaluation()  {
	delete surface_mesh;
	delete sampled_shape;
	delete sampled_mesh;
	delete kd_tree;
	delete root_finder;

	delete mesh_to_shape;
	delete shape_to_mesh;

	delete shape_tree;
	delete mesh_tree;
	annClose();
}

void ReconEvaluation::pre_process()  {
	// first thing: extract connected components
	int* vertex_ccs = new int[input_mesh->n_vertices()];
	for(int v = 0; v < input_mesh->n_vertices(); v++)
		vertex_ccs[v] = 0;
	int* face_ccs = new int[input_mesh->n_faces()];
	for(int f = 0; f < input_mesh->n_faces(); f++)
		face_ccs[f] = 0;

	int current_face_id = 0, current_vertex_id = 0;
	for(int f = 0; f < input_mesh->n_faces(); f++)  {
		if(face_ccs[f] != 0)
			continue;

		current_face_id++;
		current_vertex_id++;
		int v0, v1, v2;
		// setup connected component list for first face & vertices
		this->getVertices(input_mesh, OpenMesh::FaceHandle(f), v0, v1, v2);
		face_ccs[f] = current_face_id;

		vector<int> face_queue;
		face_queue.push_back(f);

		// breadth first search over neighboring faces
		while(face_queue.size() > 0)  {
			int next_face = face_queue[0];
			face_queue.erase(face_queue.begin());
			OpenMesh::FaceHandle f_handle(next_face);

			int num_neighbors = 0;
			for(FaceFaceIter ff_it = input_mesh->ff_iter(f_handle); ff_it; ++ff_it)  {
				num_neighbors++;
				OpenMesh::FaceHandle neighbor_face = ff_it.handle();
				if(!input_mesh->is_valid_handle(neighbor_face) || face_ccs[neighbor_face.idx()] != 0)
					continue;
				// update queue with new face
				face_queue.push_back(neighbor_face.idx());

				// update connected component list for faces & vertices
				this->getVertices(input_mesh, neighbor_face, v0, v1, v2);
				if(v0 < 0 || v0 >= input_mesh->n_vertices())
					cerr << "what is this? " << v0 << endl;
				if(v1 < 0 || v1 >= input_mesh->n_vertices())
					cerr << "what is this? " << v1 << endl;
				if(v2 < 0 || v2 >= input_mesh->n_vertices())
					cerr << "what is this? " << v2 << endl;
				face_ccs[neighbor_face.idx()] = current_face_id;
			}
		}
	}
	num_connected_components = current_face_id;

	// now, find the main connected component, the one with the most faces: this is our surface
	int* cc_sizes = new int[current_face_id+1];
	for(int i = 0; i <= current_face_id; i++)
		cc_sizes[i] = 0;
	for(int f = 0; f < input_mesh->n_faces(); f++)
		cc_sizes[face_ccs[f]]++;
	int main_cc = 0, max_faces = -1;
	for(int i = 1; i <= current_face_id; i++)  {
		if(cc_sizes[i] > max_faces)  {
			max_faces = cc_sizes[i];
			main_cc = i;
		}
	}

	// go through all of the faces, and mark the vertices belonging to them as the connected component
	for(int f = 0; f < input_mesh->n_faces(); f++)  {
		if(face_ccs[f] == main_cc)  {
			int v0, v1, v2;
			this->getVertices(input_mesh, OpenMesh::FaceHandle(f), v0, v1, v2);
			vertex_ccs[v0] = main_cc;
			vertex_ccs[v1] = main_cc;
			vertex_ccs[v2] = main_cc;
		}
	}

	// extract vertices/faces of main connected component: our surface mesh
	is_manifold = true;
	surface_mesh = new BasicTriMesh();
	int* vertex_mapping = new int[input_mesh->n_vertices()];
	int total_vertices = 0;
	for(int v = 0; v < input_mesh->n_vertices(); v++)  {
		if(vertex_ccs[v] == main_cc)  {
			surface_mesh->add_vertex(input_mesh->point(OpenMesh::VertexHandle(v)));
			vertex_mapping[v] = total_vertices++;
		}
	}

	OpenMesh::FaceHandle new_face;
	for(int f = 0; f < input_mesh->n_faces(); f++)  {
		if(face_ccs[f] == main_cc)  {
			int v0, v1, v2;
			this->getVertices(input_mesh, OpenMesh::FaceHandle(f), v0, v1, v2);
			int proper_v0 = vertex_mapping[v0], proper_v1 = vertex_mapping[v1], proper_v2 = vertex_mapping[v2];
			new_face = surface_mesh->add_face(OpenMesh::VertexHandle(proper_v0), OpenMesh::VertexHandle(proper_v1), OpenMesh::VertexHandle(proper_v2));
			if(!surface_mesh->is_valid_handle(new_face))
				is_manifold = false;
		}
	}
	delete [] vertex_mapping;

	cout << "num connected components: " << num_connected_components << " is manifold? " << is_manifold << endl;

	// if it is manifold, then get genus, number of boundary components, & length of boundary components
	if(is_manifold)  {
		int* he_ids = new int[input_mesh->n_halfedges()];
		for(int n = 0; n < input_mesh->n_halfedges(); n++)
			he_ids[n] = 0;
		int boundary_id = 0;
		boundary_length = 0;

		for(int e = 0; e < input_mesh->n_halfedges(); e++)  {
			OpenMesh::HalfedgeHandle he(e);
			if(!surface_mesh->is_valid_handle(he))
				continue;
			OpenMesh::FaceHandle he_face = surface_mesh->face_handle(he);
			if(surface_mesh->is_valid_handle(he_face))
				continue;

			// found a boundary: if it exists, then continue on
			if(he_ids[he.idx()] != 0)
				continue;

			// otherwise, trace halfedges to construct boundary, until we hit the beginning
			boundary_id++;
			OpenMesh::HalfedgeHandle start_he = he;
			bool at_beginning = false;
			int num_traced = 0;
			while(!at_beginning)  {
				num_traced++;

				he_ids[he.idx()] = boundary_id;
				OpenMesh::VertexHandle vert_from = surface_mesh->from_vertex_handle(he), vert_to = surface_mesh->to_vertex_handle(he);
				boundary_length += (surface_mesh->point(vert_to)-surface_mesh->point(vert_from)).length();

				// iterate over one-ring, find incident halfedge which has invalid face handle
				int num_invalid = 0;
				for(VertexOHalfedgeIter voh_it = surface_mesh->voh_iter(vert_to); voh_it; ++voh_it)  {
					OpenMesh::HalfedgeHandle next_he = voh_it.handle();
					OpenMesh::FaceHandle face_test = surface_mesh->face_handle(next_he);
					if(!surface_mesh->is_valid_handle(face_test))  {
						he = next_he;
						num_invalid++;
					}
				}
				at_beginning = he.idx() == start_he.idx();

				if(num_invalid >= 2)  {
					cerr << "BAD MESH! NONMANIFOLD!" << " : " << num_invalid << endl;
					is_manifold = false;
					break;
				}
				if(num_traced >= input_mesh->n_halfedges())  {
					cerr << "OpenMesh didn't pick it up, but this ain't manifold..." << endl;
					is_manifold = false;
					break;
				}
			}
			if(!is_manifold)
				break;
		}

		num_boundary_components = boundary_id;
		int euler_characteristic = surface_mesh->n_vertices() - surface_mesh->n_edges() + surface_mesh->n_faces();
		genus = (2 - euler_characteristic - num_boundary_components) / 2;
		delete [] he_ids;
	}
	else  {
		genus = -1;
		boundary_length = -1;
		num_boundary_components = -1;
	}

	// finally, compute surface area on extracted surface mesh, and input surface mesh
	int num_total_degenerate;
	this->surface_area(input_mesh, total_surface_area, num_total_degenerate);
	this->surface_area(surface_mesh, mesh_surface_area, num_degenerate_tris);

	cout << "num vertices: " << surface_mesh->n_vertices() <<"/"<< input_mesh->n_vertices() <<
		" : num faces: " << surface_mesh->n_faces() <<"/"<< input_mesh->n_faces() << endl;
	cout << "is_manifold ? " << is_manifold << endl;
	cout << "number of connected components: " << num_connected_components << endl;
	cout << "genus: " << genus << endl;
	cout << "number of boundary components: " << num_boundary_components << " ; boundary length: " << boundary_length << endl;
	cout << "total degenerate tris: " << num_total_degenerate << endl;
	cout << "num degenerate tris: " << num_degenerate_tris << endl;
	cout << "total surface area: " << total_surface_area << endl;
	cout << "mesh surface area: " << mesh_surface_area << endl;

	delete [] cc_sizes;
	delete [] vertex_ccs;
	delete [] face_ccs;
}

void ReconEvaluation::getVertices(BasicTriMesh* _mesh, OpenMesh::FaceHandle _fHandle, int& v0, int& v1, int& v2)  {
	FaceVertexIter fv_it = _mesh->fv_iter(_fHandle);
	v0 = fv_it.handle().idx();
	v1 = (++fv_it).handle().idx();
	v2 = (++fv_it).handle().idx();
}

double ReconEvaluation::surface_area(BasicTriMesh* _mesh, int _tri)  {
	OpenMesh::FaceHandle face(_tri);
	if(!_mesh->is_valid_handle(face))
		return 0;
	int v0, v1, v2;

	this->getVertices(_mesh, face, v0, v1, v2);
	Point p0 = _mesh->point(OpenMesh::VertexHandle(v0));
	Point p1 = _mesh->point(OpenMesh::VertexHandle(v1));
	Point p2 = _mesh->point(OpenMesh::VertexHandle(v2));

	if(p0 == p1 || p0 == p2 || p1 == p2)
		return 0;
	Point p_normal = (p1-p0) % (p2-p0);
	return 0.5*p_normal.length();
}

void ReconEvaluation::surface_area(BasicTriMesh* _mesh, double& area, int& num_degenerate)  {
	area = 0;
	num_degenerate = 0;

	for(int t = 0; t < _mesh->n_faces(); t++)  {
		OpenMesh::FaceHandle next_face(t);
		if(!_mesh->is_valid_handle(next_face))
			continue;
		int v0, v1, v2;

		this->getVertices(_mesh, next_face, v0, v1, v2);
		Point p0 = _mesh->point(OpenMesh::VertexHandle(v0));
		Point p1 = _mesh->point(OpenMesh::VertexHandle(v1));
		Point p2 = _mesh->point(OpenMesh::VertexHandle(v2));

		if(p0 == p1 || p0 == p2 || p1 == p2)  {
			num_degenerate++;
			continue;
		}
		Point p_normal = (p1-p0) % (p2-p0);
		area += 0.5*p_normal.length();
	}
}

void ReconEvaluation::evaluate_all()  {
	double total_dot = 0;
	Vector3 shape_ll, shape_ur;
	shape->bounds(shape_ll, shape_ur);

	// determine whether to invert normals
	for(int p = 0; p < n_mesh_points; p++)  {
		if(!determine_inversion)
			break;
		Vector3 mesh_pt = sampled_mesh[p];
		if(!shape->inBound(mesh_pt))
			continue;
		int tri_ind = mesh_tree->getTriMapping(p);

		FaceVertexIter fv_it = surface_mesh->fv_iter(OpenMesh::FaceHandle(tri_ind));
		OpenMesh::VertexHandle v0 = fv_it.handle();
		OpenMesh::VertexHandle v1 = (++fv_it).handle();
		OpenMesh::VertexHandle v2 = (++fv_it).handle();
		Point p0 = surface_mesh->point(v0), p1 = surface_mesh->point(v1), p2 = surface_mesh->point(v2);
		Point p_normal = (p1-p0) % (p2-p0);

		if(p0 == p1 || p0 == p2 || p1 == p2)
			continue;
		double area = 0.5*p_normal.length();
		if(area < AREA_EPSILON)
			continue;

		Vector3 normal(p_normal[0],p_normal[1],p_normal[2]);
		normal.normalize();

		Vector3 shape_normal = shape->normal(mesh_pt);
		total_dot += shape_normal.dotProduct(normal);
	}
	to_invert = total_dot < 0;
	cout << "to invert ? " << total_dot << " : " << to_invert << endl;

	ComputationTimer timer("reconstruction evaluation");
	timer.start();
	cerr << "compute shape correspondence" << endl;
	this->compute_shape_correspondence();
	cerr << "compute mesh correspondence" << endl;
	this->compute_mesh_correspondence();
	timer.end();
	cout << timer.getComputation() << " : " << timer.getElapsedTime() << "s" << endl;
}

void ReconEvaluation::write_evaluation(string _evalBase, double _computationTime, bool _writeCorrespondences)  {
	string global_filename = _evalBase + ".recon";
	string distribution_filename = _evalBase + ".dist";

	// write out the correspondences?
	if(_writeCorrespondences)  {
		string m2i_filename = _evalBase + ".m2i";
		string i2m_filename = _evalBase + ".i2m";
		mesh_to_shape->write_to_file(m2i_filename);
		shape_to_mesh->write_to_file(i2m_filename);
	}

	ShapeDistribution shape_distribution(shape_to_mesh, mesh_to_shape);
	shape_distribution.write_to_file(distribution_filename);

	GlobalStats global_stats;
	global_stats.num_connected_components = num_connected_components;
	global_stats.is_manifold = is_manifold;
	global_stats.num_boundary_components = num_boundary_components;
	global_stats.boundary_length = boundary_length;
	global_stats.genus = genus;
	global_stats.num_degenerate_tris = num_degenerate_tris;
	global_stats.total_surface_area = total_surface_area;
	global_stats.mesh_surface_area = mesh_surface_area;
	global_stats.computation_time = _computationTime;
	global_stats.write_stats(global_filename);
}

// problem framed as: for each point p_s on the shape, cast a ray in its normal direction to
// obtain point p_m on mesh : now have correspondence < p_m -> p_s >
void ReconEvaluation::compute_shape_correspondence()  {
	int num_prints = 20;
	int num_hits = 0;
	int num_invalid = 0;
	int div_prints = n_shape_pts/num_prints;
	for(int p = 0; p < n_shape_pts; p++)  {
		Vector3 shape_pt = sampled_shape[p];
		Vector3 normal = shape->normal(shape_pt);
		Vector3 inv_normal(normal);
		inv_normal.scale(-1.0);

		// cast ray in front & back to the triangle mesh
		Vector3 front_pt, back_pt, front_normal, back_normal, best_pt, best_normal;
		bool intersected_front = kd_tree->intersect(shape_pt, normal, front_pt, front_normal);
		bool intersected_back = kd_tree->intersect(shape_pt, inv_normal, back_pt, back_normal);
		bool hit_surface = true;

		num_hits++;
		if(intersected_front && !intersected_back)  {
			best_pt = front_pt;
			best_normal = front_normal;
		}
		else if(!intersected_front && intersected_back)  {
			best_pt = back_pt;
			best_normal = back_normal;
		}
		else if(intersected_front && intersected_back)  {
			double front_dist = shape_pt.length(front_pt);
			double back_dist = shape_pt.length(back_pt);
			if(front_dist < back_dist)  {
				best_pt = front_pt;
				best_normal = front_normal;
			}
			else  {
				best_pt = back_pt;
				best_normal = back_normal;
			}
		}
		else  {
			hit_surface = false;
			num_hits--;
		}

		// if we didn't hit the surface, then take closest point on mesh, and
		// add it to shortest distance map of shape-to-mesh
		if(!hit_surface)  {
			Vector3 closest_mesh_pt, closest_mesh_normal;
			mesh_tree->closest_point_normal_query(shape_pt, closest_mesh_pt, closest_mesh_normal);
			if(to_invert)
				closest_mesh_normal.scale(-1.0);
			shape_to_mesh->add_correspondence(shape_pt, closest_mesh_pt, shape->normal(shape_pt), closest_mesh_normal);
		}
		else  {
			// check for nearest neighbor on shape
			Vector3 closest_shape_pt = shape_tree->closest_point_query(best_pt);
			double best_distance = best_pt.length(shape_pt);
			double discrete_distance = best_pt.length(closest_shape_pt);

			// if the nearest neighbor is a better distance than the "best" distance, then
			// take closest point to shape_pt on mesh, and add to shape-to-mesh shortest distance map
			if(discrete_distance < best_distance)  {
				Vector3 closest_mesh_pt, closest_mesh_normal;
				mesh_tree->closest_point_normal_query(shape_pt, closest_mesh_pt, closest_mesh_normal);
				if(to_invert)
					closest_mesh_normal.scale(-1.0);
				shape_to_mesh->add_correspondence(shape_pt, closest_mesh_pt, shape->normal(shape_pt), closest_mesh_normal);
			}

			// otherwise, add <best_pt, shape_pt> to mesh-to-shape shortest distance map
			else  {
				if(to_invert)
					best_normal.scale(-1.0);
				mesh_to_shape->add_correspondence(best_pt, shape_pt, best_normal, normal);
			}
		}

		if((p % div_prints) == 0)
			cerr << "shape correspondence["<<p<<"]" << endl;
	}

	cout << "percentage of mesh hit: " << ((double)num_hits / (double)n_shape_pts)  << endl;
}

// problem framed as: for each point p_m on the mesh, cast a ray in its normal direction to
// obtain point p_s on shape : now have correspondence < p_s -> p_m >
void ReconEvaluation::compute_mesh_correspondence()  {
	int num_prints = 50;
	int num_hits = 0;
	int num_invalid = 0;
	int div_prints = std::max(1, n_mesh_points/num_prints);

	int valid_tris = 0;
	bool is_signed_distance_inverted = shape->is_signed_distance_inverted();


	for(int p = 0; p < n_mesh_points; p++)  {
		Vector3 mesh_pt = sampled_mesh[p];
		int tri_ind = mesh_tree->getTriMapping(p);

		FaceVertexIter fv_it = surface_mesh->fv_iter(OpenMesh::FaceHandle(tri_ind));
		OpenMesh::VertexHandle v0 = fv_it.handle();
		OpenMesh::VertexHandle v1 = (++fv_it).handle();
		OpenMesh::VertexHandle v2 = (++fv_it).handle();
		Point p0 = surface_mesh->point(v0), p1 = surface_mesh->point(v1), p2 = surface_mesh->point(v2);
		Point p_normal = (p1-p0) % (p2-p0);
		if(to_invert)  {
			for(int i = 0; i < 3; i++)
				p_normal[i]*=-1;
		}
		if(p0 == p1 || p0 == p2 || p1 == p2)
			continue;
		double area = 0.5*p_normal.length();
		if(area < AREA_EPSILON)
			continue;
		valid_tris++;

		Vector3 normal = Vector3(p_normal[0],p_normal[1],p_normal[2]);
		normal.normalize();
		Vector3 inv_normal(normal);
		inv_normal.scale(-1);

		// cast ray in front & back to the implicit surface
		Vector3 best_pt;
		double shape_func = !shape->inBound(mesh_pt) ? 1 : shape->function(mesh_pt);
		bool is_outside = is_signed_distance_inverted ? shape_func <= 0 : shape_func >= 0;
		bool hit_surface = is_outside ? root_finder->intersection(mesh_pt, inv_normal, best_pt) :
														root_finder->intersection(mesh_pt, normal, best_pt);
		num_hits++;

		// if we didn't hit the surface, then take the closest point on the shape,
		// add this correspondence to shortest distance map of mesh-to-shape
		if(!hit_surface)  {
			num_hits--;
			Vector3 closest_shape_pt, closest_shape_normal;
			shape_tree->closest_point_normal_query(mesh_pt, closest_shape_pt, closest_shape_normal);
			mesh_to_shape->add_correspondence(mesh_pt, closest_shape_pt, normal, closest_shape_normal);
		}
		else  {
			// check for nearest neighbor on mesh
			Vector3 closest_mesh_pt = mesh_tree->closest_point_query(best_pt);
			double best_distance = best_pt.length(mesh_pt);
			double discrete_distance = best_pt.length(closest_mesh_pt);

			// if the nearest neighbor is a better distance than the "best" distance, then
			// take closest point to mesh_pt on the shape, and add to mesh-to-shape shortest distance map
			if(discrete_distance < best_distance)  {
				Vector3 closest_shape_pt, closest_shape_normal;
				shape_tree->closest_point_normal_query(mesh_pt, closest_shape_pt, closest_shape_normal);
				mesh_to_shape->add_correspondence(mesh_pt, closest_shape_pt, normal, closest_shape_normal);
			}

			// otherwise add <best_pt, mesh_pt> to shape-to-mesh shortest distance map
			else
				shape_to_mesh->add_correspondence(best_pt, mesh_pt, shape->normal(best_pt), normal);
		}

		if((p % div_prints) == 0)
			cerr << "mesh correspondence["<<p<<"]" << endl;
	}

	cout << "percentage of shape hit: " << ((double)num_hits / (double)valid_tris)  << endl;
}
