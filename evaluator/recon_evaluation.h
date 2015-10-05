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

#ifndef RECONEVALUATION_H
#define RECONEVALUATION_H

#include "shape_distribution.h"
#include "ComputationTimer.h"
#include "GlobalStats.h"
#include "KdOpenMesh.h"
#include "kd_tree.h"
#include "../modeling/ImplicitFunction.h"
#include "../sampler/OrientedPointCloud.h"
#include "../sampler/ImplicitNewton.h"

#include "mesh_sampler.h"

#define AREA_EPSILON 1e-8
#define BOUND_EPSILON 0.01
#define MIN_MESH_SAMPLES 400000
#define MAX_MESH_SAMPLES 600000

class ReconEvaluation  {
	private:
		typedef BasicTriMesh::Point Point;
		typedef BasicTriMesh::VertexIter VertexIter;
		typedef BasicTriMesh::EdgeIter EdgeIter;
		typedef BasicTriMesh::HalfedgeIter HalfedgeIter;

		typedef BasicTriMesh::VertexVertexIter VertexVertexIter;
		typedef BasicTriMesh::FaceFaceIter FaceFaceIter;
		typedef BasicTriMesh::FaceVertexIter FaceVertexIter;
		typedef BasicTriMesh::VertexOHalfedgeIter VertexOHalfedgeIter;
		typedef BasicTriMesh::VertexIHalfedgeIter VertexIHalfedgeIter;

	public:
		ReconEvaluation(ImplicitFunction* _shape, OrientedPointCloud* _pointShape, BasicTriMesh* _surfaceMesh);
		~ReconEvaluation();

		// batch evaluation
		void evaluate_all();
		void write_evaluation(string _evalBase, double _computationTime, bool _writeCorrespondences);

		// correspondences
		void compute_mesh_correspondence();
		void compute_shape_correspondence();

		BasicTriMesh* getInputMesh()  { return input_mesh; }
		BasicTriMesh* getSurfaceMesh()  { return surface_mesh; }

		ShortestDistanceMap* getMeshToShapeCorrespondence()  { return mesh_to_shape; }
		ShortestDistanceMap* getShapeToMeshCorrespondence()  { return shape_to_mesh; }

		void enableInversion()  { determine_inversion = true; }

	private:
		// handle meta-data
		void pre_process();
		void getVertices(BasicTriMesh* _mesh, OpenMesh::FaceHandle _fHandle, int& v0, int& v1, int& v2);
		void surface_area(BasicTriMesh* _mesh, double& area, int& num_degenerate);
		double surface_area(BasicTriMesh* _mesh, int _tri);

		// analytical shape
		ImplicitFunction* shape;

		// particle system sampling of shape
		OrientedPointCloud* point_shape;

		// input mesh & corresponding surface mesh (i.e. largest connected component)
		BasicTriMesh* input_mesh;
		BasicTriMesh* surface_mesh;

		// bounds of shape
		Vector3 ll_bound, ur_bound, diag;
		double bound_eps;
		double lipschitz;

		// representative points for shape & mesh
		Vector3* sampled_shape;
		Vector3* sampled_mesh;

		// correspondences
		ShortestDistanceMap* mesh_to_shape;
		ShortestDistanceMap* shape_to_mesh;

		// utilities
		KdOpenMesh* kd_tree;
		ImplicitNewton* root_finder;
		ImplicitKdTree* shape_tree;
		MeshKdTree* mesh_tree;

		// meta-data
		int num_connected_components;
		bool is_manifold;
		int num_boundary_components;
		double boundary_length;
		int genus;
		int num_degenerate_tris;
		double total_surface_area;
		double mesh_surface_area;
		double computation_time;

		int n_shape_pts;
		int n_mesh_vertices;
		int n_mesh_tris;
		int n_mesh_points;

		bool determine_inversion;
		bool to_invert;
};

#endif
