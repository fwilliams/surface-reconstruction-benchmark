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

#ifndef KDTREE_H
#define KDTREE_H

#include "mesh_types.h"
#include "../modeling/Vector3.h"
#include <ANN/ANN.h>

#include <iostream>

#include "../modeling/ImplicitFunction.h"

class PointKdTree  {
	public:
		PointKdTree(int _numPoints);
		~PointKdTree();

		void add_point(int _ind, Vector3 _point);
		void finalize_tree();

		vector<int> ball_query(Vector3 _query, double _radius);
		Vector3 closest_point_query(Vector3 _query);
		int closest_index_query(Vector3 _query);
		void closest_point_query(Vector3 _query, int _numPoints, double _eps, Vector3* answer, int* inds, int& num_found);
		virtual void closest_point_normal_query(Vector3 _query, Vector3& point, Vector3& normal) = 0;

	protected:
		int num_points;
		ANNpointArray point_array;

		bool is_tree_created;
		ANNkd_tree* kd_tree;
};

class SimpleKdTree : public PointKdTree  {
	public:
		SimpleKdTree(int _numPoints) : PointKdTree(_numPoints)  {}
		~SimpleKdTree()  {}

		virtual void closest_point_normal_query(Vector3 _query, Vector3& point, Vector3& normal)  {}

	private:
};

class ImplicitKdTree : public PointKdTree  {
	public:
		ImplicitKdTree(int _numPoints, ImplicitFunction* _shape);
		~ImplicitKdTree();

		virtual void closest_point_normal_query(Vector3 _query, Vector3& point, Vector3& normal);

	private:
		ImplicitFunction* shape;
};

class MeshKdTree : public PointKdTree  {
	private:
		typedef BasicTriMesh::Point Point;
		typedef BasicTriMesh::FaceVertexIter FaceVertexIter;

	public:
		MeshKdTree(int _numPoints, BasicTriMesh* _mesh);
		~MeshKdTree();

		void add_tri_point(int _ind, int _triInd, Vector3 _point);
		int getTriMapping(int _ind)  { return tri_mapping[_ind]; }
		virtual void closest_point_normal_query(Vector3 _query, Vector3& point, Vector3& normal);

	private:
		BasicTriMesh* mesh;
		int* tri_mapping;
};

#endif
