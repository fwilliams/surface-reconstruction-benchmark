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

#include "kd_tree.h"

PointKdTree::PointKdTree(int _numPoints)  {
	num_points = _numPoints;
	point_array = annAllocPts(num_points, 3);
	is_tree_created = false;
}

PointKdTree::~PointKdTree()  {
	if(is_tree_created)
		delete kd_tree;
	for(int i = 0; i < num_points; i++)
		delete point_array[i];
	delete [] point_array;
}

void PointKdTree::add_point(int _ind, Vector3 _point)  {
	point_array[_ind] = annAllocPt(3);
	point_array[_ind][0] = _point.x;
	point_array[_ind][1] = _point.y;
	point_array[_ind][2] = _point.z;
}

void PointKdTree::finalize_tree()  {
	std::cout << "finalizing kd tree..." << std::endl;
	is_tree_created = true;
	kd_tree = new ANNkd_tree(point_array, num_points, 3);
}

vector<int> PointKdTree::ball_query(Vector3 _query, double _radius)  {
	ANNpoint query_pt = annAllocPt(3);
	query_pt[0] = _query.x;
	query_pt[1] = _query.y;
	query_pt[2] = _query.z;
	double sqd_radius = _radius*_radius;

	int the_k = kd_tree->annkFRSearch(query_pt, sqd_radius, 0);
	ANNidxArray nnIdx = new ANNidx[the_k];
	ANNdistArray dists = new ANNdist[the_k];
	kd_tree->annkFRSearch(query_pt, sqd_radius, the_k, nnIdx, dists);

	vector<int> ball_indices;
	for(int i = 0; i < the_k; i++)
		ball_indices.push_back(nnIdx[i]);

	delete [] nnIdx;
	delete [] dists;

	return ball_indices;
}

Vector3 PointKdTree::closest_point_query(Vector3 _query)  {
	ANNidxArray nnIdx = new ANNidx[1];
	ANNdistArray dists = new ANNdist[1];

	ANNpoint query_pt = annAllocPt(3);
	query_pt[0] = _query.x;
	query_pt[1] = _query.y;
	query_pt[2] = _query.z;

	kd_tree->annkSearch(query_pt, 1, nnIdx, dists);
	ANNpoint found_pt = point_array[nnIdx[0]];

	delete [] nnIdx;
	delete [] dists;

	return Vector3(found_pt[0], found_pt[1], found_pt[2]);
}

int PointKdTree::closest_index_query(Vector3 _query)  {
	ANNidxArray nnIdx = new ANNidx[1];
	ANNdistArray dists = new ANNdist[1];

	ANNpoint query_pt = annAllocPt(3);
	query_pt[0] = _query.x;
	query_pt[1] = _query.y;
	query_pt[2] = _query.z;

	kd_tree->annkSearch(query_pt, 1, nnIdx, dists);
	int found_ind = nnIdx[0];

	delete [] nnIdx;
	delete [] dists;

	return found_ind;
}

void PointKdTree::closest_point_query(Vector3 _query, int _numPoints, double _eps, Vector3* answer, int* inds, int& num_found)  {
	ANNidxArray nnIdx = new ANNidx[_numPoints];
	ANNdistArray dists = new ANNdist[_numPoints];

	ANNpoint query_pt = annAllocPt(3);
	query_pt[0] = _query.x;
	query_pt[1] = _query.y;
	query_pt[2] = _query.z;

	kd_tree->annkSearch(query_pt, _numPoints, nnIdx, dists);

	num_found = 0;
	for(int i = 0; i < _numPoints; i++)  {
		int next_id = nnIdx[i];
		if(next_id == ANN_NULL_IDX)
			break;

		ANNpoint found_pt = point_array[next_id];
		Vector3 next_pt = Vector3(found_pt[0], found_pt[1], found_pt[2]);
		if((next_pt-_query).length() > _eps)
			break;

		answer[i] = Vector3(found_pt[0], found_pt[1], found_pt[2]);
		inds[i] = next_id;
		num_found++;
	}

	delete [] nnIdx;
	delete [] dists;
}

// --- ImplicitKdTree --- //

ImplicitKdTree::ImplicitKdTree(int _numPoints, ImplicitFunction* _shape) : PointKdTree(_numPoints)  {
	shape = _shape;
}

ImplicitKdTree::~ImplicitKdTree()  {
}

void ImplicitKdTree::closest_point_normal_query(Vector3 _query, Vector3& point, Vector3& normal)  {
	ANNidxArray nnIdx = new ANNidx[1];
	ANNdistArray dists = new ANNdist[1];

	ANNpoint query_pt = annAllocPt(3);
	query_pt[0] = _query.x;
	query_pt[1] = _query.y;
	query_pt[2] = _query.z;

	kd_tree->annkSearch(query_pt, 1, nnIdx, dists);
	ANNpoint found_pt = point_array[nnIdx[0]];

	point = Vector3(found_pt[0], found_pt[1], found_pt[2]);
	normal = shape->normal(point);

	delete [] nnIdx;
	delete [] dists;
}

// --- MeshKdTree --- //

MeshKdTree::MeshKdTree(int _numPoints, BasicTriMesh* _mesh) : PointKdTree(_numPoints)  {
	mesh = _mesh;
	tri_mapping = new int[_numPoints];
}

MeshKdTree::~MeshKdTree()  {
	delete [] tri_mapping;
}

void MeshKdTree::add_tri_point(int _ind, int _triInd, Vector3 _point)  {
	point_array[_ind] = annAllocPt(3);
	point_array[_ind][0] = _point.x;
	point_array[_ind][1] = _point.y;
	point_array[_ind][2] = _point.z;
	tri_mapping[_ind] = _triInd;
}

void MeshKdTree::closest_point_normal_query(Vector3 _query, Vector3& point, Vector3& normal)  {
	ANNidxArray nnIdx = new ANNidx[1];
	ANNdistArray dists = new ANNdist[1];

	ANNpoint query_pt = annAllocPt(3);
	query_pt[0] = _query.x;
	query_pt[1] = _query.y;
	query_pt[2] = _query.z;

	kd_tree->annkSearch(query_pt, 1, nnIdx, dists);
	ANNpoint found_pt = point_array[nnIdx[0]];

	point = Vector3(found_pt[0], found_pt[1], found_pt[2]);

	int tri = tri_mapping[nnIdx[0]];
	FaceVertexIter fv_it = mesh->fv_iter(OpenMesh::FaceHandle(tri));
	OpenMesh::VertexHandle v0 = fv_it.handle();
	OpenMesh::VertexHandle v1 = (++fv_it).handle();
	OpenMesh::VertexHandle v2 = (++fv_it).handle();

	Point p0 = mesh->point(v0), p1 = mesh->point(v1), p2 = mesh->point(v2);
	Point raw_normal = (p1-p0) % (p2-p0);
	normal = Vector3(raw_normal[0], raw_normal[1], raw_normal[2]);
	normal.normalize();

	delete [] nnIdx;
	delete [] dists;
}
