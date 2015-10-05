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

#include "PointCloud.h"

PointCloud::PointCloud()  {
	has_kdtree = false;
}

PointCloud::~PointCloud()  {
	if(has_kdtree)  {
		delete kd_tree;
		for(unsigned i = 0; i < points.size(); i++)
			delete point_array[i];
	}
}

OrientedPointCloud* PointCloud::orient_points(int _knn)  {
	// perform PCA to obtain centers & unoriented normals
	Vector3* centers = new Vector3[points.size()];
	Vector3* unorg_normals = new Vector3[points.size()];
	map<knn_edge, bool> edge_map;

	cout << "performing pca, gathering knn ..." << endl;
	for(unsigned v = 0; v < points.size(); v++)  {
		Vector3 point = points[v];
		double c_xx, c_xy, c_xz, c_yy, c_yz, c_zz;
		c_xx = c_xy = c_xz = c_yy = c_yz = c_zz = 0;

		int* neighbor_inds = new int[_knn];
		this->gather_neighbors(point, _knn, neighbor_inds);
		Vector3 center_point;
		for(int k = 0; k < _knn; k++)  {
			center_point = center_point + points[ neighbor_inds[k] ];
			knn_edge new_edge(v, neighbor_inds[k]), inv_edge(neighbor_inds[k], v);
			if(edge_map.find(inv_edge) != edge_map.end())
				continue;
			edge_map[new_edge] = true;
		}
		center_point.scale((1.0 / (double)_knn));

		for(int k = 0; k < _knn; k++)  {
			Vector3 centered_pt = points[ neighbor_inds[k] ] - center_point;
			c_xx += centered_pt.x*centered_pt.x;
			c_xy += centered_pt.x*centered_pt.y;
			c_xz += centered_pt.x*centered_pt.z;
			c_yy += centered_pt.y*centered_pt.y;
			c_yz += centered_pt.y*centered_pt.z;
			c_zz += centered_pt.z*centered_pt.z;
		}

		DenseMatrix covariance(c_xx, c_xy, c_xz, c_yy, c_yz, c_zz);
		DenseEigenSystem pca(&covariance, 0, 2);
		pca.solve();

		BLinAlg::Vector* eigenvector0 = pca.getEigenvector(0);
		BLinAlg::Vector* eigenvector1 = pca.getEigenvector(1);
		BLinAlg::Vector* eigenvector2 = pca.getEigenvector(2);
		Vector3 unorg_normal = Vector3(eigenvector0->getEntry(0), eigenvector0->getEntry(1), eigenvector0->getEntry(2));
		Vector3 tangent1 = Vector3(eigenvector1->getEntry(0), eigenvector1->getEntry(1), eigenvector1->getEntry(2));
		Vector3 tangent2 = Vector3(eigenvector2->getEntry(0), eigenvector2->getEntry(1), eigenvector2->getEntry(2));

		centers[v] = center_point;
		unorg_normals[v] = unorg_normal;

		delete [] neighbor_inds;
	}
	cout << "... done performing pca" << endl;

	// now construct MST - first gather weighted edges of knn graph 
	cout << "constructing orientation edges..." << endl;
	OrientationEdge* orientation_edges = new OrientationEdge[edge_map.size()];
	int num_edges = 0;
	for(std::map<knn_edge, bool>::iterator edge_iter = edge_map.begin() ; edge_iter != edge_map.end(); edge_iter++)  {
		knn_edge next_edge = (*edge_iter).first;
		double normal_dot = unorg_normals[next_edge.first].dotProduct(unorg_normals[next_edge.second]);
		double edge_weight = normal_dot < 0 ? (1.0+normal_dot) : (1.0-normal_dot);
		if(edge_weight < 0)
			edge_weight = 0;
		orientation_edges[num_edges++] = OrientationEdge(next_edge.first, next_edge.second, edge_weight);
	}
	edge_map.clear();
	cout << "... done constructing orientation edges..." << endl;
	std::sort(orientation_edges, orientation_edges+num_edges);
	cout << "... done sorting orientation edges..." << endl;

	// setup union-find set data structures, for MST algorithm
	UFEntity<int>** all_sets = new UFEntity<int>*[points.size()];
	for(unsigned i = 0; i < points.size(); i++)
		all_sets[i] = new UFEntity<int>(i);
	cout << "... done initializing union-find" << endl;

	// construct MST
	BasicGraph mst(points.size());
	bool tree_complete = false;
	int edge_iter = 0, tree_size = 0;
	UFSet<int> set_operations;
	cout << "constructing mst ..." << endl;
	while(!tree_complete)  {
		OrientationEdge next_edge = orientation_edges[edge_iter++];
		UFEntity<int>* node1 = all_sets[next_edge.v1];
		UFEntity<int>* node2 = all_sets[next_edge.v2];
		if(set_operations.perform_union(node1, node2))  {
			mst.add_edge(next_edge.v1, next_edge.v2);
			mst.add_edge(next_edge.v2, next_edge.v1);
			tree_size++;
		}
		tree_complete = (tree_size == (points.size()-1)) || (edge_iter == num_edges);
	}
	cout << "... done constructing mst" << endl;
	delete [] orientation_edges;
	for(unsigned i = 0; i < points.size(); i++)
		delete all_sets[i];
	delete [] all_sets;

	bool* is_visited = new bool[points.size()];
	Vector3* oriented_normals = new Vector3[points.size()];
	for(unsigned i = 0; i < points.size(); i++)
		is_visited[i] = false;

	// perform MST on each separate component of knn graph
	int total_visited = 0;
	while(total_visited != points.size())  {

		// start the propagation by first initializing the normal whose center is max positive z to be positive
		int init_pt = -1;
		double max_z = -1e10;
		for(unsigned i = 0; i < points.size(); i++)  {
			if(points[i].z > max_z && !is_visited[i])  {
				init_pt = i;
				max_z = points[i].z;
			}
		}
		oriented_normals[init_pt] = unorg_normals[init_pt];
		if(oriented_normals[init_pt].z < 0)
			oriented_normals[init_pt].scale(-1);

		// perform propagation
		is_visited[init_pt] = true;
		total_visited++;
		vector<int> parent_stack;
		parent_stack.push_back(init_pt);
		while(parent_stack.size() > 0)  {
			// try to find an unvisited neighbor
			int parent_pt = parent_stack[parent_stack.size()-1];
			AdjacencyList* neighbors = mst.get_adjacencies(parent_pt);
			int child_pt = -1;
			for(int n = 0; n < neighbors->size(); n++)  {
				int next_pt = neighbors->get_adjacency(n);
				if(!is_visited[next_pt])  {
					child_pt = next_pt;
					break;
				}
			}

			// if a child point is not found, then pop the parent from the stack
			if(child_pt == -1)  {
				parent_stack.erase(parent_stack.begin()+(parent_stack.size()-1));
				continue;
			}

			// otherwise, depending on parent and child orientation, propagate, and add child to stack
			Vector3 parent_normal = oriented_normals[parent_pt];
			Vector3 child_normal = unorg_normals[child_pt];
			oriented_normals[child_pt] = unorg_normals[child_pt];
			if(oriented_normals[child_pt].dotProduct(oriented_normals[parent_pt]) < 0)
				oriented_normals[child_pt].scale(-1);
			is_visited[child_pt] = true;
			total_visited++;
			parent_stack.push_back(child_pt);
		}
	}

	// clean up
	delete [] centers;
	delete [] unorg_normals;

	OrientedPointCloud* oriented_pc = new OrientedPointCloud();
	for(unsigned v = 0; v < points.size(); v++)
		oriented_pc->addOrientedPoint(points[v], oriented_normals[v]);

	// ... further clean up
	delete [] is_visited;
	delete [] oriented_normals;

	return oriented_pc;
}

int PointCloud::gather_neighbors(Vector3 _query, int _numPts, int* neighbor_inds)  {
	if(!has_kdtree)  {
		point_array = annAllocPts(points.size(), 3);
		for(unsigned v = 0; v < points.size(); v++)  {
			Vector3 point = points[v];
			point_array[v] = annAllocPt(3);
			point_array[v][0] = point.x;
			point_array[v][1] = point.y;
			point_array[v][2] = point.z;
		}
		kd_tree = new ANNkd_tree(point_array, points.size(), 3);
		has_kdtree = true;
	}

	int num_pts = _numPts;
	if(_numPts > points.size())  {
		cerr << "kdtree query error: more neighbors than points! using number of points instead ..." << endl;
		num_pts = points.size();
	}

	ANNidxArray nnIdx = new ANNidx[num_pts];
	ANNdistArray dists = new ANNdist[num_pts];

	ANNpoint query_pt = annAllocPt(3);
	query_pt[0] = _query.x;
	query_pt[1] = _query.y;
	query_pt[2] = _query.z;

	kd_tree->annkSearch(query_pt, num_pts, nnIdx, dists);
	for(int i = 0; i < num_pts; i++)
		neighbor_inds[i] = nnIdx[i];

	delete [] nnIdx;
	delete [] dists;
	annDeallocPt(query_pt);

	return num_pts;
}
