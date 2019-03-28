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

#ifndef POINTCLOUD_H
#define POINTCLOUD_H

#include <map>
#include <algorithm>

#include "OrientedPointCloud.h"
#include "union_find.h"
#include "basic_graph.h"
#include "../modeling/DenseEigenSystem.h"
#include <ANN/ANN.h>

using namespace std;

struct OrientationEdge  {
	OrientationEdge()  {}
	OrientationEdge(int _v1, int _v2, double _weight) : v1(_v1), v2(_v2), weight(_weight)  {}

	inline bool operator<(const OrientationEdge& _otherEdge) const  {
		return weight < _otherEdge.weight;
	}

	int v1, v2;
	double weight;
};

typedef std::pair<int, int> knn_edge;

class PointCloud  {
	public:
		PointCloud();
		~PointCloud();

		void addPoint(Vector3 _point)  {
			points.push_back(_point);
		}

		Vector3 getPoint(int _i)  { return points[_i]; }
		void clear()  {
			points.clear();
		}

		int size()  { return points.size(); }

		void write_to_file(string _filename)  {
			FILE* file = fopen(_filename.c_str(), "w");
			for(int i = 0; i < points.size(); i++)  {
				Vector3 pt = points[i];
				if(i == (points.size()-1))
					fprintf(file, "%.7f %.7f %.7f", pt.x,pt.y,pt.z);
				else
					fprintf(file, "%.7f %.7f %.7f\n", pt.x,pt.y,pt.z);
			}
			fclose(file);
		}

		OrientedPointCloud* orient_points(int _knn);
		int gather_neighbors(Vector3 _query, int _numPts, int* neighbor_inds);

	private:
		vector<Vector3> points;

		bool has_kdtree;
		ANNpointArray point_array;
		ANNkd_tree* kd_tree;
};

#endif
