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

#ifndef SHORTESTDISTANCEMAP_H
#define SHORTESTDISTANCEMAP_H

#include <vector>
#include <string.h>
#include "../modeling/Vector3.h"

#include <stdio.h>
using namespace std;

struct PointPair  {
	PointPair() {}
	PointPair(Vector3 _p1, Vector3 _p2) : p1(_p1) , p2(_p2) {}

	bool has_nan()  {
		return isnan(p1.x) || isnan(p1.y) || isnan(p1.z) || isnan(p2.x) || isnan(p2.y) || isnan(p2.z);
	}

	double distance()  {
		return p1.length(p2);
	}

	Vector3 p1, p2;
};

struct NormalPair  {
	NormalPair() {}
	NormalPair(Vector3 _norm1, Vector3 _norm2) : normal1(_norm1) , normal2(_norm2) {}

	bool has_nan()  {
		return isnan(normal1.x) || isnan(normal1.y) || isnan(normal1.z) || isnan(normal2.x) || isnan(normal2.y) || isnan(normal2.z);
	}

	double angle()  {
		double dot_product = normal1.dotProduct(normal2);
		if(dot_product > 1)
			dot_product=1;
		if(dot_product < -1)
			dot_product=-1;

		double angle_radians = acos(dot_product);
		return 180.0*(angle_radians/acos(-1));
	}

	Vector3 normal1, normal2;
};

class ShortestDistanceMap  {
	public:
		ShortestDistanceMap();
		ShortestDistanceMap(string _infile);
		~ShortestDistanceMap();

		void add_correspondence(Vector3 _p1, Vector3 _p2, Vector3 _norm1, Vector3 _norm2);
		void add_correspondence(Vector3 _p1, Vector3 _p2);

		int num_correspondences()  { return point_correspondences.size(); }
		PointPair getPointCorrespondence(int _i)  { return point_correspondences[_i]; }
		NormalPair getNormalCorrespondence(int _i)  { return normal_correspondences[_i]; }

		void write_to_file(string _filename);

	private:
		std::vector<PointPair> point_correspondences;
		std::vector<NormalPair> normal_correspondences;
};

#endif
