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

#ifndef ORIENTEDPOINTCLOUD_H
#define ORIENTEDPOINTCLOUD_H

#include <stdio.h>
#include <vector>
#include "../modeling/Vector3.h"

using namespace std;

class OrientedPointCloud  {
	public:
		OrientedPointCloud();
		~OrientedPointCloud();

		void addOrientedPoint(Vector3 _point, Vector3 _normal)  {
			points.push_back(_point);
			normals.push_back(_normal);
		}

		Vector3 getPoint(int _i)  { return points[_i]; }
		Vector3 getNormal(int _i)  { return normals[_i]; }
		void getOrientedPoint(int _i, Vector3& point, Vector3& normal)  {
			point = points[_i];
			normal = normals[_i];
		}
		void clear()  {
			points.clear();
			normals.clear();
		}

		int size()  { return points.size(); }

		void write_to_file(string _filename)  {
			FILE* pts_file = fopen(_filename.c_str(), "w");
			for(int i = 0; i < points.size(); i++)  {
				Vector3 pt = points[i];
				Vector3 normal = normals[i];
				fprintf(pts_file, "%.7f %.7f %.7f %.7f %.7f %.7f\n", pt.x, pt.y, pt.z, normal.x, normal.y, normal.z);
			}
			fclose(pts_file);
		}

	private:
		vector<Vector3> points;
		vector<Vector3> normals;
};

#endif
