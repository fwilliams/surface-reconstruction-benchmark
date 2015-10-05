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

#include "ShortestDistanceMap.h"

ShortestDistanceMap::ShortestDistanceMap()  {
}

ShortestDistanceMap::ShortestDistanceMap(string _infile)  {
	FILE* file_in = fopen(_infile.c_str(), "rb");
	int num;
	size_t num_read;
	num_read = fread(&num, sizeof(int), 1, file_in);

	double* full_buffer = new double[12*num];
	num_read = fread(full_buffer, sizeof(double), (12*num), file_in);
	fclose(file_in);

	double p_buf[6];
	double n_buf[6];
	for(int p = 0; p < num; p++)  {
		double x1 = full_buffer[12*p], y1 = full_buffer[12*p+1], z1 = full_buffer[12*p+2];
		double x2 = full_buffer[12*p+3], y2 = full_buffer[12*p+4], z2 = full_buffer[12*p+5];
		double nx1 = full_buffer[12*p+6], ny1 = full_buffer[12*p+7], nz1 = full_buffer[12*p+8];
		double nx2 = full_buffer[12*p+9], ny2 = full_buffer[12*p+10], nz2 = full_buffer[12*p+11];
		point_correspondences.push_back(PointPair(Vector3(x1,y1,z1), Vector3(x2,y2,z2)));
		normal_correspondences.push_back(NormalPair(Vector3(nx1,ny1,nz1), Vector3(nx2,ny2,nz2)));
	}

	delete [] full_buffer;
}

ShortestDistanceMap::~ShortestDistanceMap()  {
	point_correspondences.clear();
	normal_correspondences.clear();
}

void ShortestDistanceMap::add_correspondence(Vector3 _p1, Vector3 _p2, Vector3 _norm1, Vector3 _norm2)  {
	point_correspondences.push_back(PointPair(_p1,_p2));
	normal_correspondences.push_back(NormalPair(_norm1,_norm2));
}

void ShortestDistanceMap::add_correspondence(Vector3 _p1, Vector3 _p2)  {
	point_correspondences.push_back(PointPair(_p1,_p2));
	normal_correspondences.push_back(NormalPair());
}

void ShortestDistanceMap::write_to_file(string _filename)  {
	FILE* file_out = fopen(_filename.c_str(), "wb");

	int num = this->num_correspondences();
	fwrite(&num, sizeof(int), 1, file_out);
	for(int i = 0; i < num; i++)  {
		PointPair point_pair = point_correspondences[i];
		NormalPair normal_pair = normal_correspondences[i];

		double all_points[] = { point_pair.p1.x, point_pair.p1.y, point_pair.p1.z,
										point_pair.p2.x, point_pair.p2.y, point_pair.p2.z };
		double all_normals[] = { normal_pair.normal1.x, normal_pair.normal1.y, normal_pair.normal1.z,
										normal_pair.normal2.x, normal_pair.normal2.y, normal_pair.normal2.z };
		fwrite(all_points, sizeof(double), 6, file_out);
		fwrite(all_normals, sizeof(double), 6, file_out);
	}

	fclose(file_out);
}
