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

#include "RangeImage.h"

RangeImage::RangeImage(int _resX, int _resY)  {
	all_points = new Vector3[_resX*_resY];
	valid_hits = new bool[_resX*_resY];

	x_res = _resX;
	y_res = _resY;

	for(int y = 0; y < y_res; y++)  {
		for(int x = 0; x < x_res; x++)  {
			all_points[x + x_res*y] = Vector3();
			valid_hits[x + x_res*y] = false;
		}
	}
}

RangeImage::~RangeImage()  {
	delete [] all_points;
	delete [] valid_hits;
}

bool RangeImage::isValidPoint(int _x, int _y)  {
	return valid_hits[_x + x_res*_y];
}

Vector3 RangeImage::getRangePoint(int _x, int _y)  {
	return all_points[_x + x_res*_y];
}

void RangeImage::setRangePoint(int _x, int _y, Vector3 _point)  {
	int pt_ind = _x + x_res*_y;
	valid_hits[pt_ind] = true;
	all_points[pt_ind] = _point;
}

void RangeImage::dump_to_ply(string _filename)  {
	int image_res = x_res*y_res;
	int num_vertices = 0;
	for(int y = 0; y < y_res; y++)  {
		for(int x = 0; x < x_res; x++)  {
			if(this->isValidPoint(x,y))
				num_vertices++;
		}
	}

	FILE* ply_out = fopen(_filename.c_str(), "w");

	// header info
	fprintf(ply_out, "ply\n");
	fprintf(ply_out, "format ascii 1.0\n");
	fprintf(ply_out, "obj_info num_cols %u\n", x_res);
	fprintf(ply_out, "obj_info num_rows %u\n", y_res);
	fprintf(ply_out, "element vertex %u\n", num_vertices);
	fprintf(ply_out, "property float x\n");
	fprintf(ply_out, "property float y\n");
	fprintf(ply_out, "property float z\n");
	fprintf(ply_out, "element range_grid %u\n", image_res);
	fprintf(ply_out, "property list uchar int vertex_indices\n");
	fprintf(ply_out, "end_header\n");

	// write out vertices
	for(int y = 0; y < y_res; y++)  {
		for(int x = 0; x < x_res; x++)  {
			if(this->isValidPoint(x,y))  {
				Vector3 range_pt = this->getRangePoint(x,y);
				fprintf(ply_out, "%.7f %.7f %.7f\n", range_pt.x, range_pt.y, range_pt.z);
			}
		}
	}

	// write out (linearized) grid
	int cur_vert = 0;
	for(int y = 0; y < y_res; y++)  {
		for(int x = 0; x < x_res; x++)  {
			if(this->isValidPoint(x,y))  {
				fprintf(ply_out, "1 %u\n", cur_vert);
				cur_vert++;
			}
			else
				fprintf(ply_out, "0\n");
		}
	}

	fclose(ply_out);
}

int RangeImage::xRes()  {
	return x_res;
}

int RangeImage::yRes()  {
	return y_res;
}
