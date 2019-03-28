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

#ifndef SPARSERANGESCAN_H
#define SPARSERANGESCAN_H

#include "PointCloud.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

struct PixelCoord  {
	PixelCoord()  {}
	PixelCoord(int _x, int _y) : x(_x) , y(_y) {}
	int x, y;
};

class SparseRangeScan  {
	public:
		SparseRangeScan(string _file);
		SparseRangeScan()  {}
		~SparseRangeScan();

		int size()  { return pixels.size(); }
		PixelCoord getPixel(int _i)  { return pixels[_i]; }
		Vector3 getPt(int _i)  { return pts[_i]; }

		void addPoint(PixelCoord _pix, Vector3 _pt)  {
			pixels.push_back(_pix);
			pts.push_back(_pt);
		}

		void write_to_file(string _filename)  {
			FILE* range_out = fopen(_filename.c_str(), "w");
			for(unsigned i = 0; i < this->size(); i++)  {
				PixelCoord pix = pixels[i];
				Vector3 pt = pts[i];
				fprintf(range_out, "%u %u %.6f %.6f %.6f\n", pix.x, pix.y, pt.x, pt.y, pt.z);
			}
			fclose(range_out);
		}

	private:
		vector<PixelCoord> pixels;
		vector<Vector3> pts;

		void tokenize_line(vector<string>* tokens, const string& input, string sep);
};

#endif
