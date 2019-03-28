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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

int main(int argc, char** argv)  {
	if(argc != 3)  {
		std::cerr << "usage: " << argv[0] << " dist_file out" << std::endl;
		return 1;
	}
	int arg_num = 1;
	string dist_filename = argv[arg_num++];
	string out_filename = argv[arg_num++];

	FILE* dist_file = fopen(dist_filename.c_str(), "rb");
	size_t num_read = 0;
	double shape_dist[8];
	num_read = fread(shape_dist, sizeof(double), 8, dist_file);
	double angle_dist[8];
	num_read = fread(angle_dist, sizeof(double), 8, dist_file);
	fclose(dist_file);

	FILE* out_file = fopen(out_filename.c_str(), "w");
	fprintf(out_file, "%.7f %.7f %.7f %.7f %.7f %.7f %.7f %.7f %.7f %.7f\n", shape_dist[0], shape_dist[1], shape_dist[2], shape_dist[3], shape_dist[4],
																									  angle_dist[0], angle_dist[1], angle_dist[2], angle_dist[3], angle_dist[4]);
	fclose(out_file);
}
