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

#include "GlobalStats.h"

GlobalStats::GlobalStats()  {
	is_valid = true;
}

GlobalStats::GlobalStats(string _infile)  {
	is_valid = true;
	FILE* file_in = fopen(_infile.c_str(), "rb");
	if(file_in == 0)  {
		is_valid = false;
		return;
	}
	int manifold_int;
	size_t num_read;

	num_read = fread(&num_connected_components, sizeof(int), 1, file_in);
	num_read = fread(&manifold_int, sizeof(int), 1, file_in);
	num_read = fread(&num_boundary_components, sizeof(int), 1, file_in);
	num_read = fread(&boundary_length, sizeof(double), 1, file_in);
	num_read = fread(&genus, sizeof(int), 1, file_in);
	num_read = fread(&num_degenerate_tris, sizeof(int), 1, file_in);
	num_read = fread(&total_surface_area, sizeof(double), 1, file_in);
	num_read = fread(&mesh_surface_area, sizeof(double), 1, file_in);
	num_read = fread(&computation_time, sizeof(double), 1, file_in);

	is_manifold = manifold_int == 1;

	/*
	cout << "num connected components: " << num_connected_components << endl;
	cout << "is manifold ? " << is_manifold << endl;
	cout << "num boundary components: " << num_boundary_components << endl;
	cout << "boundary length: " << boundary_length << endl;
	cout << "genus: " << genus << endl;
	cout << "num degenerate tris: " << num_degenerate_tris << endl;
	cout << "total surface area: " << total_surface_area << endl;
	cout << "mesh surface area: " << mesh_surface_area << endl;
	cout << "computation time: " << computation_time << endl;
	*/

	fclose(file_in);
}

GlobalStats::~GlobalStats()  {
}

void GlobalStats::write_stats(string _outfile)  {
	FILE* file_out = fopen(_outfile.c_str(), "wb");
	int manifold_int = is_manifold ? 1 : 0;

	fwrite(&num_connected_components, sizeof(int), 1, file_out);
	fwrite(&manifold_int, sizeof(int), 1, file_out);
	fwrite(&num_boundary_components, sizeof(int), 1, file_out);
	fwrite(&boundary_length, sizeof(double), 1, file_out);
	fwrite(&genus, sizeof(int), 1, file_out);
	fwrite(&num_degenerate_tris, sizeof(int), 1, file_out);
	fwrite(&total_surface_area, sizeof(double), 1, file_out);
	fwrite(&mesh_surface_area, sizeof(double), 1, file_out);
	fwrite(&computation_time, sizeof(double), 1, file_out);

	fclose(file_out);
}
