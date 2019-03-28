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

#include "../ShortestDistanceMap.h"
#include "../GlobalStats.h"

#include <algorithm>

int main(int argc, char **argv) {
	if(argc < 6)  {
		std::cerr << "usage: " << argv[0] << " eval_base num_pcs table shape genus [shape genus]*" << std::endl;
		return 1;
	}
	int arg_num = 1;
	string eval_base = argv[arg_num++];
	int num_pcs = atoi(argv[arg_num++]);
	string table_filename = argv[arg_num++];
	FILE* table_file = fopen(table_filename.c_str(), "w");

	vector<string> all_algorithms;
	all_algorithms.push_back("apss");
	all_algorithms.push_back("fourier");
	all_algorithms.push_back("imls");
	all_algorithms.push_back("mpu");
	all_algorithms.push_back("mpusmooth");
	all_algorithms.push_back("poisson");
	all_algorithms.push_back("rbf");
	all_algorithms.push_back("scattered");
	all_algorithms.push_back("spss");
	all_algorithms.push_back("wavelet");

	vector<string> all_shapes;
	vector<int> all_genus;
	while(arg_num < argc)  {
		string next_shape = argv[arg_num++];
		int next_genus = atoi(argv[arg_num++]);
		all_shapes.push_back(next_shape);
		all_genus.push_back(next_genus);
	}
	string exp = "laser";

	fprintf(table_file, "\\begin{tabular}{|l|l|l|l|l|l|l|}\n");
	fprintf(table_file, "\\hline\n");
	fprintf(table_file, "algorithm & comps & bndry & manifold & genus & time \\\\ \\hline\n");
	for(unsigned a = 0; a < all_algorithms.size(); a++)  {
		string algorithm = all_algorithms[a];

		double ave_ccs = 0;
		double ave_boundary = 0;
		double ave_is_manifold = 0; 
		double ave_genus = 0; 
		double ave_comp = 0;

		int valid_shapes = 0;
		int valid_manifold = 0;

		for(unsigned s = 0; s < all_shapes.size(); s++)  {
			int shape_genus = all_genus[s];
			string shape = all_shapes[s];
			for(unsigned p = 0; p < num_pcs; p++)  {
				char pc_num[30];
				sprintf(pc_num, "%u", p);
				string next_shape = eval_base+"/"+shape+"/"+exp+"/"+algorithm+"/"+shape+"_"+pc_num+".recon";
				GlobalStats global_stats(next_shape);
				if(!global_stats.is_valid_file())  {
					cerr << "badness: " << next_shape << endl;
					continue;
				}

				valid_shapes++;
				ave_ccs += global_stats.num_connected_components;
				ave_comp += global_stats.computation_time;
				ave_is_manifold += global_stats.is_manifold ? 1 : 0;

				if(global_stats.is_manifold)  {
					valid_manifold++;
					ave_boundary += global_stats.boundary_length;
					ave_genus += fabs(global_stats.genus-shape_genus);
				}
			}
		}

		ave_ccs /= (double)valid_shapes;
		ave_comp /= (double)valid_shapes;
		ave_is_manifold /= (double)valid_shapes;

		ave_boundary /= (double)valid_manifold;
		ave_genus /= (double)valid_manifold;

		fprintf(table_file, "%s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\ \\hline\n",
				algorithm.c_str(), ave_ccs, ave_boundary, ave_is_manifold, ave_genus, ave_comp);
	}
	fprintf(table_file, "\\hline\n");
	fprintf(table_file, "\\end{tabular}\n");

	fclose(table_file);
}
