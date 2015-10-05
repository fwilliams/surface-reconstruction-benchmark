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

int main(int argc, char** argv)  {
	if(argc != 6)  {
		std::cerr << "usage: " << argv[0] << " evaluation_base shape experiment num_shapes data_base" << std::endl;
		return 1;
	}
	int arg_num = 1;
	string eval_base = argv[arg_num++];
	string shape = argv[arg_num++];
	string exp = argv[arg_num++];
	int num_shapes = atoi(argv[arg_num++]);

	string data_base = argv[arg_num++];
	string distmean_filename = data_base + "_distmean.txt";
	string distmax_filename = data_base + "_distmax.txt";
	string anglemean_filename = data_base + "_anglemean.txt";
	string anglemax_filename = data_base + "_anglemax.txt";

	FILE* distmean_file = fopen(distmean_filename.c_str(), "w");
	FILE* distmax_file = fopen(distmax_filename.c_str(), "w");
	FILE* anglemean_file = fopen(anglemean_filename.c_str(), "w");
	FILE* anglemax_file = fopen(anglemax_filename.c_str(), "w");

	vector<string> all_algorithms;
	all_algorithms.push_back("apss");
	all_algorithms.push_back("fourier");
	all_algorithms.push_back("imls");
	all_algorithms.push_back("mpu");
	//all_algorithms.push_back("mpusmooth");
	all_algorithms.push_back("poisson");
	all_algorithms.push_back("rbf");
	all_algorithms.push_back("scattered");
	all_algorithms.push_back("spss");
	all_algorithms.push_back("wavelet");

	//double x_start = 5;
	//double x_inc = 1;
	//double x_start = 40;
	//double x_inc = (90.0-x_start)/(double)(num_shapes-1);
	double x_start = 5.0;
	double x_end = 30.0;
	double x_inc = (x_end-x_start)/(double)(num_shapes-1);

	for(int i = 0; i < num_shapes; i++)  {
		if(i == 0)  {
			x_start += x_inc;
			continue;
		}

		fprintf(distmean_file, "%.10f ", x_start);
		fprintf(distmax_file, "%.10f ", x_start);
		fprintf(anglemean_file, "%.10f ", x_start);
		fprintf(anglemax_file, "%.10f ", x_start);

		for(unsigned a = 0; a < all_algorithms.size(); a++)  {
			string algorithm = all_algorithms[a];

			char* eval_num = new char[30];
			sprintf(eval_num, "%u", i);
			string base_alg_file = eval_base+"/"+shape+"/"+exp+"/"+algorithm+"/"+shape+"_"+eval_num;
			string stats_filename = base_alg_file+".recon";
			string dist_filename = base_alg_file+".dist";
			FILE* stats_file = fopen(stats_filename.c_str(), "rb");
			FILE* dist_file = fopen(dist_filename.c_str(), "rb");

			if(stats_file == 0 || dist_file == 0)  {
				cerr << "bad stat file or dist file! " << base_alg_file << endl;
				continue;
			}

			size_t num_read = 0;
			double shape_dist[8];
			num_read = fread(shape_dist, sizeof(double), 8, dist_file);
			double angle_dist[8];
			num_read = fread(angle_dist, sizeof(double), 8, dist_file);
			fclose(dist_file);
			fclose(stats_file);

			double mean_dist = shape_dist[5] > shape_dist[6] ? shape_dist[5] : shape_dist[6];
			//double mean_dist = shape_dist[7];
			double max_dist = shape_dist[4];
			double mean_angle = angle_dist[5] > angle_dist[6] ? angle_dist[5] : angle_dist[6];
			//double mean_angle = angle_dist[7];
			double max_angle = angle_dist[4];

			if(a == (all_algorithms.size()-1))  {
				fprintf(distmean_file, "%.10f\n", mean_dist);
				fprintf(distmax_file, "%.10f\n", max_dist);
				fprintf(anglemean_file, "%.10f\n", mean_angle);
				fprintf(anglemax_file, "%.10f\n", max_angle);
			}
			else  {
				fprintf(distmean_file, "%.10f ", mean_dist);
				fprintf(distmax_file, "%.10f ", max_dist);
				fprintf(anglemean_file, "%.10f ", mean_angle);
				fprintf(anglemax_file, "%.10f ", max_angle);
			}

			delete [] eval_num;
		}

		x_start += x_inc;
	}

	fclose(distmean_file);
	fclose(distmax_file);
	fclose(anglemean_file);
	fclose(anglemax_file);
}
