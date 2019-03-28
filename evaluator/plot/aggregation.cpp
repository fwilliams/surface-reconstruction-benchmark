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

struct Measure  {
	Measure()  {}
	Measure(int _index, double _measure) : i(_index), m(_measure)  {}
	inline bool operator<(const Measure& _otherMeasure) const {
		return m < _otherMeasure.m;
	}
	int i;
	double m;
};

void quartiles(vector<Measure>* _distribution, double* stats, int* stat_inds)  {
	int size = _distribution->size();
	stats[2] = 0.5*(_distribution->at((size-1)/2).m+_distribution->at(size/2).m);
	stat_inds[2] = _distribution->at((size-1)/2).i;

	double lq_perc = 0.25*(size-1);
	int lq_low = lq_perc;
	int lq_high = lq_low+1;
	stats[1] = (1.0-(lq_perc-lq_low))*_distribution->at(lq_low).m + (lq_perc-lq_low)*_distribution->at(lq_high).m;
	stat_inds[1] = _distribution->at(lq_low).i;

	double uq_perc = 0.75*(size-1);
	int uq_low = uq_perc;
	int uq_high = uq_low+1;
	stats[3] = (1.0-(uq_perc-uq_low))*_distribution->at(uq_low).m + (uq_perc-uq_low)*_distribution->at(uq_high).m;
	stat_inds[3] = _distribution->at(uq_low).i;

	stats[0] = _distribution->at(0).m;
	stat_inds[0] = _distribution->at(0).i;
	stats[4] = _distribution->at(size-1).m;
	stat_inds[4] = _distribution->at(size-1).i;
}

int main(int argc, char** argv)  {
	if(argc < 5)  {
		std::cerr << "usage: " << argv[0] << " evaluation_base shape num_shapes data_base [exp]" << std::endl;
		return 1;
	}
	int arg_num = 1;
	string eval_base = argv[arg_num++];
	string shape = argv[arg_num++];
	int num_shapes = atoi(argv[arg_num++]);

	string data_base = argv[arg_num++];
	string distmean_filename = data_base + "_distmean.txt";
	string distmax_filename = data_base + "_distmax.txt";
	string anglemean_filename = data_base + "_anglemean.txt";
	string anglemax_filename = data_base + "_anglemax.txt";

	string exp = "";
	if(argc == 6)
		exp = argv[5];

	FILE* distmean_file = fopen(distmean_filename.c_str(), "w");
	FILE* distmax_file = fopen(distmax_filename.c_str(), "w");
	FILE* anglemean_file = fopen(anglemean_filename.c_str(), "w");
	FILE* anglemax_file = fopen(anglemax_filename.c_str(), "w");

	vector<FILE*> all_files;
	all_files.push_back(distmean_file);
	all_files.push_back(distmax_file);
	all_files.push_back(anglemean_file);
	all_files.push_back(anglemax_file);

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

	for(unsigned a = 0; a < all_algorithms.size(); a++)  {
		string algorithm = all_algorithms[a];

		vector<double> mean_dist;
		vector<double> max_dist;
		vector<double> mean_angle;
		vector<double> max_angle;

		vector<double> computation_time;

		for(int i = 0; i < num_shapes; i++)  {
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

			double pc_mean_dist = shape_dist[5] > shape_dist[6] ? shape_dist[5] : shape_dist[6];
			//double pc_mean_dist = shape_dist[7];
			double pc_max_dist = shape_dist[4];
			double pc_mean_angle = angle_dist[5] > angle_dist[6] ? angle_dist[5] : angle_dist[6];
			//double pc_mean_angle = angle_dist[7];
			double pc_max_angle = angle_dist[4];

			mean_dist.push_back(pc_mean_dist);
			max_dist.push_back(pc_max_dist);
			mean_angle.push_back(pc_mean_angle);
			max_angle.push_back(pc_max_angle);

			delete [] eval_num;
		}

		vector< vector<double> > distributions;
		distributions.push_back(mean_dist);
		distributions.push_back(max_dist);
		distributions.push_back(mean_angle);
		distributions.push_back(max_angle);

		for(int i = 0; i < 4; i++)  {
			vector<double> distribution = distributions[i];
			vector<Measure> full_distribution;
			for(unsigned d = 0; d < distribution.size(); d++)
				full_distribution.push_back(Measure(d, distribution[d]));
			FILE* dist_file = all_files[i];

			std::sort(full_distribution.begin(),full_distribution.end());
			double quartile_info[5];
			int quartile_indices[5];
			quartiles(&full_distribution, quartile_info, quartile_indices);
			for(int i = 0; i < 5; i++)
				fprintf(dist_file, "%.10f ", quartile_info[i]);
			for(int i = 0; i < 4; i++)
				fprintf(dist_file, "%u ", quartile_indices[i]);
			fprintf(dist_file, "%u\n", quartile_indices[4]);
		}
	}

	fclose(distmean_file);
	fclose(distmax_file);
	fclose(anglemean_file);
	fclose(anglemax_file);
}
