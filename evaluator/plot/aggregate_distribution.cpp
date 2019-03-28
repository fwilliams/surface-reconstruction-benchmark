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

void dump_init(FILE* _file, string _xlabel, string _ylabel)  {
	fprintf(_file, "set xrange [0:1]\n");
	fprintf(_file, "set boxwidth 0.35\n");
	fprintf(_file, "set style fill empty\n");
	fprintf(_file, "set key top rm vertical\n");
	fprintf(_file, "set xtics nomirror\n");
	fprintf(_file, "set size 1,2\n");
	//fprintf(_file, "set logscale y\n");
	fprintf(_file, "set xlabel \"%s\" font \"Helvetica,50\"\n", _xlabel.c_str());
	fprintf(_file, "set ylabel \"%s\" font \"Helvetica,50\"\n", _ylabel.c_str());

	fprintf(_file, "unset xtics\n");
	fprintf(_file, "unset xrange\n");

	fprintf(_file, "set style line 401 lc rgbcolor \"#E0CA79\" lw 16 lt 1\n");
	fprintf(_file, "set term postscript eps enhanced color dashlength 2.5 blacktext \"Helvetica\" 33\n");
	fprintf(_file, "set style rect fc rgb \"black\" fs solid 0.07 border\n");
}

void dump_distribution(FILE* _file, string _dataFile, string _epsFile, double _yMin, double _yMax)  {
	double plot_y_min = 0.8*_yMin;
	double plot_y_max = 1.05*_yMax;

	fprintf(_file, "set yrange [%.7f:%.7f]\n", plot_y_min, plot_y_max);

	fprintf(_file, "set output \'%s\'\n", _epsFile.c_str());

	fprintf(_file, "plot \\\n");
	fprintf(_file, "\'%s\' using (0.5):2:1:5:4 every ::0::0 with candlesticks ls 401 fs solid 0.3 notitle whiskerbars 0.9, \\\n", _dataFile.c_str());
	fprintf(_file, "\'\' using (0.5):3:3:3:3 every ::0::0 with candlesticks lt -1 lw 5 notitle\n");
	fprintf(_file, "set output\n");
	fprintf(_file, "!epstopdf %s\n", _epsFile.c_str());
}

int main(int argc, char** argv)  {
	if(argc != 4)  {
		std::cerr << "usage: " << argv[0] << " dist_base num_pcs data_base" << std::endl;
		return 1;
	}
	int arg_num = 1;
	string input_base = argv[arg_num++];
	int num_pcs = atoi(argv[arg_num++]);
	string data_base = argv[arg_num++];

	string distmean_filename = data_base+"_distmean.dat";
	string distmax_filename = data_base+"_distmax.dat";
	string anglemean_filename = data_base+"_anglemean.dat";
	string anglemax_filename = data_base+"_anglemax.dat";

	FILE* distmean_file = fopen(distmean_filename.c_str(), "w");
	FILE* distmax_file = fopen(distmax_filename.c_str(), "w");
	FILE* anglemean_file = fopen(anglemean_filename.c_str(), "w");
	FILE* anglemax_file = fopen(anglemax_filename.c_str(), "w");

	vector<FILE*> all_dat_files;
	all_dat_files.push_back(distmean_file);
	all_dat_files.push_back(distmax_file);
	all_dat_files.push_back(anglemean_file);
	all_dat_files.push_back(anglemax_file);

	vector<double> mean_dist;
	vector<double> max_dist;
	vector<double> mean_angle;
	vector<double> max_angle;

	for(int i = 0; i < num_pcs; i++)  {
		char* eval_num = new char[30];
		sprintf(eval_num, "%u", i);
		string dist_filename = input_base+"_"+eval_num+".dist";
		FILE* dist_file = fopen(dist_filename.c_str(), "rb");

		if(dist_file == 0)  {
			cerr << "bad dist file! " << dist_filename << endl;
			continue;
		}

		size_t num_read = 0;
		double shape_dist[8];
		num_read = fread(shape_dist, sizeof(double), 8, dist_file);
		double angle_dist[8];
		num_read = fread(angle_dist, sizeof(double), 8, dist_file);
		fclose(dist_file);

		double pc_mean_dist = shape_dist[5] > shape_dist[6] ? shape_dist[5] : shape_dist[6];
		double pc_max_dist = shape_dist[4];
		double pc_mean_angle = angle_dist[5] > angle_dist[6] ? angle_dist[5] : angle_dist[6];
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

	vector<double> min_values;
	vector<double> max_values;

	for(int i = 0; i < distributions.size(); i++)  {
		vector<double> distribution = distributions[i];
		vector<Measure> full_distribution;
		for(unsigned d = 0; d < distribution.size(); d++)
			full_distribution.push_back(Measure(d, distribution[d]));
		FILE* dist_file = all_dat_files[i];

		std::sort(full_distribution.begin(),full_distribution.end());
		double quartile_info[5];
		int quartile_indices[5];
		quartiles(&full_distribution, quartile_info, quartile_indices);
		for(int i = 0; i < 5; i++)
			fprintf(dist_file, "%.10f ", quartile_info[i]);
		for(int i = 0; i < 4; i++)
			fprintf(dist_file, "%u ", quartile_indices[i]);
		fprintf(dist_file, "%u\n", quartile_indices[4]);

		min_values.push_back(quartile_info[0]);
		max_values.push_back(quartile_info[4]);
	}
	fclose(distmean_file);
	fclose(distmax_file);
	fclose(anglemean_file);
	fclose(anglemax_file);

	string distmean_gnu_filename = data_base+"_distmean.gnu";
	string distmean_eps_filename = data_base+"_distmean.eps";
	FILE* distmean_gnu_file = fopen(distmean_gnu_filename.c_str(), "w");
	string distmax_gnu_filename = data_base+"_distmax.gnu";
	string distmax_eps_filename = data_base+"_distmax.eps";
	FILE* distmax_gnu_file = fopen(distmax_gnu_filename.c_str(), "w");

	string anglemean_gnu_filename = data_base+"_anglemean.gnu";
	string anglemean_eps_filename = data_base+"_anglemean.eps";
	FILE* anglemean_gnu_file = fopen(anglemean_gnu_filename.c_str(), "w");
	string anglemax_gnu_filename = data_base+"_anglemax.gnu";
	string anglemax_eps_filename = data_base+"_anglemax.eps";
	FILE* anglemax_gnu_file = fopen(anglemax_gnu_filename.c_str(), "w");

	dump_init(distmean_gnu_file, "Position", "Mean Distance (mm)");
	dump_distribution(distmean_gnu_file, distmean_filename, distmean_eps_filename, min_values[0], max_values[0]);

	dump_init(distmax_gnu_file, "Position", "Hausdorff Distance (mm)");
	dump_distribution(distmax_gnu_file, distmax_filename, distmax_eps_filename, min_values[1], max_values[1]);

	dump_init(anglemean_gnu_file, "Normal", "Mean Angle Deviation ({/Symbol \\260})");
	dump_distribution(anglemean_gnu_file, anglemean_filename, anglemean_eps_filename, min_values[2], max_values[2]);

	dump_init(anglemax_gnu_file, "Normal", "Max Angle Deviation ({/Symbol \\260})");
	dump_distribution(anglemax_gnu_file, anglemax_filename, anglemax_eps_filename, min_values[3], max_values[3]);

	fclose(distmean_gnu_file);
	fclose(distmax_gnu_file);
	fclose(anglemean_gnu_file);
	fclose(anglemax_gnu_file);

	char exec_gnu_distmean[300];
	sprintf(exec_gnu_distmean, "gnuplot %s", distmean_gnu_filename.c_str());
	char exec_gnu_distmax[300];
	sprintf(exec_gnu_distmax, "gnuplot %s", distmax_gnu_filename.c_str());
	char exec_gnu_anglemean[300];
	sprintf(exec_gnu_anglemean, "gnuplot %s", anglemean_gnu_filename.c_str());
	char exec_gnu_anglemax[300];
	sprintf(exec_gnu_anglemax, "gnuplot %s", anglemax_gnu_filename.c_str());

	int result = system(exec_gnu_distmean);
	result = system(exec_gnu_distmax);
	result = system(exec_gnu_anglemean);
	result = system(exec_gnu_anglemax);

	char rm_exec[500];
	sprintf(rm_exec, "rm %s %s %s %s %s %s %s %s %s %s %s %s", distmean_filename.c_str(), distmax_filename.c_str(), anglemean_filename.c_str(), anglemax_filename.c_str(), distmean_gnu_filename.c_str(), distmax_gnu_filename.c_str(), anglemean_gnu_filename.c_str(), anglemax_gnu_filename.c_str(), distmean_eps_filename.c_str(), distmax_eps_filename.c_str(), anglemean_eps_filename.c_str(), anglemax_eps_filename.c_str());
	result = system(rm_exec);
}
