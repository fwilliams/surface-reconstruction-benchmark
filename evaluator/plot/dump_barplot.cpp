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

void dump_init(FILE* _file, int _numShapes, string _title)  {
	fprintf(_file, "set logscale y\n");
	fprintf(_file, "set xtics\n");
	fprintf(_file, "set title \"%s\" font \"Helvetica,35\"\n", _title.c_str());
	fprintf(_file, "set size 1.15,1\n");
	fprintf(_file, "set origin 0,0\n");
	fprintf(_file, "set boxwidth 0.75\n");

	fprintf(_file, "set style line 101 lc rgbcolor \"#7FBF67\" lw 8 lt 1\n");
	fprintf(_file, "set style line 102 lc rgbcolor \"#D37281\" lw 8 lt 1\n");
	fprintf(_file, "set style line 103 lc rgbcolor \"#4E7B8C\" lw 8 lt 1\n");
	fprintf(_file, "set style line 104 lc rgbcolor \"#E0CA79\" lw 8 lt 1\n");
	fprintf(_file, "set style line 105 lc rgbcolor \"#655D98\" lw 8 lt 1\n");
	fprintf(_file, "set style line 106 lc rgbcolor \"#7FBF67\" lw 8 lt 2\n");
	fprintf(_file, "set style line 107 lc rgbcolor \"#D37281\" lw 8 lt 2\n");
	fprintf(_file, "set style line 108 lc rgbcolor \"#4E7B8C\" lw 8 lt 2\n");
	fprintf(_file, "set style line 109 lc rgbcolor \"#E0CA79\" lw 8 lt 2\n");
	fprintf(_file, "set style line 110 lc rgbcolor \"#655D98\" lw 8 lt 2\n");

	fprintf(_file, "set term postscript eps enhanced color dashlength 2 blacktext \"Helvetica\" 20\n");

	fprintf(_file, "set key off\n");
	fprintf(_file, "set xrange [-0.55:%.7f]\n", (_numShapes-0.45));
}

void dump_dist_plot(FILE* _file, string _dataFile, string _epsFile, vector<double> _maxDists, vector<double> _meanDists, vector<int> _styles, vector<string> _titles)  {
	double y_min = 1e10, y_max = -1e10;
	int n_shapes = _maxDists.size();
	for(int i = 0; i < n_shapes; i++)  {
		y_min = _maxDists[i] < y_min ? _maxDists[i] : y_min;
		y_min = _meanDists[i] < y_min ? _meanDists[i] : y_min;
		y_max = _maxDists[i] > y_max ? _maxDists[i] : y_max;
		y_max = _meanDists[i] > y_max ? _meanDists[i] : y_max;
	}

	double plot_y_min = 0.8*y_min;
	double plot_y_max = 1.05*y_max;

	fprintf(_file, "set output \'%s\'\n", _epsFile.c_str());

	fprintf(_file, "set yrange [%.7f:%.7f]\n", plot_y_min, plot_y_max);
	fprintf(_file, "plot \\\n");

	for(int i = 0; i < n_shapes; i++)  {
		string next_title = _titles[i];
		if(i == 0)
			fprintf(_file, " \'%s\' ", _dataFile.c_str());
		else
			fprintf(_file, " \'\' ");

		fprintf(_file, "using (%u):%u every ::1::1 with boxes fs solid 0.3 ls %u notitle , \\\n",
				  i, (i+1), _styles[i]);
		fprintf(_file, " \'\' ");
		fprintf(_file, "using (%u):(%u):xticlabel(\"%s\") with points pointsize 0 notitle , \\\n",
				i, i, next_title.c_str());
		fprintf(_file, " \'\' ");
		fprintf(_file, "using (%u):%u:%u:%u:%u every ::0::0 with candlesticks lt -1 lw 4 notitle",
				i, (i+1), (i+1), (i+1), (i+1), i, i);

		if(i == (n_shapes-1))
			fprintf(_file, "\n");
		else
			fprintf(_file, " , \\\n");
	}

	fprintf(_file, "set output\n");
}

void dump_angle_plot(FILE* _file, string _dataFile, string _epsFile, vector<double> _meanDists, vector<int> _styles, vector<string> _titles)  {
	double y_min = 1e10, y_max = -1e10;
	int n_shapes = _meanDists.size();
	for(int i = 0; i < n_shapes; i++)  {
		y_min = _meanDists[i] < y_min ? _meanDists[i] : y_min;
		y_max = _meanDists[i] > y_max ? _meanDists[i] : y_max;
	}

	double plot_y_min = 0.8*y_min;
	double plot_y_max = 1.05*y_max;

	fprintf(_file, "set output \'%s\'\n", _epsFile.c_str());

	fprintf(_file, "set yrange [%.7f:%.7f]\n", plot_y_min, plot_y_max);
	fprintf(_file, "plot \\\n");

	for(int i = 0; i < n_shapes; i++)  {
		string next_title = _titles[i];
		if(i == 0)
			fprintf(_file, " \'%s\' ", _dataFile.c_str());
		else
			fprintf(_file, " \'\' ");

		fprintf(_file, "using (%u):%u every ::0::0 with boxes fs solid 0.3 ls %u notitle , \\\n",
				  i, (i+1), _styles[i]);
		fprintf(_file, " \'\' ");
		fprintf(_file, "using (%u):(%u):xticlabel(\"%s\") with points pointsize 0 notitle",
				i, i, next_title.c_str());

		if(i == (n_shapes-1))
			fprintf(_file, "\n");
		else
			fprintf(_file, " , \\\n");
	}

	fprintf(_file, "set output\n");
}
int main(int argc, char** argv)  {
	if(argc != 5)  {
		std::cerr << "usage: " << argv[0] << " evaluation_base shape pc_num data_base" << std::endl;
		return 1;
	}
	int arg_num = 1;
	string eval_base = argv[arg_num++];
	string shape = argv[arg_num++];
	int pc_num = atoi(argv[arg_num++]);

	string data_base = argv[arg_num++];
	string dist_filename = data_base + "_bardist.gnu";
	string dist_dataname = data_base + "_bardist.dat";
	string dist_eps = data_base + "_dist.eps";
	string angle_filename = data_base + "_barangle.gnu";
	string angle_dataname = data_base + "_barangle.dat";
	string angle_eps = data_base + "_angle.eps";

	FILE* dist_file = fopen(dist_filename.c_str(), "w");
	FILE* dist_data_file = fopen(dist_dataname.c_str(), "w");
	FILE* angle_file = fopen(angle_filename.c_str(), "w");
	FILE* angle_data_file = fopen(angle_dataname.c_str(), "w");

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

	vector<string> algorithm_titles;
	algorithm_titles.push_back("APSS");
	algorithm_titles.push_back("Fouri");
	algorithm_titles.push_back("IMLS");
	algorithm_titles.push_back("MPU");
	algorithm_titles.push_back("MPUS");
	algorithm_titles.push_back("Pois");
	algorithm_titles.push_back("RBF");
	algorithm_titles.push_back("Scat");
	algorithm_titles.push_back("SPSS");
	algorithm_titles.push_back("Wave");

	vector<double> mean_dists;
	vector<double> max_dists;
	vector<double> mean_angles;

	vector<string> titles;
	vector<int> styles;
	int line_styles = 101;

	for(unsigned a = 0; a < all_algorithms.size(); a++)  {
		string algorithm = all_algorithms[a];

		char* eval_num = new char[30];
		sprintf(eval_num, "%u", pc_num);
		string base_alg_file = eval_base+"/"+algorithm+"/"+shape+"_"+eval_num;
		string dist_filename = base_alg_file+".dist";
		FILE* dist_file = fopen(dist_filename.c_str(), "rb");

		if(dist_file == 0)  {
			cerr << "bad dist file! " << base_alg_file << endl;
			line_styles++;
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

		mean_dists.push_back(pc_mean_dist);
		max_dists.push_back(pc_max_dist);
		mean_angles.push_back(pc_mean_angle);

		delete [] eval_num;

		titles.push_back(algorithm_titles[a]);
		styles.push_back(line_styles);

		line_styles++;
	}

	for(unsigned i = 0; i < mean_dists.size(); i++)  {
		if(i != (mean_dists.size()-1))
			fprintf(dist_data_file, "%.7f ", mean_dists[i]);
		else
			fprintf(dist_data_file, "%.7f\n", mean_dists[i]);
	}
	for(unsigned i = 0; i < max_dists.size(); i++)  {
		if(i != (max_dists.size()-1))
			fprintf(dist_data_file, "%.7f ", max_dists[i]);
		else
			fprintf(dist_data_file, "%.7f", max_dists[i]);
	}

	for(unsigned i = 0; i < mean_angles.size(); i++)  {
		if(i != (mean_angles.size()-1))
			fprintf(angle_data_file, "%.7f ", mean_angles[i]);
		else
			fprintf(angle_data_file, "%.7f", mean_angles[i]);
	}

	dump_init(dist_file, titles.size(), "Distance (mm)");
	dump_dist_plot(dist_file, dist_dataname, dist_eps, max_dists, mean_dists, styles, titles);
	dump_init(angle_file, titles.size(), "Angle ({/Symbol \\260})");
	dump_angle_plot(angle_file, angle_dataname, angle_eps, mean_angles, styles, titles);

	fclose(dist_file);
	fclose(angle_file);
	fclose(dist_data_file);
	fclose(angle_data_file);
}
