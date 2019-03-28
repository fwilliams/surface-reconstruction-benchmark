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

void dump_dist_plot(FILE* _file, string _dataFile, string _epsFile, double _yMin, double _yMax)  {
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

void dump_angle_plot(FILE* _file, string _dataFile, string _epsFile, double _yMin, double _yMax)  {
	double plot_y_min = 0.8*_yMin;
	double plot_y_max = 1.05*_yMax;

	fprintf(_file, "set yrange [%.7f:%.7f]\n", plot_y_min, plot_y_max);

	fprintf(_file, "set output \'%s\'\n", _epsFile.c_str());

	fprintf(_file, "plot \\\n");
	fprintf(_file, "\'%s\' using (0.5):7:6:10:9 every ::0::0 with candlesticks ls 401 fs solid 0.3 notitle whiskerbars 0.9, \\\n", _dataFile.c_str());
	fprintf(_file, "\'\' using (0.5):8:8:8:8 every ::0::0 with candlesticks lt -1 lw 5 notitle\n");
	fprintf(_file, "set output\n");
	fprintf(_file, "!epstopdf %s\n", _epsFile.c_str());
}

int main(int argc, char** argv)  {
	if(argc != 3)  {
		std::cerr << "usage: " << argv[0] << " dist_file data_base" << std::endl;
		return 1;
	}
	int arg_num = 1;
	string input_filename = argv[arg_num++];

	string data_base = argv[arg_num++];
	string data_filename = data_base + "_data.dat";
	string dist_filename = data_base + "_dist.gnu";
	string angle_filename = data_base + "_angle.gnu";

	string dist_eps = data_base + "_dist.eps";
	string angle_eps = data_base + "_angle.eps";

	FILE* data_file = fopen(data_filename.c_str(), "w");
	FILE* dist_file = fopen(dist_filename.c_str(), "w");
	FILE* angle_file = fopen(angle_filename.c_str(), "w");

	FILE* input_file = fopen(input_filename.c_str(), "rb");
	if(input_file == 0)  {
		cerr << "bad dist file! " << input_filename << endl;
		return -1;
	}

	size_t num_read = 0;
	double shape_dist[8];
	num_read = fread(shape_dist, sizeof(double), 8, input_file);
	double angle_dist[8];
	num_read = fread(angle_dist, sizeof(double), 8, input_file);
	fclose(input_file);

	fprintf(data_file, "%.7f %.7f %.7f %.7f %.7f %.7f %.7f %.7f %.7f %.7f\n", shape_dist[0], shape_dist[1], shape_dist[2], shape_dist[3], shape_dist[4],
																									  angle_dist[0], angle_dist[1], angle_dist[2], angle_dist[3], angle_dist[4]);
	fclose(data_file);

	dump_init(dist_file, "Position", "error (mm)");
	dump_dist_plot(dist_file, data_filename, dist_eps, shape_dist[0], shape_dist[4]);
	dump_init(angle_file, "Normal", "error ({/Symbol \\260})");
	dump_angle_plot(angle_file, data_filename, angle_eps, angle_dist[0], angle_dist[4]);

	fclose(dist_file);
	fclose(angle_file);

	char exec_gnu_dist[300];
	sprintf(exec_gnu_dist, "gnuplot %s", dist_filename.c_str());
	char exec_gnu_angle[300];
	sprintf(exec_gnu_angle, "gnuplot %s", angle_filename.c_str());
	int result = system(exec_gnu_dist);
	result = system(exec_gnu_angle);

	char rm_exec[500];
	sprintf(rm_exec, "rm %s %s %s %s %s", data_filename.c_str(), dist_filename.c_str(), angle_filename.c_str(), dist_eps.c_str(), angle_eps.c_str());
	result = system(rm_exec);
}
