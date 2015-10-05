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

#include "UniformSampler.h"
#include "SphericalImager.h"
#include "../modeling/shape_loader.h"

#include <cstdio>
#include<cstdlib>
#include <iostream>
using namespace std;

#define NUM_ARGS 13

void update_indices(int* indices, int* index_sizes)  {
	int max_index = -1;
	for(int i = (NUM_ARGS-1); i >= 0; i--)  {
		if(index_sizes[i] == 0)
			continue;
		if((indices[i]+1) == index_sizes[i])
			continue;

		indices[i] = indices[i]+1;
		max_index = i;
		break;
	}

	if(max_index == -1)
		return;

	for(int i = (max_index+1); i < NUM_ARGS; i++)
		indices[i] = 0;
}

bool configs_remain(int* indices, int* index_sizes)  {
	int total_configs = 0, current_configs = 0;
	for(int i = 0; i < NUM_ARGS; i++)  {
		if(index_sizes[i] > 0)
			current_configs += (indices[i]+1);
		total_configs += index_sizes[i];
	}
	return current_configs < total_configs;
}

int my_round(double _val)  {
	int rounded_val = (int)(_val+0.5);
	return rounded_val;
}

vector<int> extract_int_params(int _argNum, char** argv, int& arg_num)  {
	int arg_place = _argNum;
	vector<int> params;
	string arg_string = argv[arg_place++];

	if(arg_string.compare("range") == 0)  {
		int min = atoi(argv[arg_place++]);
		int max = atoi(argv[arg_place++]);
		int num = atoi(argv[arg_place++]);

		for(int i = 0; i < num; i++)  {
			double i_perc = (double)i / (double)(num-1);
			int next_val = my_round(min*(1.0-i_perc) + max*i_perc);
			if(params.size() > 0 && params[params.size()-1] == next_val)
				continue;
			params.push_back(next_val);
		}
	}
	else
		params.push_back(atoi(arg_string.c_str()));

	arg_num = arg_place;
	return params;
}

vector<double> extract_double_params(int _argNum, char** argv, int& arg_num)  {
	int arg_place = _argNum;
	vector<double> params;
	string arg_string = argv[arg_place++];

	if(arg_string.compare("range") == 0)  {
		double min = atof(argv[arg_place++]);
		double max = atof(argv[arg_place++]);
		int num = atoi(argv[arg_place++]);

		for(int i = 0; i < num; i++)  {
			double i_perc = (double)i / (double)(num-1);
			double next_val = min*(1.0-i_perc) + max*i_perc;
			if(params.size() > 0 && params[params.size()-1] == next_val)
				continue;
			params.push_back(next_val);
		}
	}
	else
		params.push_back(atof(arg_string.c_str()));

	arg_num = arg_place;
	return params;
}

int main(int argc, char** argv)  {
	bool dump_images = false;
	if(dump_images)  {
		if(argc != 6)  {
			cerr << "usage: " << argv[0] << " shape res num_scans min_range base_file" << endl;
			return 1;
		}

		ImplicitFunction* shape = ShapeLoader::load_shape(argv[1]);
		int res = atoi(argv[2]);
		int num_scans = atoi(argv[3]);
		double min_range = atof(argv[4]);
		string base_file = argv[5];

		SphericalImager imager(shape, num_scans, res, res);
		imager.setMinRange(min_range);
		imager.setRandomSampleRotation(false);
		imager.sample(base_file);

		return 0;
	}

	if(argc < 3)  {
		cerr << "usage: " << argv[0] << " shape_file config_base ([camera_param range min max num])* ([camera_param value])*" << endl;
		return 1;
	}

	int arg_num = 1;
	string shape_filename = argv[arg_num++];
	string config_base = argv[arg_num++];

	vector<int> resolutions;
	vector<int> scans;

	vector<double> min_ranges;
	vector<double> max_ranges;
	vector<int> num_stripeses;
	vector<double> laser_fovs;
	vector<double> peak_thresholds;
	vector<double> std_thresholds;
	vector<double> additive_noises;
	vector<double> laser_smoothers;

	vector<double> registration_noises;

	vector<int> to_random_sample_rotation;
	vector<int> normal_types;

	while(arg_num < argc)  {
		string next_arg = argv[arg_num++];
		int temp_arg = arg_num;

		if(next_arg.compare("res") == 0)  {
			vector<int> params = extract_int_params(temp_arg, argv, arg_num);
			for(unsigned i = 0; i < params.size(); i++)
				resolutions.push_back(params[i]);
		}
		else if(next_arg.compare("scans") == 0)  {
			vector<int> params = extract_int_params(temp_arg, argv, arg_num);
			for(unsigned i = 0; i < params.size(); i++)
				scans.push_back(params[i]);
		}

		// scanner
		else if(next_arg.compare("min_range") == 0)  {
			vector<double> params = extract_double_params(temp_arg, argv, arg_num);
			for(unsigned i = 0; i < params.size(); i++)
				min_ranges.push_back(params[i]);
		}
		else if(next_arg.compare("max_range") == 0)  {
			vector<double> params = extract_double_params(temp_arg, argv, arg_num);
			for(unsigned i = 0; i < params.size(); i++)
				max_ranges.push_back(params[i]);
		}

		else if(next_arg.compare("num_stripes") == 0)  {
			vector<int> params = extract_int_params(temp_arg, argv, arg_num);
			for(unsigned i = 0; i < params.size(); i++)
				num_stripeses.push_back(params[i]);
		}
		else if(next_arg.compare("laser_fov") == 0)  {
			vector<double> params = extract_double_params(temp_arg, argv, arg_num);
			for(unsigned i = 0; i < params.size(); i++)
				laser_fovs.push_back(params[i]);
		}

		else if(next_arg.compare("peak_threshold") == 0)  {
			vector<double> params = extract_double_params(temp_arg, argv, arg_num);
			for(unsigned i = 0; i < params.size(); i++)
				peak_thresholds.push_back(params[i]);
		}
		else if(next_arg.compare("std_threshold") == 0)  {
			vector<double> params = extract_double_params(temp_arg, argv, arg_num);
			for(unsigned i = 0; i < params.size(); i++)
				std_thresholds.push_back(params[i]);
		}

		else if(next_arg.compare("additive_noise") == 0)  {
			vector<double> params = extract_double_params(temp_arg, argv, arg_num);
			for(unsigned i = 0; i < params.size(); i++)
				additive_noises.push_back(params[i]);
		}
		else if(next_arg.compare("laser_smoother") == 0)  {
			vector<double> params = extract_double_params(temp_arg, argv, arg_num);
			for(unsigned i = 0; i < params.size(); i++)
				laser_smoothers.push_back(params[i]);
		}

		// registration
		else if(next_arg.compare("registration_noise") == 0)  {
			vector<double> params = extract_double_params(temp_arg, argv, arg_num);
			for(unsigned i = 0; i < params.size(); i++)
				registration_noises.push_back(params[i]);
		}

		// misc
		else if(next_arg.compare("random_rotation") == 0)  {
			vector<int> params = extract_int_params(temp_arg, argv, arg_num);
			for(unsigned i = 0; i < params.size(); i++)
				to_random_sample_rotation.push_back(params[i]);
		}

		else if(next_arg.compare("normal_type") == 0)  {
			vector<int> params = extract_int_params(temp_arg, argv, arg_num);
			for(unsigned i = 0; i < params.size(); i++)
				normal_types.push_back(params[i]);
		}

		else  {
			cerr << "unknown option: " << next_arg << endl;
		}
	}

	if(scans.size() == 0)  {
		cerr << "need to specify some number of scans!" << endl;
		return 1;
	}

	if(resolutions.size() == 0)  {
		cerr << "need to specify some resolutions!" << endl;
		return 1;
	}

	int* indices = new int[NUM_ARGS];
	int index_sizes[] = { resolutions.size() , scans.size() , /* necessary stuff */
								 min_ranges.size() , max_ranges.size() , num_stripeses.size() , laser_fovs.size() ,
								 peak_thresholds.size() , std_thresholds.size() , additive_noises.size() , laser_smoothers.size() , /* scanner */
								 registration_noises.size() , to_random_sample_rotation.size() , normal_types.size() /* misc */
							  };
	for(int i = 0; i < NUM_ARGS; i++)
		indices[i] = 0;

	int file_ind = 0;
	while(true)  {
		ostringstream file_ind_str;
		file_ind_str << file_ind;
		string config_filename = config_base + "_" + file_ind_str.str() + ".cnf";
		string npts_filename = config_base + "_" + file_ind_str.str() + ".npts";
		FILE* config_file = fopen(config_filename.c_str(), "w");

		fprintf(config_file, "[sampler]\n");
		fprintf(config_file, "pathdir:bin\n");
		fprintf(config_file, "sampling_type:uniform\n");
		fprintf(config_file, "infile:%s\n", shape_filename.c_str());
		fprintf(config_file, "outfile:%s\n", npts_filename.c_str());

		// necessities
		fprintf(config_file, "[uniform]\n");
		fprintf(config_file, "exec_name:uniform_sampler\n");
		fprintf(config_file, "camera_res_x:%u\n", resolutions[ indices[0] ]);
		fprintf(config_file, "camera_res_y:%u\n", resolutions[ indices[0] ]);
		fprintf(config_file, "scan_res:%u\n", scans[ indices[1] ]);

		// scanner
		if(min_ranges.size() > 0)
			fprintf(config_file, "min_range:%.7f\n", min_ranges[ indices[2] ]);
		if(max_ranges.size() > 0)
			fprintf(config_file, "max_range:%.7f\n", max_ranges[ indices[3] ]);

		if(num_stripeses.size() > 0)
			fprintf(config_file, "num_stripes:%u\n", num_stripeses[ indices[4] ]);
		if(laser_fovs.size() > 0)
			fprintf(config_file, "laser_fov:%.7f\n", laser_fovs[ indices[5] ]);

		if(peak_thresholds.size() > 0)
			fprintf(config_file, "peak_threshold:%.7f\n", peak_thresholds[ indices[6] ]);
		if(std_thresholds.size() > 0)
			fprintf(config_file, "std_threshold:%.7f\n", std_thresholds[ indices[7] ]);

		if(additive_noises.size() > 0)
			fprintf(config_file, "additive_noise:%.7f\n", additive_noises[ indices[8] ]);
		if(laser_smoothers.size() > 0)
			fprintf(config_file, "laser_smoother:%.7f\n", laser_smoothers[ indices[9] ]);

		// misc
		if(registration_noises.size() > 0)
			fprintf(config_file, "registration_error:%.7f\n", registration_noises[ indices[10] ]);
		if(to_random_sample_rotation.size() > 0)
			fprintf(config_file, "random_sample_rotation:%u\n", to_random_sample_rotation[ indices[11] ]);
		if(normal_types.size() > 0)
			fprintf(config_file, "normal_type:%u\n", normal_types[ indices[12] ]);

		file_ind++;
		fclose(config_file);

		if(!configs_remain(indices, index_sizes))
			break;

		update_indices(indices, index_sizes);
	}

	return 0;
}
