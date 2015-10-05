/*
Copyright (c) 2010, Matt Berger and Bigyan Ankur Mukherjee
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

#ifndef UNIFORMSAMPLER_H
#define UNIFORMSAMPLER_H

#include "../evaluator/ComputationTimer.h"

#include "RayCaster.h"
#include "OrientedRangeImage.h"
#include "RangeImage.h"

#include "NPTSReader.h"
#include "PointCloud.h"

#include "RangeScanner.h"
#include "../modeling/ImplicitSphere.h"
#include "../modeling/RotationMatrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

#define NORMALS_ANALYTICAL   0
#define NORMALS_PCA_ORIENTED 1
#define NORMALS_PCA_MST      2

class UniformSampler  {
	public:
		UniformSampler(ImplicitFunction* _implicit, int _numScans, int _resX, int _resY);
		~UniformSampler();

		// all angle parameters in degrees

		void setMinRange(double _minRange)  { min_range = _minRange; }
		void setMaxRange(double _maxRange)  { max_range = _maxRange; }

		void set_num_stripes(int _numStripes)  { num_stripes = _numStripes; }
		void set_laser_fov(double _laserFOV)  { laser_fov = acos(-1.0)*(_laserFOV/180.0); }

		void set_noise(double _noise)  { additive_noise = _noise; }
		void set_laser_smoother(double _smoother)  { laser_smoother = _smoother; }
		void set_peak_threshold(double _threshold)  { peak_threshold = _threshold; }
		void set_std_threshold(double _threshold)  { std_threshold = _threshold; }

		void set_registration_error(double _registerScale)  {
			to_register = _registerScale > 0;
			registration_error = acos(-1.0)*(_registerScale/180.0);
		}

		void setRandomSampleRotation(bool _toRotate)  {
			to_random_sample_rotation = _toRotate;
		}

		void setNormalType(int _normalType)  {
			normal_type = _normalType;
		}

		void setPCAKNN(int _num)  {
			pca_knn = _num;
		}

		void set_stripe_dump(string _base)  {
			stripe_base = _base;
		}

		void sample();
		void dump_to_file(string _filename);

		bool getRigidTransformation(string _filename, DenseMatrix& rotation, Vector3& translation);

	private:
		ImplicitFunction* implicit_function;
		vector<OrientedRangeImage*> range_images;

		double my_pi;

		// required
		int num_scans;
		int res_x, res_y;

		// optional

		// range scanner parameters
		double min_range, max_range;
		int num_stripes;
		double laser_fov;
		double additive_noise;
		double laser_smoother;
		double peak_threshold, std_threshold;

		// registration parameters
		bool to_register;
		double registration_error;

		// normal parameters
		int normal_type;
		int pca_knn;

		// misc
		double baseline;
		bool to_random_sample_rotation;

		bool dump_range_images;
		string stripe_base;

		void dump_to_movie();
		void clear_range_images();
		void tokenize_line(vector<string>* tokens, const string& input, string sep);
};

#endif
