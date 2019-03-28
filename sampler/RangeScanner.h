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

#ifndef RANGESCANNER_H
#define RANGESCANNER_H

#include "RayCaster.h"
#include "OrientedRangeImage.h"
#include "ImplicitNewton.h"
#include "../modeling/Rigid.h"
#include "SparseRangeScan.h"
#include "PNGImage.h"

#include "../modeling/RotationMatrix.h"
#include "GaussianFit.h"

#define MY_FLOOR(x) ((int)x)
#define MY_CEIL(x) (x == ((int)x)) ? x : (int)(x+1)

struct RangePeak  {
	RangePeak()  {}
	RangePeak(double _x, double _t) : x(_x) , t(_t)  {}
	void set_peak(double _x, double _t)  {
		x = _x;
		t = _t;
	}
	double x, t;
};

class RangeSubset  {
	public:
		RangeSubset(int _resx, int _resy, double _gaussianSigma);
		~RangeSubset();
		void expand_range(int _x);
		void finalize_range();

		void contribute_point(int _x, int _y, double _intensity);
		void write_to_file(string _filename, double** _r, double** _g, double** _b);
		void smooth();

		int startX()  { return start_x; }
		int endX()  { return end_x; }
		int imageWidth() { return image_width; }
		bool containsNothing()  { return contains_nothing; }

		double get_intensity(int _x, int _y)  {
			int true_x = _x-start_x;
			return intensities[_y][true_x];
		}

	private:
		double** intensities;
		int start_x, end_x, image_width;
		int res_x, res_y;

		int scan_line;

		double gauss_sigma;
		double* gaussian_filter;
		int filter_size, filter_half;

		bool contains_nothing;
};

struct LaserStripe  {
	LaserStripe()  {}
	LaserStripe(Vector3 _nMinus, double _dMinus, Vector3 _n, double _d, Vector3 _nPlus, double _dPlus)  {
		n_minus = _nMinus;
		d_minus = _dMinus;
		n = _n;
		d = _d;
		n_plus = _nPlus;
		d_plus = _dPlus;
	}

	double minus_signed_distance(Vector3 _pt)  { return n_minus.dotProduct(_pt) + d_minus; }
	double plus_signed_distance(Vector3 _pt)  { return n_plus.dotProduct(_pt) + d_plus; }
	double pulse_abs_distance(Vector3 _pt)  { return fabs(n.dotProduct(_pt) + d); }

	Vector3 n_minus, n, n_plus;
	double d_minus, d, d_plus;
};

class RangeOctree  {
	public:
		RangeOctree(SparseRangeScan* _scan);
		RangeOctree(SparseRangeScan* _scan, int _depth);
		~RangeOctree();

		void finalize();
		void construct_tree();
		void add_point(int _p);

		void intersect_stripe(LaserStripe _stripe, vector<int>* overlapping_inds);
		bool stripe_aabb_intersection(LaserStripe _stripe);
		bool aabb_plane_intersection(Vector3 _n, double _d);

		bool is_leaf()  { return isLeaf; }
		bool is_root()  { return depth == 0; }

		int size()  { return point_inds.size(); }

	private:
		SparseRangeScan* scan;
		RangeOctree** children;

		Vector3 ll, ur;

		vector<int> point_inds;
		int depth;
		bool isLeaf;
};

class RangeScanner  {
	public:
		RangeScanner(ImplicitFunction* _shape, RayCaster _sensor, Vector3 _laserPos);
		RangeScanner(ImplicitFunction* _shape, RayCaster _sensor, Vector3 _laserPos, double _lipschitz);
		~RangeScanner();

		SparseRangeScan* optical_triangulation();

		// user-defined parameters
		void set_range(double _minRange, double _maxRange)  { min_range = _minRange; max_range = _maxRange; }
		void set_num_stripes(int _numStripes)  { num_stripes = _numStripes; }
		void set_laser_fov(double _laserFOV)  { laser_fov = _laserFOV; }

		void set_noise(double _noise)  { additive_noise = _noise; }
		void set_peak_threshold(double _threshold)  { peak_threshold = _threshold; }
		void set_std_threshold(double _threshold)  { std_threshold = _threshold; }
		void set_laser_smoother(double _smoother)  { laser_smoother = _smoother; }

		void set_base_stripe(string _baseStripe, int _stripeStart)  {
			base_stripe = _baseStripe;
			stripe_start = _stripeStart;
			dump_stripes = true;
		}

		int num_stripes_processed()  { return stripes_processed; }

		void ray_tracer(string _filename);
		void depth_image(SparseRangeScan* _scan, string _filename);

	private:
		ImplicitFunction* shape;
		double lipschitz;
		RayCaster sensor;
		Vector3 laser_pos;

		// user-defined parameters
		double min_range, max_range;
		int num_stripes;
		double laser_fov;

		double additive_noise;
		double peak_threshold, std_threshold;
		double laser_smoother;

		SparseRangeScan sparse_scan;
		SparseRangeScan image_scan;

		bool dump_stripes;
		string base_stripe;
		int stripes_processed, stripe_start;

		double sampleNormalDistribution(double _sigma, double _magnitude);
		double gaussian_beam(double _x, double _sigma);

		void interpolate_stripe(LaserStripe _stripe, Vector3 _pt, Vector3& n, double& d, double& w);
		void interpolate_stripe(Vector3 _startPt, double _shift, double _laserFOV, Vector3 _pt, Vector3& n, double& d, double& w);
		void construct_laser_stripe(Vector3 _startPt, double _shift, double _laserFOV, LaserStripe& stripe);
		void image_sensor();
};

#endif
