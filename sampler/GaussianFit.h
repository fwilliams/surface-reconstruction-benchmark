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

#ifndef GAUSSIANFIT_H
#define GAUSSIANFIT_H

#include "../modeling/lm_method.h"
#include <vector>
using namespace std;

class GaussianFit : public LevenbergMarquardt  {
	public:
		GaussianFit();
		~GaussianFit();

		void clear_data();
		void add_point(double _x, double _o);
		bool do_fit(double _initPeak, double _initMean, double _initSigma);
		bool do_fit();

		double get_peak()  { return peak; }
		double get_mean()  { return mean; }
		double get_std()  { return fabs(sigma); }
		int num_observations()  { return observations.size(); }
		bool hasBeenFit()  { return has_been_fit; }

		virtual BLinAlg::Vector* initial_guess();
		virtual BLinAlg::Vector* next_function(BLinAlg::Vector* _approximation);
		virtual DenseMatrix* next_jacobian(BLinAlg::Vector* _approximation);

		double evaluation(double _x);
		double local_evaluation(double _x, BLinAlg::Vector* _approximation);

	private:
		vector<double> x_coords;
		vector<double> observations;
		double peak, mean, sigma;

		bool has_been_fit;
};

#endif
