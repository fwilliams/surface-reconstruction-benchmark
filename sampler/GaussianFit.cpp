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

#include "GaussianFit.h"

GaussianFit::GaussianFit() : LevenbergMarquardt(50, 0.1)  {
	peak = 0;
	mean = 0;
	sigma = 0;

	has_been_fit = false;
}

GaussianFit::~GaussianFit()  {
	x_coords.clear();
	observations.clear();
}

void GaussianFit::clear_data()  {
	x_coords.clear();
	observations.clear();
}

void GaussianFit::add_point(double _x, double _o)  {
	if(has_been_fit)
		cout << "warning: already fit this guy!" << endl;
	x_coords.push_back(_x);
	observations.push_back(_o);
}

bool GaussianFit::do_fit(double _initPeak, double _initMean, double _initSigma)  {
	peak = _initPeak;
	mean = _initMean;
	sigma = _initSigma;

	if(x_coords.size() < 3)
		return false;

	BLinAlg::Vector* gaussian_solution = this->nonlinear_least_squares();
	double scaled_energy = (2.0*this->sqd_energy) / (double)(x_coords.size());
	if(this->squared_energy() > 0.1)  {
		delete gaussian_solution;
		return false;
	}

	peak = gaussian_solution->getEntry(0);
	mean = gaussian_solution->getEntry(1);
	sigma = gaussian_solution->getEntry(2);
	delete gaussian_solution;

	return true;
}

bool GaussianFit::do_fit()  {
	has_been_fit = true;
	if(x_coords.size() < 3)
		return false;

	double init_mean = 0, init_peak = 0, total_weight = 0;
	for(unsigned i = 0; i < x_coords.size(); i++)  {
		init_mean += observations[i]*x_coords[i];
		total_weight += observations[i];
		init_peak = observations[i] > init_peak ? observations[i] : init_peak;
	}
	init_mean /= total_weight;

	double init_std = 0;
	for(unsigned i = 0; i < x_coords.size(); i++)
		init_std += (x_coords[i]-init_mean)*(x_coords[i]-init_mean);
	init_std = sqrt(init_std / (double)x_coords.size());
	return this->do_fit(init_peak,init_mean,init_std);
}

BLinAlg::Vector* GaussianFit::initial_guess()  {
	BLinAlg::Vector* my_guess = new BLinAlg::Vector(3);
	my_guess->setEntry(0,peak);
	my_guess->setEntry(1,mean);
	my_guess->setEntry(2,sigma);
	return my_guess;
}

BLinAlg::Vector* GaussianFit::next_function(BLinAlg::Vector* _approximation)  {
	BLinAlg::Vector* my_function = new BLinAlg::Vector(x_coords.size());
	for(unsigned s = 0; s < x_coords.size(); s++)  {
		double gauss_func = this->local_evaluation(x_coords[s], _approximation);
		my_function->setEntry(s, (gauss_func-observations[s]));
	}
	return my_function;
}

DenseMatrix* GaussianFit::next_jacobian(BLinAlg::Vector* _approximation)  {
	DenseMatrix* my_jacobian = new DenseMatrix(x_coords.size(), 3);
	double local_peak = _approximation->getEntry(0), local_mean = _approximation->getEntry(1), local_sigma = _approximation->getEntry(2);

	for(unsigned s = 0; s < x_coords.size(); s++)  {
		double x_diff = x_coords[s]-local_mean;
		double exp_factor = exp(-(x_diff*x_diff)/(local_sigma*local_sigma));
		double peak_derivative = exp_factor;
		double mean_derivative = ((2.0*local_peak*x_diff)/(local_sigma*local_sigma)) * exp_factor;
		double sigma_derivative = ((2.0*local_peak*x_diff*x_diff)/(local_sigma*local_sigma*local_sigma)) * exp_factor;

		my_jacobian->setEntry(0,s, peak_derivative);
		my_jacobian->setEntry(1,s, mean_derivative);
		my_jacobian->setEntry(2,s, sigma_derivative);
	}

	return my_jacobian;
}

double GaussianFit::evaluation(double _x)  {
	double x_sqd = (_x-mean)*(_x-mean);
	return peak*exp(-x_sqd/(sigma*sigma));
}

double GaussianFit::local_evaluation(double _x, BLinAlg::Vector* _approximation)  {
	double local_peak = _approximation->getEntry(0), local_mean = _approximation->getEntry(1), local_sigma = _approximation->getEntry(2);
	double x_sqd = (_x-local_mean)*(_x-local_mean);
	return local_peak*exp(-x_sqd/(local_sigma*local_sigma));
}
