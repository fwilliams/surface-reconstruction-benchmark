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

#include "lm_method.h"

LevenbergMarquardt::LevenbergMarquardt(int _solverIterations, float _energyEps)  {
	solver_iterations = _solverIterations;
	energy_eps = _energyEps;
}

LevenbergMarquardt::LevenbergMarquardt()  {
	solver_iterations = 25;
	energy_eps = 0.1;
}

LevenbergMarquardt::~LevenbergMarquardt()  {
}

BLinAlg::Vector* LevenbergMarquardt::nonlinear_least_squares()  {
	// get initial guess
	BLinAlg::Vector* solution = this->initial_guess();
	double damping = 0.0;

	double velta = 2.0, beta = 2.0, gamma = 3.0;

	// for each iteration ...
	for(int i = 0; i < solver_iterations; i++)  {
		// construct function & jacobian at this next guess
		BLinAlg::Vector* next_func = this->next_function(solution);
		sqd_energy = 0.5*next_func->dot(next_func);

		// check to see if we are satisfied
		double scaled_energy = (2.0*sqd_energy) / (double)(next_func->dim());
		if( scaled_energy < (energy_eps*energy_eps) )  {
			delete next_func;
			break;
		}

		DenseMatrix* next_jacob = this->next_jacobian(solution);
		DenseMatrix* jacob_transpose = next_jacob->transpose();

		// setup system
		BLinAlg::Vector* rhs = jacob_transpose->mult(next_func);
		rhs->scale(-1.0);

		DenseMatrix* normal_matrix = next_jacob->normalMatrix();

		// for first iteration, follow [Eigensatz et al. 08]: largest diagonal of normal matrix * 10^-6
		if(i == 0)  {
			double largest_diagonal = -1e10;
			for(int r = 0; r < normal_matrix->rows(); r++)  {
				double diagonal = normal_matrix->getEntry(r,r);
				largest_diagonal = diagonal > largest_diagonal ? diagonal : largest_diagonal;
			}
			damping = largest_diagonal*10e-6;
		}

		DenseMatrix* damping_matrix = new DenseMatrix(normal_matrix->rows());
		damping_matrix->scale(damping);
		//damping_matrix->scale(0.0);
		normal_matrix->add(damping_matrix);

		// solve
		BLinAlg::Vector* delta = normal_matrix->linearSystem(rhs);

		// now follow schedule of [Nielsen et al. 99] for updates
		{
			BLinAlg::Vector* energy_gradient = jacob_transpose->mult(next_func);
			BLinAlg::Vector* dampened_delta = new BLinAlg::Vector(delta);
			dampened_delta->scale(damping);
			dampened_delta->sub(energy_gradient);

			BLinAlg::Vector* new_solution = new BLinAlg::Vector(solution);
			new_solution->add(delta);
			BLinAlg::Vector* new_next_func = this->next_function(new_solution);
			double new_sqd_energy = 0.5*new_next_func->dot(new_next_func);

			double gain = (sqd_energy - new_sqd_energy) / (0.5*delta->dot(dampened_delta));
			delete energy_gradient;
			delete dampened_delta;
			delete new_next_func;
			delete new_solution;

			// if gain > 0, then we have decreased the energy: proceed ...
			if(gain > 0.0)  {
				solution->add(delta);
				double inv_gamma = 1.0 / gamma;
				double mess = 1.0 - (beta-1.0)*pow((2.0*gain-1.0), 3);
				damping = inv_gamma > mess ? inv_gamma : mess;
				velta = beta;
			}

			// otherwise, up the gradient descent portion for stability ...
			else  {
				damping *= velta;
				velta *= 2.0;
			}
		}

		delete next_func;
		delete next_jacob;
		delete jacob_transpose;
		delete rhs;
		delete normal_matrix;
		delete damping_matrix;
		delete delta;
	}

	// return us the solution!
	return solution;
}
