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

#include "ImplicitNewton.h"

#include <iostream>
using namespace std;

ImplicitNewton::ImplicitNewton(ImplicitFunction* _implicit, double _newtonEps, double _hitEps)  {
	implicit_function = _implicit;
	newton_eps = _newtonEps;
	hit_eps = _hitEps;
	init(100);

	has_grid = false;
}

ImplicitNewton::ImplicitNewton(ImplicitFunction* _implicit, double _newtonEps, double _hitEps, double _lipschitz)  {
	implicit_function = _implicit;
	newton_eps = _newtonEps;
	hit_eps = _hitEps;
	lipschitz = _lipschitz;

	has_grid = false;
}

ImplicitNewton::ImplicitNewton(ImplicitFunction* _implicit, UniformOccupancyGrid* _grid, double _newtonEps, double _hitEps)  {
	implicit_function = _implicit;
	newton_eps = _newtonEps;
	hit_eps = _hitEps;
	init(100);

	binary_grid = _grid;
	has_grid = true;
}

ImplicitNewton::ImplicitNewton(ImplicitFunction* _implicit, UniformOccupancyGrid* _grid, double _newtonEps, double _hitEps, double _lipschitz)  {
	implicit_function = _implicit;
	newton_eps = _newtonEps;
	hit_eps = _hitEps;
	lipschitz = _lipschitz;

	binary_grid = _grid;
	has_grid = true;
}

ImplicitNewton::~ImplicitNewton()  {
	if(has_grid)
		delete binary_grid;
}

void ImplicitNewton::init(int _res)  {
	lipschitz = implicit_function->lipschitz_constant(_res);
	cout << "lipschitz: " << lipschitz << endl;
}

bool ImplicitNewton::intersection(Vector3 _rayPos, Vector3 _rayDir, Vector3& point)  {
	double t_hit = 0;
	if(!this->intersection(_rayPos, _rayDir, t_hit))
		return false;
	point = _rayPos + _rayDir*t_hit;
	return true;
}

bool ImplicitNewton::intersection(Vector3 _rayPos, Vector3 _rayDir, double& thit)  {
	bool debug = false;
	bool verbose_debug = false;
	bool tolerance_met = false;

	int max_iterations = 5000;
	int current_iteration = 0;

	double t_near, t_far;
	bool is_intersection = implicit_function->rayBoundIntersection(_rayPos, _rayDir, t_near, t_far);
	if(!is_intersection)  {
		if(debug) cout << "can't intersect!" << endl;
		return false;
	}

	if(debug) cout << "t near: " << t_near << " , t far: " << t_far << endl;

	double t_val = t_near;

	int num_grid_marching = 0;

	while(!tolerance_met)  {
		// first do check on if we should give up...
		if(current_iteration == max_iterations)  {
			if(debug) cout << "stopping at max iterations." << endl;
			break;
		}

		// if we escape the bounding box, then get out...
		if(t_val > t_far)  {
			if(debug) cout << "going outside of bound: " << current_iteration << endl;
			break;
		}

		Vector3 next_pt = _rayPos + _rayDir*t_val;

		// at this point, if we are using the uniform occupancy grid, then use it!
		if(has_grid && !binary_grid->isOccupied(next_pt))  {
			double t_delta = binary_grid->grid_march(next_pt, _rayDir)+binary_grid->getCellEps();
			t_val += t_delta;
			current_iteration++;
			num_grid_marching++;
			continue;
		}

		double next_func = implicit_function->function(next_pt);

		if(next_func > 0.0 && verbose_debug)
			cout << "we are INSIDE of the shape..."  << endl;
		if(MY_ABS(next_func) < hit_eps)  {
			if(debug) cout << "HIT! " << current_iteration << endl;
			tolerance_met = true;
			break;
		}	

		double abs_func = next_func < 0 ? -next_func : next_func;
		t_val += newton_eps*(abs_func / lipschitz);

		current_iteration++;
	}

	/*
	if(num_grid_marching > 0)
		cout << "num grid marching: " << num_grid_marching << endl;
		*/

	if(tolerance_met)  {
		if(debug) cout << "we got this shit!" << endl;
		thit = t_val;
		return true;
	}

	if(debug) cout << "failed..." << endl;
	return false;
}

void ImplicitNewton::intersection(Vector3 _rayPos, Vector3 _rayDir, double& thit, int& return_val)  {
	bool debug = false;
	bool verbose_debug = false;
	bool tolerance_met = false;

	int max_iterations = 5000;
	int current_iteration = 0;

	double t_near, t_far;
	bool is_intersection = implicit_function->rayBoundIntersection(_rayPos, _rayDir, t_near, t_far);
	if(!is_intersection)  {
		if(debug) cout << "can't intersect!" << endl;
		return_val = 0;
		return;
	}

	if(debug) cout << "t near: " << t_near << " , t far: " << t_far << endl;

	double t_val = t_near;

	int num_grid_marching = 0;

	while(!tolerance_met)  {
		// first do check on if we should give up...
		if(current_iteration == max_iterations)  {
			if(debug) cout << "stopping at max iterations." << endl;
			break;
		}

		// if we escape the bounding box, then get out...
		if(t_val > t_far)  {
			if(debug) cout << "going outside of bound: " << current_iteration << endl;
			break;
		}

		Vector3 next_pt = _rayPos + _rayDir*t_val;

		// at this point, if we are using the uniform occupancy grid, then use it!
		if(has_grid && !binary_grid->isOccupied(next_pt))  {
			double t_delta = binary_grid->grid_march(next_pt, _rayDir)+binary_grid->getCellEps();
			t_val += t_delta;
			current_iteration++;
			num_grid_marching++;
			continue;
		}

		double next_func = implicit_function->function(next_pt);

		if(next_func > 0.0 && verbose_debug)
			cout << "we are INSIDE of the shape..."  << endl;
		if(MY_ABS(next_func) < hit_eps)  {
			if(debug) cout << "HIT! " << current_iteration << endl;
			tolerance_met = true;
			break;
		}	

		double abs_func = next_func < 0 ? -next_func : next_func;
		t_val += newton_eps*(abs_func / lipschitz);

		current_iteration++;
	}

	/*
	if(num_grid_marching > 0)
		cout << "num grid marching: " << num_grid_marching << endl;
		*/

	if(tolerance_met)  {
		if(debug) cout << "we got this shit!" << endl;
		thit = t_val;
		return_val = 1;
		return;
	}

	if(debug) cout << "failed..." << endl;
	return_val = 0;
}
