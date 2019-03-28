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

#include "SphericalImager.h"

#include <iostream>
using namespace std;

SphericalImager::SphericalImager(ImplicitFunction* _implicit, int _numScans, int _resX, int _resY)  {
	implicit_function = _implicit;

	num_scans = _numScans;
	res_x = _resX;
	res_y = _resY;

	my_pi = acos(-1.0);

	min_range = 60;

	to_random_sample_rotation = false;
}

SphericalImager::~SphericalImager()  {
}

void SphericalImager::sample(string _basefile)  {
	cout << "------ scanner parameters ------" << endl;
	cout << "   num scans : " << num_scans << endl;
	cout << "   image resolution : " << res_x << ", " << res_y << endl;
	cout << "   sensor range : " << min_range << endl;
	cout << "   do random sample rotation ? " << to_random_sample_rotation << endl;

	Vector3 ll_bound, ur_bound;
	implicit_function->bounds(ll_bound, ur_bound);

	Vector3 diag_bound = ur_bound-ll_bound;
	Vector3 center = (ur_bound+ll_bound)*0.5;
	double sampling_radius = min_range + 0.5*diag_bound.length();

	time_t cur_time;
	time(&cur_time);
	srand(cur_time);
	for(int i = 0; i < 20; i++)
		rand();

	vector<Vector3> uniform_sampling;
	ostringstream num_samples;
	num_samples << num_scans;
	string particle_file = "./particle_sampler/particles"+num_samples.str()+".npts";
	OrientedPointCloud sphere_pc;
	NPTSReader npts_reader(particle_file);

	if(!npts_reader.read_pc(&sphere_pc))  {
		cout << "couldn't open file " << particle_file << " ; instead, going to sample ..." << endl;
		ImplicitSphere* sampling_sphere = new ImplicitSphere(1);
		uniform_sampling = sampling_sphere->sample_surface(num_scans, 0, false);
	}
	else  {
		double rx = (double)rand() / (double)RAND_MAX;
		double ry = (double)rand() / (double)RAND_MAX;
		double rz = (double)rand() / (double)RAND_MAX;
		double random_angle = 2.0*acos(-1)*((double)rand() / (double)RAND_MAX);
		Vector3 random_axis(rx, ry, rz);
		random_axis.normalize();
		RotationMatrix sphere_rotation(random_axis, random_angle);

		for(int i = 0; i < sphere_pc.size(); i++)  {
			Vector3 sphere_pt = sphere_pc.getPoint(i);
			if(to_random_sample_rotation)  {
				Vector3 rotated_pt = sphere_rotation.rotate(sphere_pt);
				uniform_sampling.push_back(rotated_pt);
			}
			else
				uniform_sampling.push_back(sphere_pt);
		}
	}
	if(uniform_sampling.size() != num_scans)
		cout << "ERROR: uniform sampling not at prescribed size!" << endl;

	double lipschitz = implicit_function->lipschitz_constant(100);

	for(int s = 0; s < num_scans; s++)  {
		ComputationTimer timer("spherical imager");
		timer.start();

		Vector3 sphere_sample = uniform_sampling[s];
		Vector3 look_from = center + sphere_sample*sampling_radius;
		Vector3 look_at = center;
		Vector3 view_dir = look_at-look_from;
		view_dir.normalize();

		Vector3 up;
		if(view_dir.x == 0 && view_dir.y == 1 && view_dir.z == 0)
			up = Vector3(0,0,1);
		else  {
			Vector3 fake_up(0,1,0);
			Vector3 crossed = view_dir.cross(fake_up);
			crossed.normalize();
			up = crossed.cross(view_dir);
			up.normalize();
		}

		double negate_up = (double)rand() / (double)RAND_MAX;
		up = negate_up < 0.5 ? up*-1.0 : up;

		cout << "scan " << s << " / " << num_scans << endl;
		RayCaster ray_caster(look_from, view_dir, up, res_x, res_y, 50, my_pi/4.0);
		Vector3 laser_pos = look_from + ray_caster.getUBasis()*70.0;

		char* next_file = new char[300];
		sprintf(next_file, "%s-%u.png", _basefile.c_str(), s);

		RangeScanner range_scanner(implicit_function, ray_caster, laser_pos, lipschitz);
		range_scanner.set_range(min_range, 1e10);
		range_scanner.ray_tracer(next_file);
		delete [] next_file;

		timer.end();
		cout << timer.getComputation() << " : " << timer.getElapsedTime() << "s" << endl;
	}

}
