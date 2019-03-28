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

#include "RayCaster.h"

#include <iostream>
using namespace std;

RayCaster::RayCaster(Vector3 _position, Vector3 _view, Vector3 _up, int _resX, int _resY, double _near, double _fov) :
position(_position), view(_view), up(_up), res_x(_resX), res_y(_resY), near(_near), fov(_fov) {
	init_img_plane();
}

RayCaster::~RayCaster()  {
}

void RayCaster::map_to_image(Vector3 _pos, double& x, double& y)  {
	/*
	double s = 1.0;
	double fc_1 = 2813.36663/s, fc_2 = 2814.00245/s;
	double cc_x = 758.7654/s, cc_y = 1015.69320/s;

	DenseMatrix world_frame(3);
	world_frame.setEntry(0,0,cross.x); world_frame.setEntry(0,1,up.x); world_frame.setEntry(0,2,view.x);
	world_frame.setEntry(1,0,cross.y); world_frame.setEntry(1,1,up.y); world_frame.setEntry(1,2,view.y);
	world_frame.setEntry(2,0,cross.z); world_frame.setEntry(2,1,up.z); world_frame.setEntry(2,2,view.z);

	Vector3 camera_pos = world_frame.mult(_pos-position);
	double unit_x = camera_pos.x/camera_pos.z, unit_y = camera_pos.y/camera_pos.z;
	x = cc_x + fc_1*unit_x;
	y = cc_y + fc_2*unit_y;
	*/
}

void RayCaster::generatePerspectiveRay(double _x, double _y, Vector3& ray_pos, Vector3& ray_dir)  {
	double s = 1.0;

	double fc_1 = 2813.36269/s, fc_2 = 2813.99589/s;
	double cc_x = 758.75227/s, cc_y = 1015.65442/s;
	double z_c = 10.0;
	double x_c = (z_c/fc_1)*(_x-cc_x), y_c = (z_c/fc_2)*(_y-cc_y);

	Vector3 c_pos = Vector3(0,0,0);
	Vector3 c_ray_pos = Vector3(x_c,y_c,z_c);
	//cout << "<" << _x << ", " << _y << " > : " << c_ray_pos << endl;

	DenseMatrix camera_frame(3);
	camera_frame.setEntry(0,0,cross.x); camera_frame.setEntry(1,0,up.x); camera_frame.setEntry(2,0,view.x);
	camera_frame.setEntry(0,1,cross.y); camera_frame.setEntry(1,1,up.y); camera_frame.setEntry(2,1,view.y);
	camera_frame.setEntry(0,2,cross.z); camera_frame.setEntry(1,2,up.z); camera_frame.setEntry(2,2,view.z);

	Vector3 w_pos = camera_frame.mult(c_pos) + position;
	Vector3 w_ray_pos = camera_frame.mult(c_ray_pos) + position;
	Vector3 test_dir = w_ray_pos-w_pos;
	test_dir.normalize();

	ray_pos = position;
	Vector3 img_pos = ll_start + u_param*_x + v_param*_y;
	ray_dir = img_pos-ray_pos;
	ray_dir.normalize();

	//cout << "test dir: " << test_dir << " : ray dir: " << ray_dir << endl;
	//ray_dir = test_dir;
}

void RayCaster::generateOrthographicRay(double _x, double _y, Vector3& ray_pos, Vector3& ray_dir)  {
	Vector3 img_pos = ll_start + u_param*_x + v_param*_y;
	ray_pos = img_pos - view*near;
	ray_dir = img_pos-ray_pos;
	ray_dir.normalize();
}

void RayCaster::generateOrthoPerspectiveRay(double _x, double _y, Vector3& ray_pos, Vector3& ray_dir)  {
}

/*
void RayCaster::convert_local_system(Vector3 _pt, Vector3 _normal, RotationMatrix* _perturbation, Vector3& perturbedPt, Vector3& perturbedNormal)  {
	DenseMatrix local_basis(3);
	local_basis.setEntry(0,0, u_basis.x); local_basis.setEntry(1,0, u_basis.y); local_basis.setEntry(2,0, u_basis.z);
	local_basis.setEntry(0,1, v_basis.x); local_basis.setEntry(1,1, v_basis.y); local_basis.setEntry(2,1, v_basis.z);
	local_basis.setEntry(0,2, z_basis.x); local_basis.setEntry(1,2, z_basis.y); local_basis.setEntry(2,2, z_basis.z);

	Vector3 perturbed_u = _perturbation->rotate(u_basis);
	Vector3 perturbed_v = _perturbation->rotate(v_basis);
	Vector3 perturbed_z = _perturbation->rotate(z_basis);
	DenseMatrix inv_perturbed_basis(3);
	inv_perturbed_basis.setEntry(0,0, perturbed_u.x); inv_perturbed_basis.setEntry(1,0, perturbed_v.x); inv_perturbed_basis.setEntry(2,0, perturbed_z.x);
	inv_perturbed_basis.setEntry(0,1, perturbed_u.y); inv_perturbed_basis.setEntry(1,1, perturbed_v.y); inv_perturbed_basis.setEntry(2,1, perturbed_z.y);
	inv_perturbed_basis.setEntry(0,2, perturbed_u.z); inv_perturbed_basis.setEntry(1,2, perturbed_v.z); inv_perturbed_basis.setEntry(2,2, perturbed_z.z);

	// bring point,normal into local coordinate system
	Vector3 local_pt = local_basis.mult(_pt-position);
	Vector3 local_normal = local_basis.mult(_normal);

	// and project back, perturbed
	perturbedPt = inv_perturbed_basis.mult(local_pt) + position;
	perturbedNormal = inv_perturbed_basis.mult(local_normal);
}
*/

void RayCaster::init_img_plane()  {
	// normalize view and up, just in case they weren't already
	view.normalize();
	up.normalize();

	// get other crossed vector on image plane
	cross = view.cross(up);

	// image plane stuff
	double aspect_ratio = (double)res_x / (double)res_y;
	double plane_width = 2.0*near*tan(0.5*fov);
	double plane_height = plane_width / aspect_ratio;

	Vector3 img_center = position + view*near;
	Vector3 ll_corner = img_center + cross*(-0.5*plane_width) + up*(-0.5*plane_height);

	u_param = cross * plane_width * (1.0 / (double)res_x);
	v_param = up * plane_height * (1.0 / (double)res_y);
	ll_start = ll_corner + u_param*0.5 + v_param*0.5;

	u_basis = u_param;
	u_basis.normalize();
	v_basis = v_param;
	v_basis.normalize();
	z_basis = v_basis.cross(u_basis);
}
