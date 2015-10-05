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

#ifndef RAYCASTER_H
#define RAYCASTER_H

#include "../modeling/Vector3.h"
#include "../modeling/RotationMatrix.h"

class RayCaster  {
	public:
		RayCaster()  {}
		RayCaster(Vector3 _position, Vector3 _view, Vector3 _up, int _resX, int _resY, double _near, double _fov);
		~RayCaster();

		void map_to_image(Vector3 _pos, double& x, double& y);
		void generatePerspectiveRay(double _x, double _y, Vector3& ray_pos, Vector3& ray_dir);
		void generateOrthographicRay(double _x, double _y, Vector3& ray_pos, Vector3& ray_dir);

		// ASSUMPTION: orthographic in x direction
		void generateOrthoPerspectiveRay(double _x, double _y, Vector3& ray_pos, Vector3& ray_dir);

		//void convert_local_system(Vector3 _pt, Vector3 _normal, RotationMatrix* _perturbation, Vector3& perturbedPt, Vector3& perturbedNormal);

		Vector3 getUBasis()  { return u_basis; }
		Vector3 getVBasis()  { return v_basis; }

		Vector3 camera_pos()  { return position; }
		int resx()  { return res_x; }
		int resy()  { return res_y; }

	private:
		Vector3 position, view, up, cross;
		int res_x, res_y;
		double near, fov;

		Vector3 u_param, v_param;
		Vector3 ll_start;

		Vector3 u_basis, v_basis, z_basis;

		void init_img_plane();
};

#endif
