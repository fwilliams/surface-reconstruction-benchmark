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

#include "ImplicitFunction.h"
#include <cstdlib>
#include <cstdio>

#include "../evaluator/ComputationTimer.h"

using namespace std;

void ImplicitFunction::generateTangentPlane(Vector3 _pt, Vector3& vec1, Vector3& vec2)  {
	// generate a random point on the sphere centered at _pt
	double rand_num1 = (double)rand() / (double)RAND_MAX;
	double rand_num2 = (double)rand() / (double)RAND_MAX;
	double z = (1.0-rand_num1)*-1.0 + (rand_num1)*1.0;
	double t = (rand_num2)*2.0*acos(-1.0);
	double x = sqrt(1.0-z*z)*cos(t);
	double y = sqrt(1.0-z*z)*sin(t);
	Vector3 rand_pt = Vector3(x,y,z) + _pt;

	// project point onto tangent plane
	Vector3 pt_normal = this->normal(_pt);
	Vector3 projected_pt = rand_pt + pt_normal*(pt_normal.dotProduct(_pt-rand_pt));

	// and form orthogonal basis of tangent plane
	vec1 = projected_pt-_pt;
	vec1.normalize();
	vec2 = vec1.cross(pt_normal);
}

DenseMatrix ImplicitFunction::weingarten_map(Vector3 _pt)  {
	// generate tangent plane
	Vector3 tangent_vec1, tangent_vec2;
	this->generateTangentPlane(_pt, tangent_vec1, tangent_vec2);

	// compute gradient & hessian at this point
	Vector3 pt_gradient = this->gradient(_pt);
	DenseMatrix pt_hessian = this->hessian(_pt);

	// and now construct each entry of weingarten matrix
	double inv_gradient = 1.0 / pt_gradient.length();
	double w_00 = inv_gradient*tangent_vec1.dotProduct( pt_hessian.mult(tangent_vec1) );
	double w_01 = inv_gradient*tangent_vec1.dotProduct( pt_hessian.mult(tangent_vec2) );
	double w_10 = w_01;
	double w_11 = inv_gradient*tangent_vec2.dotProduct( pt_hessian.mult(tangent_vec2) );

	DenseMatrix weingarten_matrix(2);
	weingarten_matrix.setEntry(0,0, w_00);
	weingarten_matrix.setEntry(0,1, w_01);
	weingarten_matrix.setEntry(1,0, w_10);
	weingarten_matrix.setEntry(1,1, w_11);

	return weingarten_matrix;
}

ImplicitCurvatures ImplicitFunction::compute_curvatures(Vector3 _pt)  {
	// generate tangent plane
	Vector3 tangent_vec1, tangent_vec2;
	this->generateTangentPlane(_pt, tangent_vec1, tangent_vec2);

	// compute gradient & hessian at this point
	Vector3 pt_gradient = this->gradient(_pt);
	DenseMatrix pt_hessian = this->hessian(_pt);

	// and now construct each entry of weingarten matrix
	double inv_gradient = 1.0 / pt_gradient.length();
	double w_00 = inv_gradient*tangent_vec1.dotProduct( pt_hessian.mult(tangent_vec1) );
	double w_01 = inv_gradient*tangent_vec1.dotProduct( pt_hessian.mult(tangent_vec2) );
	double w_10 = w_01;
	double w_11 = inv_gradient*tangent_vec2.dotProduct( pt_hessian.mult(tangent_vec2) );

	DenseMatrix weingarten_matrix(2);
	weingarten_matrix.setEntry(0,0, w_00);
	weingarten_matrix.setEntry(0,1, w_01);
	weingarten_matrix.setEntry(1,0, w_10);
	weingarten_matrix.setEntry(1,1, w_11);
	DenseEigenSystem eigen_system(&weingarten_matrix, 0, 1);
	eigen_system.solve();

	// TODO: get proper principal directions
	return ImplicitCurvatures(eigen_system.getEigenvalue(0), eigen_system.getEigenvalue(1));
}

Vector3 ImplicitFunction::projection(Vector3 _pt, double _eps, int _maxIterations, double _lambda)  {
	Vector3 ll, ur;
	this->bounds(ll,ur);
	double eps_bound = (ll-ur).length()*_eps;
	Vector3 pt = _pt;

	for(int i = 0; i < _maxIterations; i++)  {
		double pt_function = this->function(pt);
		Vector3 pt_gradient = this->gradient(pt);

		double approx_distance = fabs( (pt_function / pt_gradient.length()) );
		if(approx_distance < eps_bound)
			break;

		Vector3 pt_normal = pt_gradient;
		pt_normal.normalize();
		pt_normal.scale(-(pt_function*_lambda)/pt_gradient.length());
		pt = pt + pt_normal;
	}

	return pt;
}

float ImplicitFunction::spatial_derivative(float*** _volume, Vector3 _voxel, int _dim, float _step, int _res)  {
	int x = _voxel.x, y = _voxel.y, z = _voxel.z;

	if(_dim == 0)  {
		if(x == 0)
			return (_volume[z][y][1]-_volume[z][y][0])/_step;
		else if(x == _res-1)
			return (_volume[z][y][_res-1]-_volume[z][y][_res-2])/_step;
		else
			return (_volume[z][y][x+1]-_volume[z][y][x-1])/(2.0*_step);
	}
	else if(_dim == 1)  {
		if(y == 0)
			return (_volume[z][1][x]-_volume[z][0][x])/_step;
		else if(y == _res-1)
			return (_volume[z][_res-1][x]-_volume[z][_res-2][x])/_step;
		else
			return (_volume[z][y+1][x]-_volume[z][y-1][x])/(2.0*_step);
	}
	else  {
		if(z == 0)
			return (_volume[1][y][x]-_volume[0][y][x])/_step;
		else if(z == _res-1)
			return (_volume[_res-1][y][x]-_volume[_res-2][y][x])/_step;
		else
			return (_volume[z+1][y][x]-_volume[z-1][y][x])/(2.0*_step);
	}
}

void ImplicitFunction::write_to_volume(string _fileout, int _maxRes, double _epsBound)  {
	this->write_to_volume(_fileout, 64, _maxRes, _epsBound);
}

double ImplicitFunction::lipschitz_constant(int _res)  {
	Vector3 minBound, maxBound;
	this->bounds(minBound, maxBound);

	int num_dots = 10;
	int thresh_dot = _res / num_dots;
	double lipschitz = -1e10;

	cout << "computing lipschitz constant " << flush;
	for(int z = 0; z < _res; z++)  {
		if((z % thresh_dot) == 0)
			cout << "." << flush;
		double z_interp = (double)z/(double)(_res-1);
		double new_z = minBound.z*(1.0-z_interp) + maxBound.z*z_interp;

		for(int y = 0; y < _res; y++)  {
			double y_interp = (double)y/(double)(_res-1);
			double new_y = minBound.y*(1.0-y_interp) + maxBound.y*y_interp;

			for(int x = 0; x < _res; x++)  {
				double x_interp = (double)x/(double)(_res-1);
				double new_x = minBound.x*(1.0-x_interp) + maxBound.x*x_interp;
				Vector3 new_pt(new_x,new_y,new_z);
				Vector3 grad_value = this->gradient(new_pt);
				double grad_mag = grad_value.length();
				lipschitz = grad_mag > lipschitz ? grad_mag : lipschitz;
			}
		}
	}
	cout << " : " << lipschitz << endl;

	return lipschitz;
}

void ImplicitFunction::write_to_volume(string _fileout, int _minRes, int _maxRes, double _epsBound)  {
	Vector3 min, max;
	this->bounds(min, max);

	printf("shape center: %.10f %.10f %.10f\n", min.x, min.y, min.z);

	// expand bb by the epsilon
	Vector3 epsilon = (max-min)*(0.5*_epsBound);
	min = min - epsilon;
	max = max + epsilon;

	// ranges of bounds
	Vector3 bounds = max - min;

	// from the bounds, compute resolutions
	int res[3];
	if(bounds.x > bounds.y && bounds.x > bounds.z)  {
		res[0] = _maxRes;
		res[1] = _maxRes*(bounds.y / bounds.x);
		res[2] = _maxRes*(bounds.z / bounds.x);
	}
	else if(bounds.y > bounds.z)  {
		res[1] = _maxRes;
		res[0] = _maxRes*(bounds.x / bounds.y);
		res[2] = _maxRes*(bounds.z / bounds.y);
	}
	else  {
		res[2] = _maxRes;
		res[0] = _maxRes*(bounds.x / bounds.z);
		res[1] = _maxRes*(bounds.y / bounds.z);
	}

	for(int i = 0; i < 3; i++)  {
		res[i] = res[i] < _minRes ? _minRes : res[i];
	}

	Vector3 inv_res(1.0/res[0], 1.0/res[1], 1.0/res[2]);

	cout << "resolution: <" << res[0] << ", " << res[1] << ", " << res[2] << ">" << endl;

	Vector3 voxel_res(bounds);
	voxel_res.x *= inv_res.x;
	voxel_res.y *= inv_res.y;
	voxel_res.z *= inv_res.z;
	printf("voxel resolution: %.10f %.10f %.10f\n", voxel_res.x, voxel_res.y, voxel_res.z);

	// compute interval length and sub-center start
	Vector3 intervals = Vector3(bounds.x*inv_res.x, bounds.y*inv_res.y, bounds.z*inv_res.z);
	Vector3 subcenter = min + intervals*0.5;

	char buffer[256];
   sprintf(buffer, "%s-%dx%dx%d.raw", _fileout.c_str(), res[0], res[1], res[2]);
   FILE* raw_file = fopen(buffer, "wb");
	for(int z = 0; z < res[2]; z++)  {
		ComputationTimer timer("volume");
		timer.start();
		double z_value = subcenter.z + (z*intervals.z);
      cout << "Progress: " << z << "/" << res[2] << "\r";
      cout.flush();
		for(int y = 0; y < res[1]; y++)  {
			double y_value = subcenter.y + (y*intervals.y);
			for(int x = 0; x < res[0]; x++)  {
				double x_value = subcenter.x + (x*intervals.x);
				double func_value = this->function(Vector3(x_value, y_value, z_value));
				//cout << Vector3(x_value, y_value, z_value) << " : " << func_value << endl;
				 fwrite(&func_value, sizeof(double), 1, raw_file);
			}
		}	
		timer.end();
		cout << timer.getComputation() << " : " << timer.getElapsedTime() << "s" << endl;
	}
   cout << endl;

  

   fclose(raw_file);

   cout << "Dumping finished to file: " << buffer << endl;
}

void ImplicitFunction::write_kitten_kaboodle(string _base, int _maxRes, double _epsBound)  {
	Vector3 min, max;
	this->bounds(min, max);

	// expand bb by the epsilon
	Vector3 epsilon = (max-min)*(0.5*_epsBound);
	min = min - epsilon;
	max = max + epsilon;

	// ranges of bounds
	Vector3 bounds = max - min;

	// from the bounds, compute resolutions
	int res[3];
	if(bounds.x > bounds.y && bounds.x > bounds.z)  {
		res[0] = _maxRes;
		res[1] = _maxRes*(bounds.y / bounds.x);
		res[2] = _maxRes*(bounds.z / bounds.x);
	}
	else if(bounds.y > bounds.z)  {
		res[1] = _maxRes;
		res[0] = _maxRes*(bounds.x / bounds.y);
		res[2] = _maxRes*(bounds.z / bounds.y);
	}
	else  {
		res[2] = _maxRes;
		res[0] = _maxRes*(bounds.x / bounds.z);
		res[1] = _maxRes*(bounds.y / bounds.z);
	}

	Vector3 inv_res(1.0/res[0], 1.0/res[1], 1.0/res[2]);

	Vector3 voxel_res(bounds);
	voxel_res.x *= inv_res.x;
	voxel_res.y *= inv_res.y;
	voxel_res.z *= inv_res.z;

	// compute interval length and sub-center start
	Vector3 intervals = Vector3(bounds.x*inv_res.x, bounds.y*inv_res.y, bounds.z*inv_res.z);
	Vector3 subcenter = min + intervals*0.5;

	// first dump out meta data
	string meta_vol_name = _base + ".meta";
	FILE* meta_vol_file = fopen(meta_vol_name.c_str(), "wb");
	cout << "writing out meta data to " << meta_vol_name << " ..." << endl;
	float min_x = min.x, min_y = min.y, min_z = min.z;
	float max_x = max.x, max_y = max.y, max_z = max.z;
	int res_x = res[0], res_y = res[1], res_z = res[2];
	fwrite(&min_x, sizeof(float), 1, meta_vol_file);
	fwrite(&min_y, sizeof(float), 1, meta_vol_file);
	fwrite(&min_z, sizeof(float), 1, meta_vol_file);
	fwrite(&max_x, sizeof(float), 1, meta_vol_file);
	fwrite(&max_y, sizeof(float), 1, meta_vol_file);
	fwrite(&max_z, sizeof(float), 1, meta_vol_file);
	fwrite(&res_x, sizeof(int), 1, meta_vol_file);
	fwrite(&res_y, sizeof(int), 1, meta_vol_file);
	fwrite(&res_z, sizeof(int), 1, meta_vol_file);
	fclose(meta_vol_file);

	// setup function & derivative files
	string func_vol_name = _base + ".vol";
	string dx_vol_name = _base + "_dx.vol";
	string dy_vol_name = _base + "_dy.vol";
	string dz_vol_name = _base + "_dz.vol";

	FILE* func_vol_file = fopen(func_vol_name.c_str(), "wb");
	FILE* dx_vol_file = fopen(dx_vol_name.c_str(), "wb");
	FILE* dy_vol_file = fopen(dy_vol_name.c_str(), "wb");
	FILE* dz_vol_file = fopen(dz_vol_name.c_str(), "wb");

	// write out the function & gradient
	float*** grad_x_volume = new float**[ res[2] ];
	float*** grad_y_volume = new float**[ res[2] ];
	float*** grad_z_volume = new float**[ res[2] ];
	cout << "creating gradient volumes ..." << endl;
	for(int z = 0; z < res[2]; z++)  {
		grad_x_volume[z] = new float*[ res[1] ];
		grad_y_volume[z] = new float*[ res[1] ];
		grad_z_volume[z] = new float*[ res[1] ];
		for(int y = 0; y < res[1]; y++)  {
			grad_x_volume[z][y] = new float[ res[0] ];
			grad_y_volume[z][y] = new float[ res[0] ];
			grad_z_volume[z][y] = new float[ res[0] ];
		}
	}
	cout << "... done creating gradient volumes ..." << endl;

	for(int z = 0; z < res[2]; z++)  {
		double z_value = subcenter.z + (z*intervals.z);
      cout << "Progress: " << z << "/" << res[2] << "\r";
      cout.flush();
		for(int y = 0; y < res[1]; y++)  {
			double y_value = subcenter.y + (y*intervals.y);
			for(int x = 0; x < res[0]; x++)  {
				double x_value = subcenter.x + (x*intervals.x);
				Vector3 pt(x_value, y_value, z_value);
				float func_value = this->function(pt);
				Vector3 grad_value = this->gradient(pt);
				float grad_x = grad_value.x, grad_y = grad_value.y, grad_z = grad_value.z;
				grad_x_volume[z][y][x] = grad_x;
				grad_y_volume[z][y][x] = grad_y;
				grad_z_volume[z][y][x] = grad_z;

				fwrite(&func_value, sizeof(float), 1, func_vol_file);
				fwrite(&grad_x, sizeof(float), 1, dx_vol_file);
				fwrite(&grad_y, sizeof(float), 1, dy_vol_file);
				fwrite(&grad_z, sizeof(float), 1, dz_vol_file);
			}
		}	
	}
   cout << endl;

	// close up function & derivative
	fclose(func_vol_file);
	fclose(dx_vol_file);
	fclose(dy_vol_file);
	fclose(dz_vol_file);

	// now take care of second-order derivatives
	string dxx_vol_name = _base + "_dxx.vol", dxy_vol_name = _base + "_dxy.vol", dxz_vol_name = _base + "_dxz.vol";
	string dyy_vol_name = _base + "_dyy.vol", dyz_vol_name = _base + "_dyz.vol", dzz_vol_name = _base + "_dzz.vol";
	FILE* dxx_vol_file = fopen(dxx_vol_name.c_str(), "wb");
	FILE* dxy_vol_file = fopen(dxy_vol_name.c_str(), "wb");
	FILE* dxz_vol_file = fopen(dxz_vol_name.c_str(), "wb");
	FILE* dyy_vol_file = fopen(dyy_vol_name.c_str(), "wb");
	FILE* dyz_vol_file = fopen(dyz_vol_name.c_str(), "wb");
	FILE* dzz_vol_file = fopen(dzz_vol_name.c_str(), "wb");

	for(int z = 0; z < res[2]; z++)  {
		double z_value = subcenter.z + (z*intervals.z);
      cout << "Progress: " << z << "/" << res[2] << "\r";
      cout.flush();
		for(int y = 0; y < res[1]; y++)  {
			double y_value = subcenter.y + (y*intervals.y);
			for(int x = 0; x < res[0]; x++)  {
				double x_value = subcenter.x + (x*intervals.x);
				Vector3 voxel(x,y,z);

				float dxx = this->spatial_derivative(grad_x_volume, voxel, 0, intervals.x, res[0]);
				float dxy = this->spatial_derivative(grad_x_volume, voxel, 1, intervals.y, res[1]);
				float dxz = this->spatial_derivative(grad_x_volume, voxel, 2, intervals.z, res[2]);
				float dyy = this->spatial_derivative(grad_y_volume, voxel, 1, intervals.y, res[1]);
				float dyz = this->spatial_derivative(grad_y_volume, voxel, 2, intervals.z, res[2]);
				float dzz = this->spatial_derivative(grad_z_volume, voxel, 2, intervals.z, res[2]);

				fwrite(&dxx, sizeof(float), 1, dxx_vol_file);
				fwrite(&dxy, sizeof(float), 1, dxy_vol_file);
				fwrite(&dxz, sizeof(float), 1, dxz_vol_file);
				fwrite(&dyy, sizeof(float), 1, dyy_vol_file);
				fwrite(&dyz, sizeof(float), 1, dyz_vol_file);
				fwrite(&dzz, sizeof(float), 1, dzz_vol_file);
			}
		}	
	}
   cout << endl;

	// close up second-order derivatives
	fclose(dxx_vol_file);
	fclose(dxy_vol_file);
	fclose(dxz_vol_file);
	fclose(dyy_vol_file);
	fclose(dyz_vol_file);
	fclose(dzz_vol_file);
}

BasicTriMesh* ImplicitFunction::isosurface(int _minRes, int _maxRes, double _epsBound)  {
	Vector3 min, max;
	this->bounds(min, max);

	// expand bb by the epsilon
	Vector3 epsilon = (max-min)*(0.5*_epsBound);
	min = min - epsilon;
	max = max + epsilon;

	// ranges of bounds
	Vector3 bounds = max - min;

	// from the bounds, compute resolutions
	int res[3];
	if(bounds.x > bounds.y && bounds.x > bounds.z)  {
		res[0] = _maxRes;
		res[1] = _maxRes*(bounds.y / bounds.x);
		res[2] = _maxRes*(bounds.z / bounds.x);
	}
	else if(bounds.y > bounds.z)  {
		res[1] = _maxRes;
		res[0] = _maxRes*(bounds.x / bounds.y);
		res[2] = _maxRes*(bounds.z / bounds.y);
	}
	else  {
		res[2] = _maxRes;
		res[0] = _maxRes*(bounds.x / bounds.z);
		res[1] = _maxRes*(bounds.y / bounds.z);
	}

	Vector3 inv_res(1.0/res[0], 1.0/res[1], 1.0/res[2]);

	Vector3 voxel_res(bounds);
	voxel_res.x *= inv_res.x;
	voxel_res.y *= inv_res.y;
	voxel_res.z *= inv_res.z;

	// compute interval length and sub-center start
	Vector3 intervals = Vector3(bounds.x*inv_res.x, bounds.y*inv_res.y, bounds.z*inv_res.z);
	Vector3 subcenter = min + intervals*0.5;

	float*** scalar_field = new float**[ res[2] ];
	for(int z = 0; z < res[2]; z++)  {
		scalar_field[z] = new float*[ res[1] ];
		for(int y = 0; y < res[1]; y++)
			scalar_field[z][y] = new float[ res[0] ];
	}

	int num_dots = 10;
	int thresh_dot = res[2] / num_dots;
	for(int z = 0; z < res[2]; z++)  {
		if((z % thresh_dot) == 0)
			cout << "." << flush;
		double z_value = subcenter.z + (z*intervals.z);
		for(int y = 0; y < res[1]; y++)  {
			double y_value = subcenter.y + (y*intervals.y);
			for(int x = 0; x < res[0]; x++)  {
				double x_value = subcenter.x + (x*intervals.x);
				Vector3 pt(x_value, y_value, z_value);
				scalar_field[z][y][x] = this->function(pt);
			}
		}	
	}
	cout << endl;

	std::map<mc_edge, int> vertex_hash;
	BasicTriMesh* mc_mesh = new BasicTriMesh();

	for(int z = 0; z < res[2]-1; z++)  {
		double z_value = subcenter.z + (z*intervals.z);
		double next_z_value = subcenter.z + ((z+1)*intervals.z);
		for(int y = 0; y < res[1]-1; y++)  {
			double y_value = subcenter.y + (y*intervals.y);
			double next_y_value = subcenter.y + ((y+1)*intervals.y);
			for(int x = 0; x < res[0]-1; x++)  {
				double x_value = subcenter.x + (x*intervals.x);
				double next_x_value = subcenter.x + ((x+1)*intervals.x);

				double f_0 = scalar_field[z][y][x]      , f_1 = scalar_field[z][y][x+1];
				double f_2 = scalar_field[z+1][y][x+1]  , f_3 = scalar_field[z+1][y][x];
				double f_4 = scalar_field[z][y+1][x]    , f_5 = scalar_field[z][y+1][x+1];
				double f_6 = scalar_field[z+1][y+1][x+1], f_7 = scalar_field[z+1][y+1][x];
				int v_0 = x+y*res[0]+z*(res[0]*res[1]), v_1 = (x+1)+y*res[0]+z*(res[0]*res[1]);
				int v_2 = (x+1)+y*res[0]+(z+1)*(res[0]*res[1]), v_3 = x+y*res[0]+(z+1)*(res[0]*res[1]);
				int v_4 = x+(y+1)*res[0]+z*(res[0]*res[1]), v_5 = (x+1)+(y+1)*res[0]+z*(res[0]*res[1]);
				int v_6 = (x+1)+(y+1)*res[0]+(z+1)*(res[0]*res[1]), v_7 = x+(y+1)*res[0]+(z+1)*(res[0]*res[1]);

				// compute the mask
				int mask = 0;
				if(f_0 <= 0) mask |= 1;
				if(f_1 <= 0) mask |= 2;
				if(f_2 <= 0) mask |= 4;
				if(f_3 <= 0) mask |= 8;
				if(f_4 <= 0) mask |= 16;
				if(f_5 <= 0) mask |= 32;
				if(f_6 <= 0) mask |= 64;
				if(f_7 <= 0) mask |= 128;

				if(msEdgeTable[mask] == 0)
					continue;

				mc_edge symbolic_edges[] = { mc_edge(v_0,v_1), mc_edge(v_1,v_2), mc_edge(v_2,v_3), mc_edge(v_3,v_0),
													  mc_edge(v_4,v_5), mc_edge(v_5,v_6), mc_edge(v_6,v_7), mc_edge(v_7,v_4),
													  mc_edge(v_0,v_4), mc_edge(v_1,v_5), mc_edge(v_2,v_6), mc_edge(v_3,v_7) };
				Vector3 corner_pts[] = { Vector3(x_value,y_value,z_value), Vector3(next_x_value,y_value,z_value),
												 Vector3(next_x_value,y_value,next_z_value), Vector3(x_value, y_value, next_z_value),
												 Vector3(x_value,next_y_value,z_value), Vector3(next_x_value,next_y_value,z_value),
												 Vector3(next_x_value,next_y_value,next_z_value), Vector3(x_value, next_y_value, next_z_value) };

				Vector3 edges[12];
				if(msEdgeTable[mask] & 1)
					edges[0] = interpolEdge(f_0,f_1, corner_pts[0], corner_pts[1]);
				if(msEdgeTable[mask] & 2)
					edges[1] = interpolEdge(f_1,f_2, corner_pts[1], corner_pts[2]);
				if(msEdgeTable[mask] & 4)
					edges[2] = interpolEdge(f_2,f_3, corner_pts[2], corner_pts[3]);
				if(msEdgeTable[mask] & 8)
					edges[3] = interpolEdge(f_3,f_0, corner_pts[3], corner_pts[0]);
				if(msEdgeTable[mask] & 16)
					edges[4] = interpolEdge(f_4,f_5, corner_pts[4], corner_pts[5]);
				if(msEdgeTable[mask] & 32)
					edges[5] = interpolEdge(f_5,f_6, corner_pts[5], corner_pts[6]);
				if(msEdgeTable[mask] & 64)
					edges[6] = interpolEdge(f_6,f_7, corner_pts[6], corner_pts[7]);
				if(msEdgeTable[mask] & 128)
					edges[7] = interpolEdge(f_7,f_4, corner_pts[7], corner_pts[4]);
				if(msEdgeTable[mask] & 256)
					edges[8] = interpolEdge(f_0,f_4, corner_pts[0], corner_pts[4]);
				if (msEdgeTable[mask] & 512)
					edges[9] = interpolEdge(f_1,f_5, corner_pts[1], corner_pts[5]);
				if(msEdgeTable[mask] & 1024)
					edges[10] = interpolEdge(f_2,f_6, corner_pts[2], corner_pts[6]);
				if(msEdgeTable[mask] & 2048)
					edges[11] = interpolEdge(f_3,f_7, corner_pts[3], corner_pts[7]);

				for(int i = 0; msTriTable[mask][i]!=-1; i+=3)  {
					mc_edge e0 = symbolic_edges[msTriTable[mask][i]];
					mc_edge inv_e0 = mc_edge(e0.second, e0.first);
					mc_edge e1 = symbolic_edges[msTriTable[mask][i+1]];
					mc_edge inv_e1 = mc_edge(e1.second, e1.first);
					mc_edge e2 = symbolic_edges[msTriTable[mask][i+2]];
					mc_edge inv_e2 = mc_edge(e2.second, e2.first);

					Vector3 p0 = edges[msTriTable[mask][i]];
					Vector3 p1 = edges[msTriTable[mask][i+1]];
					Vector3 p2 = edges[msTriTable[mask][i+2]];

					OpenMesh::VertexHandle v0h, v1h, v2h;

					if(vertex_hash.find(e0) == vertex_hash.end())  {
						BasicTriMesh::Point new_pt; new_pt[0] = p0.x; new_pt[1] = p0.y; new_pt[2] = p0.z;
						v0h = mc_mesh->add_vertex(new_pt);
						vertex_hash[e0] = v0h.idx();
						vertex_hash[inv_e0] = v0h.idx();
					}
					else
						v0h = OpenMesh::VertexHandle(vertex_hash[e0]);

					if(vertex_hash.find(e1) == vertex_hash.end())  {
						BasicTriMesh::Point new_pt; new_pt[0] = p1.x; new_pt[1] = p1.y; new_pt[2] = p1.z;
						v1h = mc_mesh->add_vertex(new_pt);
						vertex_hash[e1] = v1h.idx();
						vertex_hash[inv_e1] = v1h.idx();
					}
					else
						v1h = OpenMesh::VertexHandle(vertex_hash[e1]);

					if(vertex_hash.find(e2) == vertex_hash.end())  {
						BasicTriMesh::Point new_pt; new_pt[0] = p2.x; new_pt[1] = p2.y; new_pt[2] = p2.z;
						v2h = mc_mesh->add_vertex(new_pt);
						vertex_hash[e2] = v2h.idx();
						vertex_hash[inv_e2] = v2h.idx();
					}
					else
						v2h = OpenMesh::VertexHandle(vertex_hash[e2]);

					BasicTriMesh::Point pt1 = mc_mesh->point(v0h);
					BasicTriMesh::Point pt2 = mc_mesh->point(v1h);
					BasicTriMesh::Point pt3 = mc_mesh->point(v2h);
					if(pt1 == pt2 || pt1 == pt3 || pt2 == pt3)
						cout << "duplicate vertex!!!" << endl;
					mc_mesh->add_face(v0h, v1h, v2h);
				}
			}
		}
	}

	// clean up scalar field
	for(int z = 0; z < res[2]; z++)  {
		for(int y = 0; y < res[1]; y++)
			delete [] scalar_field[z][y];
	}
	for(int z = 0; z < res[2]; z++)
			delete [] scalar_field[z];
	delete [] scalar_field;

	return mc_mesh;
}

Vector3 ImplicitFunction::interpolEdge(double _f1, double _f2, Vector3 _p1, Vector3 _p2)  {
	if(_f2 == 0 && _f1 == 0)
		cout << "on a vertex!" << endl;
	double alpha = (0 - _f1) / (_f2 - _f1);
	return _p1 + (_p2-_p1)*alpha;
}

//************************************************************************
// Static precalculated table
//************************************************************************
int ImplicitFunction::msEdgeTable[256]={
    0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
};

int ImplicitFunction::msTriTable[256][16] = {
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
    {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
    {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
    {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
    {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
    {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
    {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
    {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
    {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
    {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
    {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
    {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
    {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
    {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
    {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
    {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
    {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
    {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
    {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
    {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
    {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
    {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
    {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
    {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
    {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
    {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
    {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
    {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
    {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
    {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
    {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
    {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
    {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
    {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
    {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
    {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
    {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
    {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
    {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
    {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
    {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
    {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
    {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
    {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
    {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
    {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
    {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
    {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
    {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
    {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
    {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
    {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
    {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
    {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
    {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
    {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
    {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
    {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
    {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
    {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
    {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
    {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
    {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
    {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
    {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
    {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
    {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
    {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
    {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
    {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
    {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
    {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
    {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
    {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
    {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
    {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
    {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
    {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
    {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
    {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
    {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
    {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
    {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
    {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
    {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
    {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
    {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
    {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
    {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
    {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
    {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
    {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
    {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
    {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
    {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};
