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

#ifndef IMPLICITFUNCTION_H
#define IMPLICITFUNCTION_H

#include "Vector3.h"
#include "DenseEigenSystem.h"
#include "mesh_types.h"

#if HAS_TEEM
#include <teem/nrrd.h>
#endif

#include <iostream>
#include <vector>
#include <stdio.h>
using namespace std;
using namespace BLinAlg;

struct ImplicitCurvatures  {
	ImplicitCurvatures()  {}
	ImplicitCurvatures(double _minCurvature, double _maxCurvature)  {
		min_curvature = _minCurvature;
		max_curvature = _maxCurvature;
	}
	ImplicitCurvatures(double _minCurvature, double _maxCurvature, Vector3 _minPC, Vector3 _maxPC)  {
		min_curvature = _minCurvature;
		max_curvature = _maxCurvature;
		pc_min = _minPC;
		pc_max = _maxPC;
	}

	double min_curvature, max_curvature;
	Vector3 pc_min, pc_max;
};

typedef std::pair<int, int> mc_edge;

class ImplicitFunction  {
	public:
		virtual double function(Vector3 _pt) = 0;
		virtual Vector3 gradient(Vector3 _pt) = 0;
		virtual void functionAndGradient(Vector3 _pt, double &function, Vector3& gradient)  {
			cerr << "ME STUPID ME DON'T KNOW HOW TO COMPUTE FUNCTION & GRADIENT" << endl;
		}

		virtual DenseMatrix hessian(Vector3 _pt) = 0;
		virtual Vector3 normal(Vector3 _pt) = 0;

		virtual void bounds(Vector3& ll_bound, Vector3& ur_bound) = 0;

		virtual Vector3 projection(Vector3 _pt, double _eps, int _maxIterations, double _lambda);
		float spatial_derivative(float*** _volume, Vector3 _voxel, int _dim, float _step, int _res);

		virtual void write_to_volume(string _fileout, int _maxRes, double _epsBound);
		virtual void write_to_volume(string _fileout, int _minRes, int _maxRes, double _epsBound);
		virtual void write_kitten_kaboodle(string _base, int _maxRes, double _epsBound);
		BasicTriMesh* isosurface(int _minRes, int _maxRes, double _epsBound);

		virtual void generatePovScene(string _filename)  {
			cerr << "ME STUPID ME DON'T KNOW HOW TO GENERATE POV SCENCE" << endl;
		}

		virtual bool is_signed_distance_inverted()  { return false; }
		virtual double lipschitz_constant(int _res);

		void generateTangentPlane(Vector3 _pt, Vector3& vec1, Vector3& vec2);
		DenseMatrix weingarten_map(Vector3 _pt);
		ImplicitCurvatures compute_curvatures(Vector3 _pt);

		virtual void dump_to_file(string _filename)  {
			cerr << "ME STUPID ME DON'T KNOW HOW TO DUMP TO FILE" << endl;
		}

		virtual vector<Vector3> sample_surface(int _numPoints, double _adaptivity, bool _dumpPts)  {
			vector<Vector3> nuttin;
			cerr << "ME STUPID ME DON'T KNOW HOW TO SAMPLE SURFACE" << endl;
			return nuttin;
		}

		virtual void sample_surface(int _numPoints, double _adaptivity, string _outfile)  {
			FILE* sampler_out = fopen(_outfile.c_str(), "w");
			vector<Vector3> sample_pts = this->sample_surface(_numPoints, _adaptivity, false);
			for(unsigned v = 0; v < sample_pts.size(); v++)  {
				Vector3 pt = sample_pts[v];
				Vector3 norm = this->normal(pt);
				fprintf(sampler_out, "%.7f %.7f %.7f %.7f %.7f %.7f\n", pt.x,pt.y,pt.z, norm.x,norm.y,norm.z);
			}
			fclose(sampler_out);
		}

		bool inBound(Vector3 _pt)  {
			Vector3 min, max;
			this->bounds(min, max);
			return _pt.x > min.x && _pt.x < max.x && _pt.y > min.y && _pt.y < max.y && _pt.z > min.z && _pt.z < max.z;
		}

		bool rayBoundIntersection(Vector3 _rayPos, Vector3 _rayDir, double& t_near, double& t_far)  {
			Vector3 llB, urB;
			this->bounds(llB, urB);

			double t0 = 0.0;
			double t1 = 1e10;

			double direction[] = { _rayDir.x, _rayDir.y, _rayDir.z };
			double origin[] = { _rayPos.x, _rayPos.y, _rayPos.z };
			double min[] = { llB.x, llB.y, llB.z };
			double max[] = { urB.x, urB.y, urB.z };

			for(int i = 0; i < 3; i++)  {
				double invRayDir = 1.0 / direction[i];
				double near_t = (min[i] - origin[i]) * invRayDir;
				double far_t = (max[i] - origin[i]) * invRayDir;

				if(near_t > far_t)  {
					double temp = near_t;
					near_t = far_t;
					far_t = temp;
				}
				t0 = near_t > t0 ? near_t : t0;
				t1 = far_t < t1 ? far_t : t1;

				if(t0 > t1) 
					return false;
			}

			t_near = t0;
			t_far = t1;
			return true;
		}

	private:
		static int msEdgeTable[256];
		static int msTriTable[256][16];

		Vector3 interpolEdge(double _f1, double _f2, Vector3 _p1, Vector3 _p2);
};

#endif
