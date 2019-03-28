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

#ifndef IMPLICITSPHERE_H
#define IMPLICITSPHERE_H

#include "ImplicitFunction.h"
#include <stdlib.h>
#include <fstream>
#include <sstream>

class ImplicitSphere : public ImplicitFunction {
	public:
		ImplicitSphere(double _radius)  {
			radius = _radius;
			center = Vector3(0,0,0);
		}
		ImplicitSphere(double _radius, Vector3 _center)  {
			radius = _radius;
			center = _center;
		}
		~ImplicitSphere();

		virtual double function(Vector3 _pt)  {
			double diff_x = _pt.x-center.x;
			double diff_y = _pt.y-center.y;
			double diff_z = _pt.z-center.z;
			return diff_x*diff_x + diff_y*diff_y + diff_z*diff_z - (radius*radius);
		}

		virtual Vector3 gradient(Vector3 _pt)  {
			return Vector3(2.0*(_pt.x-center.x), 2.0*(_pt.y-center.y), 2.0*(_pt.z-center.z));
		}

		virtual DenseMatrix hessian(Vector3 _pt)  {
			DenseMatrix hessian_matrix(3);
			hessian_matrix.setEntry(0,0, 2.0);
			hessian_matrix.setEntry(1,1, 2.0);
			hessian_matrix.setEntry(2,2, 2.0);
			return hessian_matrix;
		}

		virtual Vector3 normal(Vector3 _pt)  {
			Vector3 simple_normal(_pt.x-center.x, _pt.y-center.y, _pt.z-center.z);
			simple_normal.normalize();
			return simple_normal;
		}

		virtual void bounds(Vector3& ll_bound, Vector3& ur_bound)  {
			ll_bound = Vector3(-radius, -radius, -radius);
			ur_bound = Vector3(radius, radius, radius);
		}

		/**
		 * STOLEN FROM JOSIAH MANSON. SHEEEEIT
		 */
		void tokenize_line(vector<string>* tokens, const string& input, string sep)  {
			string comm;
			for (unsigned int i = 0; i < input.size(); i++)  {
				const char ch = input[i];
				bool added = false;
				for (unsigned int s = 0; s < sep.size(); s++)  {
					if (ch == sep[s])  {
						tokens->push_back(comm);
						comm = "";
						added = true;
						break;
					}
				}
				if (!added)
					comm += ch;
			}
			if (comm != "")
				tokens->push_back(comm);
		}

		vector<Vector3> sample_surface(int _numPoints, double _adaptivity, bool _dumpPts)  {
			string sample_out = "sphere_sampled";

			ostringstream radius_str, num_samples;
			radius_str << radius;
			num_samples << _numPoints;

			string exec_sampler = "./optimize-sphere " + radius_str.str() + " " + 
									num_samples.str() + " " + sample_out;
			int sampler_result = system(exec_sampler.c_str());

			vector<Vector3> all_points;
			string npts_filename = sample_out + ".npts";
			char line_raw[256];
			ifstream npts_file(npts_filename.c_str());
			if (npts_file.fail())
				return all_points;

			for (;;)  {
				npts_file.getline(line_raw, 256);
				if (npts_file.eof())
					break;

				string line(line_raw);
				vector<string> tokens;
				tokenize_line(&tokens, line, " \t");

				Vector3 new_position = center + Vector3(atof(tokens[0].c_str()), atof(tokens[1].c_str()), atof(tokens[2].c_str()));
				all_points.push_back(new_position);
			}
			npts_file.close();

			string exec_rm;
			if(_dumpPts)
				exec_rm = "rm " + sample_out + ".ptcl";
			else
				exec_rm = "rm " + sample_out + ".npts " + sample_out + ".ptcl";

			int rm_result = system(exec_rm.c_str());

			return all_points;
		}

	private:
		double radius;
		Vector3 center;
};

#endif
