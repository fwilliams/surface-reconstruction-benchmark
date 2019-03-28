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

#include "../modeling/shape_loader.h"
#include "NPTSReader.h"

#include <iostream>
using namespace std;

int main(int argc, char** argv)   {
	if(argc != 4)  {
		cerr << "usage : " << argv[0] << " pc shape output_pc" << endl;
		return 1;
	}

	OrientedPointCloud pc;
	NPTSReader npts_reader(argv[1]);
	npts_reader.read_pc(&pc);
	ImplicitFunction* shape = ShapeLoader::load_shape(argv[2]);
	string out_pc_filename = argv[3];

	if(false)  {
		OrientedPointCloud oriented_pc;
		for(int i = 0; i < pc.size(); i++)  {
			Vector3 pt = pc.getPoint(i);
			Vector3 pt_normal = pc.getNormal(i);
			Vector3 shape_normal = shape->normal(pt);
			if(pt_normal.dotProduct(shape_normal) < 0)
				pt_normal.scale(-1);
			oriented_pc.addOrientedPoint(pt, pt_normal);
		}
		oriented_pc.write_to_file(out_pc_filename);
		return 1;
	}

	if(out_pc_filename.compare(argv[1]) == 0)  {
		cerr << "output same as input! don't overwrite that shit!" << endl;
		return 1;
	}

	double lambda = 0.6;
	int max_iterations = 10;

	OrientedPointCloud projected_pc;
	double ave_func = 0, ave_proj_func = 0;
	for(int i = 0; i < pc.size(); i++)  {
		Vector3 pt = pc.getPoint(i);
		double approx_distance = shape->function(pt) / (shape->gradient(pt)).length();
		double abs_approx = fabs(approx_distance);
		ave_func += abs_approx;

		double prev_dist = abs_approx;
		for(int j = 0; j < max_iterations; j++)  {
			double pt_function = shape->function(pt);
			Vector3 pt_gradient = shape->gradient(pt);
			Vector3 pt_normal = pt_gradient;
			pt_normal.normalize();
			pt_normal.scale(-(pt_function*lambda)/pt_gradient.length());
			pt = pt + pt_normal;

			double next_dist = pt_function / pt_gradient.length();
			double next_abs = fabs(next_dist);
			if(next_abs > prev_dist)
				cerr << "error in projection: " << next_abs << " : " << prev_dist << endl;
			prev_dist = next_abs;
		}

		double corrected_func = shape->function(pt);
		Vector3 corrected_gradient = shape->gradient(pt);
		Vector3 pt_normal = corrected_gradient;
		pt_normal.normalize();

		approx_distance = corrected_func / corrected_gradient.length();
		abs_approx = fabs(approx_distance);
		ave_proj_func += abs_approx;

		projected_pc.addOrientedPoint(pt, pt_normal);
	}
	cout << "ave func: " << (ave_func / (double)pc.size()) << endl;
	cout << "ave proj func: " << (ave_proj_func / (double)pc.size()) << endl;

	projected_pc.write_to_file(out_pc_filename);
}
