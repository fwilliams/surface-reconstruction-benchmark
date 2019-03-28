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

#include "XYZReader.h"
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "NPTSReader.h"
#include "SparseRangeScan.h"
#include "PointCloud.h"
#include "../modeling/shape_loader.h"

#include "mesh_types.h"

int main(int argc, char** argv)  {
	if(argc != 3)  {
		cerr << "usage: " << argv[0] << " pc_in pc_out" << endl;
		return 1;
	}

	int arg_num=1;
	string pc_in = argv[arg_num++];
	string pc_out = argv[arg_num++];

	bool invert_orientation = false;
	bool center_pc = false;
	bool orient_pc = true;
	bool orient_scan = false;

	if(invert_orientation)  {
		OrientedPointCloud pc;
		NPTSReader pc_reader(pc_in);
		pc_reader.read_pc(&pc);

		OrientedPointCloud inv_pc;
		for(int i = 0; i < pc.size(); i++)  {
			Vector3 pt = pc.getPoint(i);
			Vector3 inv_normal = pc.getNormal(i);
			inv_normal.scale(-1);
			inv_pc.addOrientedPoint(pt, inv_normal);
		}
		inv_pc.write_to_file(pc_out);
		return 0;
	}

	if(center_pc)  {
		PointCloud pc;
		XYZReader xyz_reader(pc_in);
		xyz_reader.read_pc(&pc);

		Vector3 min_pc(1e10,1e10,1e10), max_pc(-1e10,-1e10,-1e10);
		for(int s = 0; s < pc.size(); s++)  {
			Vector3 pt = pc.getPoint(s);
			min_pc.expand_min_bound(pt);
			max_pc.expand_max_bound(pt);
		}
		cout << "min : " << min_pc << " max : " << max_pc << endl;

		double rot_angle = acos(-1)/2.0;
		PointCloud centered_pc;
		Vector3 center = (max_pc+min_pc)*0.5;
		for(int s = 0; s < pc.size(); s++)  {
			Vector3 pt = pc.getPoint(s);
			Vector3 centered_pt = pt-center;
			centered_pt.y *= -1;
			centered_pt.x *= -1;

			double rotated_x = centered_pt.x*cos(rot_angle) + centered_pt.z*sin(rot_angle);
			double rotated_z = -centered_pt.x*sin(rot_angle) + centered_pt.z*cos(rot_angle);
			Vector3 rotated_pt = Vector3(rotated_x, centered_pt.y, rotated_z);

			centered_pc.addPoint(rotated_pt);
		}

		cout << "center: " << center << endl;
		Vector3 camera_pt(227.604, -58.3642, -55.6916);
		Vector3 centered_pt = camera_pt-center;
		centered_pt.y *= -1;
		centered_pt.x *= -1;

		double rotated_x = centered_pt.x*cos(rot_angle) + centered_pt.z*sin(rot_angle);
		double rotated_z = -centered_pt.x*sin(rot_angle) + centered_pt.z*cos(rot_angle);
		Vector3 rotated_pt = Vector3(rotated_x, centered_pt.y, rotated_z);
		cout << "new camera pt: " << rotated_pt << endl;
		//centered_pc.write_to_file(pc_out);
		return 0;
	}

	if(orient_pc)  {
		PointCloud pc;
		XYZReader xyz_reader(pc_in);
		xyz_reader.read_pc(&pc);

		OrientedPointCloud* oriented_pc = pc.orient_points(12);
		oriented_pc->write_to_file(pc_out);
		return 0;
	}

	if(orient_scan)  {
		SparseRangeScan range_scan(pc_in);
		PointCloud pc;
		for(int i = 0; i < range_scan.size(); i++)
			pc.addPoint(range_scan.getPt(i));

		OrientedPointCloud* oriented_pc = pc.orient_points(6);
		oriented_pc->write_to_file(pc_out);
		return 0;
	}
}
