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

#ifndef UNIFORMOCCUPANCYGRID_H
#define UNIFORMOCCUPANCYGRID_H

#include "../modeling/Vector3.h"
#include <iostream>
using namespace std;

class UniformOccupancyGrid  {
	public:
		UniformOccupancyGrid(Vector3 _ll, Vector3 _ur, int _res);
		~UniformOccupancyGrid();

		int index(int _x, int _y, int _z)  {
			int proper_x = _x < 0 ? 0 : _x;
			proper_x = _x >= res ? res-1 : _x;
			int proper_y = _y < 0 ? 0 : _y;
			proper_y = _y >= res ? res-1 : _y;
			int proper_z = _z < 0 ? 0 : _z;
			proper_z = _z >= res ? res-1 : _z;
			return proper_x + res*(proper_y + res*proper_z);
		}

		int coord_to_index(Vector3 _pt)  {
			Vector3 normalized_pt = (_pt-ll_bound);
			int lattice_x = normalized_pt.x / cell_size.x;
			int lattice_y = normalized_pt.y / cell_size.y;
			int lattice_z = normalized_pt.z / cell_size.z;
			return index(lattice_x, lattice_y, lattice_z);
		}

		void coord_to_cell(Vector3 _pt, int &x, int &y, int &z)  {
			Vector3 normalized_pt = (_pt-ll_bound);
			int lattice_x = normalized_pt.x / cell_size.x;
			int lattice_y = normalized_pt.y / cell_size.y;
			int lattice_z = normalized_pt.z / cell_size.z;

			x = lattice_x < 0 ? 0 : lattice_x;
			x = lattice_x >= res ? res-1 : lattice_x;
			y = lattice_y < 0 ? 0 : lattice_y;
			y = lattice_y >= res ? res-1 : lattice_y;
			z = lattice_z < 0 ? 0 : lattice_z;
			z = lattice_z >= res ? res-1 : lattice_z;
		}

		void bounds(int _cellX, int _cellY, int _cellZ, Vector3& ll, Vector3& ur)  {
			double ll_x = ll_bound.x + cell_size.x*_cellX;
			double ll_y = ll_bound.y + cell_size.y*_cellY;
			double ll_z = ll_bound.z + cell_size.z*_cellZ;
			double ur_x = ll_bound.x + cell_size.x*(_cellX+1);
			double ur_y = ll_bound.y + cell_size.y*(_cellY+1);
			double ur_z = ll_bound.z + cell_size.z*(_cellZ+1);
			ll = Vector3(ll_x,ll_y,ll_z);
			ur = Vector3(ur_x,ur_y,ur_z);
		}

		void addPoint(Vector3 _pt)  {
			binary_volume[this->coord_to_index(_pt)] = true;
		}
		void post_process();

		bool isOccupied(Vector3 _pt);

		bool isNumericallyZero(double _val)  {
			double abs_val = _val < 0 ? -_val : _val;
			return abs_val < 1e-8;
		}

		double getCellEps()  { return cell_eps; }

		double grid_march(Vector3 _pt, Vector3 _dir);

	private:
		Vector3 ll_bound, ur_bound;
		Vector3 cell_size;
		int res;
		bool* binary_volume;

		double cell_eps;
};

#endif
