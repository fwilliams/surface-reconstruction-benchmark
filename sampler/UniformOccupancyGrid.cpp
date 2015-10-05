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

#include "UniformOccupancyGrid.h"

UniformOccupancyGrid::UniformOccupancyGrid(Vector3 _ll, Vector3 _ur, int _res)  {
	ll_bound = _ll;
	ur_bound = _ur;
	res = _res;
	binary_volume = new bool[_res*_res*_res];
	int res3 = _res*_res*_res;
	for(int i = 0; i < res3; i++)
		binary_volume[i] = false;

	double inv_res = 1.0 / _res;
	cell_size = (ur_bound-ll_bound)*inv_res;

	cell_eps = cell_size.length()*1e-10;
}

UniformOccupancyGrid::~UniformOccupancyGrid()  {
	delete binary_volume;
}

bool UniformOccupancyGrid::isOccupied(Vector3 _pt)  {
	int cx,cy,cz;
	this->coord_to_cell(_pt, cx, cy, cz);
	//cout << _pt << " lattice cell <"<<cx<<","<<cy<<","<<cz<<">" << endl;
	int occ_flag = binary_volume[this->coord_to_index(_pt)] ? 1 : 0;
	//cout << "checking for occupation with pt " << _pt << " cell: <"<<cx<<","<<cy<<","<<cz<<"> " << occ_flag << endl;
	return binary_volume[this->coord_to_index(_pt)];
}

void UniformOccupancyGrid::post_process()  {
	bool* temp_binary = new bool[res*res*res];

	// for every lattice point, if its occupied, then mark neighbors to be occupied
	for(int z = 0; z < res; z++)  {
		for(int y = 0; y < res; y++)  {
			for(int x = 0; x < res; x++)  {
				if(!binary_volume[this->index(x,y,z)])
					continue;

				for(int i = -1; i <= 1; i++)  {
					for(int j = -1; j <= 1; j++)  {
						for(int k = -1; k <= 1; k++)  {
							temp_binary[this->index(x+i,y+j,z+k)] = true;
						}
					}
				}

			}
		}
	}

	// swap
	for(int z = 0; z < res; z++)  {
		for(int y = 0; y < res; y++)  {
			for(int x = 0; x < res; x++)
				binary_volume[this->index(x,y,z)] = temp_binary[this->index(x,y,z)];
		}
	}

	delete temp_binary;
}

double UniformOccupancyGrid::grid_march(Vector3 _pt, Vector3 _dir)  {
	// get cell
	int cell_x, cell_y, cell_z;
	this->coord_to_cell(_pt, cell_x, cell_y, cell_z);

	// get cell bounds
	Vector3 ll, ur;
	this->bounds(cell_x,cell_y,cell_z, ll, ur);

	/*
	cout << ll << endl;
	cout << _pt << endl;
	cout << ur << endl;
	*/

	double border_x = _dir.x < 0 ? ll.x : ur.x;
	double border_y = _dir.y < 0 ? ll.y : ur.y;
	double border_z = _dir.z < 0 ? ll.z : ur.z;

	double t_x = isNumericallyZero(_dir.x) ? 1e10 : (border_x-_pt.x)/_dir.x;
	double t_y = isNumericallyZero(_dir.y) ? 1e10 : (border_y-_pt.y)/_dir.y;
	double t_z = isNumericallyZero(_dir.z) ? 1e10 : (border_z-_pt.z)/_dir.z;

	double t_val = 0;
	if(t_x < t_y && t_x < t_z)
		t_val = t_x;
	else if(t_y < t_z)
		t_val = t_y;
	else
		t_val = t_z;

	return t_val;
}
