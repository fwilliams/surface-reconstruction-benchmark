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

#include "UniformGrid.h"

UniformGrid::UniformGrid(int _xRes, int _yRes, int _zRes, Vector3 _ll, Vector3 _ur)  {
	x_res = _xRes;
	y_res = _yRes;
	z_res = _zRes;
	ll_bound = _ll;
	ur_bound = _ur;

	uniform_grid = new GridCell***[z_res];
	for(int z = 0; z < z_res; z++)  {
		uniform_grid[z] = new GridCell**[y_res];
		for(int y = 0; y < y_res; y++)  {
			uniform_grid[z][y] = new GridCell*[x_res];
			for(int x = 0; x < x_res; x++)
				uniform_grid[z][y][x] = new GridCell();
		}
	}
}

UniformGrid::~UniformGrid()  {
	for(int z = 0; z < z_res; z++)  {
		for(int y = 0; y < y_res; y++)  {
			for(int x = 0; x < x_res; x++)
				delete uniform_grid[z][y][x];
		}
	}
	for(int z = 0; z < z_res; z++)  {
		for(int y = 0; y < y_res; y++)
			delete [] uniform_grid[z][y];
	}
	for(int z = 0; z < z_res; z++)
		delete [] uniform_grid[z];
	delete [] uniform_grid;
}

void UniformGrid::map_to_grid(Vector3 _pt, int& x_grid, int& y_grid, int& z_grid)  {
	x_grid = x_res*((_pt.x-ll_bound.x) / (ur_bound.x-ll_bound.x));
	y_grid = y_res*((_pt.y-ll_bound.y) / (ur_bound.y-ll_bound.y));
	z_grid = z_res*((_pt.z-ll_bound.z) / (ur_bound.z-ll_bound.z));

	if(x_grid < 0)
		x_grid = 0;
	if(x_grid >= x_res)
		x_grid = x_res-1;
	if(y_grid < 0)
		y_grid = 0;
	if(y_grid >= y_res)
		y_grid = y_res-1;
	if(z_grid < 0)
		z_grid = 0;
	if(z_grid >= z_res)
		z_grid = z_res-1;
}

void UniformGrid::add_point(GridPoint _pt)  {
	int gx, gy, gz;
	this->map_to_grid(_pt.pt, gx, gy, gz);
	//xxxxxxxx
	if(gx < 0 || gy < 0 || gz < 0)
		return;
	uniform_grid[gz][gy][gx]->add_point(_pt);
}

bool UniformGrid::remove_point(Vector3 _pt)  {
	int gx, gy, gz;
	this->map_to_grid(_pt, gx, gy, gz);
	//xxxxxxxx
	if(gx < 0 || gy < 0 || gz < 0)
		return false;
	return uniform_grid[gz][gy][gx]->remove_point(_pt);
}

void UniformGrid::neighborhood_query(Vector3 _pt, double _radius, vector<GridPoint>* neighborhood)  {
	double sqd_radius = _radius*_radius;

	Vector3 pt = _pt;

	Vector3 diag_sphere = Vector3(_radius,_radius,_radius);
	Vector3 ll_sphere = pt-diag_sphere, ur_sphere = pt+diag_sphere;
	int ll_x, ll_y, ll_z, ur_x, ur_y, ur_z;
	this->map_to_grid(ll_sphere, ll_x, ll_y, ll_z);
	this->map_to_grid(ur_sphere, ur_x, ur_y, ur_z);

	//xxxxxxxxxx
	if(ll_x < 0 || ll_y < 0 || ll_z < 0)
		return;

	for(int z = ll_z; z <= ur_z; z++)  {
		for(int y = ll_y; y <= ur_y; y++)  {
			for(int x = ll_x; x <= ur_x; x++)  {
				GridCell* next_cell = uniform_grid[z][y][x];
				GridCellParams grid_cell(x, y, z, x_res, y_res, z_res, ll_bound, ur_bound);
				if(!next_cell->intersects_sphere(grid_cell, pt, _radius))
					continue;
				for(int n = 0; n < next_cell->size(); n++)  {
					GridPoint next_pt = next_cell->get_grid_point(n);
					if(next_pt.pt.sqdDist(pt) > sqd_radius)
						continue;
					neighborhood->push_back(next_pt);
				}
			}
		}
	}
}
