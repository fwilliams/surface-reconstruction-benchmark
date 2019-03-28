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

#ifndef UNIFORMGRID_H
#define UNIFORMGRID_H

#include "../modeling/Vector3.h"
#include <vector>
using namespace std;

struct GridPoint  {
	GridPoint() : pt(Vector3()) , ind(-1)  {}
	GridPoint(Vector3 _pt, int _ind) : pt(_pt) , ind(_ind)  {}
	Vector3 pt;
	int ind;
};

struct GridCellParams  {
	GridCellParams()  {}
	GridCellParams(int _x, int _y, int _z, int _resX, int _resY, int _resZ, Vector3 _ll, Vector3 _ur)  {
		x = _x;
		y = _y;
		z = _z;
		res_x = _resX;
		res_y = _resY;
		res_z = _resZ;
		llBound = _ll;
		urBound = _ur;
	}
	int x, y, z;
	int res_x, res_y, res_z;
	Vector3 llBound, urBound;
};

struct GridCell  {
	GridCell() {
		is_instantiated = false;
	}
	~GridCell()  {
		if(is_instantiated)  {
			all_points->clear();
			delete all_points;
		}
	}

	int size()  {
		if(!is_instantiated)
			return 0;
		return all_points->size();
	}

	void add_point(GridPoint _pt)   {
		if(!is_instantiated)  {
			all_points = new vector<GridPoint>();
			is_instantiated = true;
		}
		all_points->push_back(_pt);
	}

	GridPoint get_grid_point(int _i)  { return all_points->at(_i); }
	Vector3 get_point(int _i)  { return all_points->at(_i).pt; }

	void clear()  {
		if(is_instantiated)
			all_points->clear();
	}

	void compute_bounds(GridCellParams _cell, Vector3& ll, Vector3& ur)  {
		double z_alpha_start = (double)_cell.z/(double)_cell.res_z, z_alpha_end = (double)(_cell.z+1)/(double)_cell.res_z;
		double z_start = _cell.llBound.z*(1.0-z_alpha_start) + _cell.urBound.z*z_alpha_start;
		double z_end = _cell.llBound.z*(1.0-z_alpha_end) + _cell.urBound.z*z_alpha_end;

		double y_alpha_start = (double)_cell.y/(double)_cell.res_y, y_alpha_end = (double)(_cell.y+1)/(double)_cell.res_y;
		double y_start = _cell.llBound.y*(1.0-y_alpha_start) + _cell.urBound.y*y_alpha_start;
		double y_end = _cell.llBound.y*(1.0-y_alpha_end) + _cell.urBound.y*y_alpha_end;

		double x_alpha_start = (double)_cell.x/(double)_cell.res_x, x_alpha_end = (double)(_cell.x+1)/(double)_cell.res_x;
		double x_start = _cell.llBound.x*(1.0-x_alpha_start) + _cell.urBound.x*x_alpha_start;
		double x_end = _cell.llBound.x*(1.0-x_alpha_end) + _cell.urBound.x*x_alpha_end;

		ll = Vector3(x_start,y_start,z_start);
		ur = Vector3(x_end,y_end,z_end);
	}

	bool in_bounds(GridCellParams _cell, Vector3 _pt)  {
		Vector3 ll, ur;
		this->compute_bounds(_cell, ll, ur);
		return _pt.x >= ll.x && _pt.y >= ll.y && _pt.z >= ll.z && _pt.x <= ur.x && _pt.y <= ur.y && _pt.z <= ur.z;
	}

	bool remove_point(GridPoint _pt)  {
		if(!is_instantiated)
			return false;
		for(unsigned i = 0; i < all_points->size(); i++)  {
			if(_pt.pt == all_points->at(i).pt)  {
				all_points->erase(all_points->begin()+i,all_points->begin()+(i+1));
				return true;
			}
		}
		return false;
	}

	bool remove_point(Vector3 _pt)  {
		if(!is_instantiated)
			return false;
		for(unsigned i = 0; i < all_points->size(); i++)  {
			if(_pt == all_points->at(i).pt)  {
				all_points->erase(all_points->begin()+i,all_points->begin()+(i+1));
				return true;
			}
		}
		return false;
	}

	bool intersects_sphere(GridCellParams _cell, Vector3 _pt, double _r)  {
		Vector3 ll, ur;
		this->compute_bounds(_cell, ll, ur);
		double r2 = _r*_r;
		double bmin[] = { ll.x, ll.y, ll.z };
		double bmax[] = { ur.x, ur.y, ur.z };
		double c[] = { _pt.x, _pt.y, _pt.z };
		double dmin = 0, diff = 0;

		for(int i = 0; i < 3; i++)  {
			if(c[i] < bmin[i])  {
				diff = c[i]-bmin[i];
				dmin += diff*diff;
			}
			else if(c[i] > bmax[i])  {
				diff = c[i]-bmax[i];
				dmin += diff*diff;
			}
		}

		return dmin <= r2;
	}

	vector<GridPoint>* all_points;
	bool is_instantiated;
};

class UniformGrid  {
	public:
		UniformGrid(int _xRes, int _yRes, int _zRes, Vector3 _ll, Vector3 _ur);
		~UniformGrid();

		void map_to_grid(Vector3 _pt, int& x_grid, int& y_grid, int& z_grid);

		void add_point(GridPoint _pt);
		bool remove_point(Vector3 _pt);

		void neighborhood_query(Vector3 _pt, double _radius, vector<GridPoint>* neighborhood);

	private:
		int x_res, y_res, z_res;
		Vector3 ll_bound, ur_bound;
		GridCell**** uniform_grid;
};

#endif
