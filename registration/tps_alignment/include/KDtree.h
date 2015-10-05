#ifndef KDTREE_H
#define KDTREE_H
/*
Copyright 2007 Szymon Rusinkiewicz
Princeton University

KDtree.h
A K-D tree for points, with limited capabilities (find nearest point to 
a given point, or to a ray). 

-----

This file is part of tps_alignment.

tps_alignment is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.

tps_alignment is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <vector>

class KDtree {
private:
	mutable const float *cache;
	class Node;
	Node *root;
	void build(const float *ptlist, int n);

public:
	struct CompatFunc { virtual bool operator () (const float *p) const = 0; };

	KDtree(const float *ptlist, int n)
		{ build(ptlist, n); }
	template <class T> KDtree(std::vector<T> &v)
		{ build((const float *) &v[0], v.size()); }
	~KDtree();

	const float *closest_to_pt(const float *p,
				   float maxdist2,
				   const CompatFunc *iscompat = NULL) const;
	const float *closest_to_ray(const float *p, const float *dir,
				    float maxdist2,
				    const CompatFunc *iscompat = NULL) const;

};

#endif
