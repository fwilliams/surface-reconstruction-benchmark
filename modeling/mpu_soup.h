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

#ifndef MPU_SOUP
#define MPU_SOUP

#include "collision.h"

#include "tri_integrator.h"
#include "DenseMatrix.h"
#include "DenseEigenSystem.h"

#include <ANN/ANN.h>
#include <stdio.h>

#define MAX_OCTREE_DEPTH 11

static int int_triples[3][3] = { {0, 1, 2} , {1, 0, 2} , {2, 0, 1} };
static int int_quadruples[12][4] = { {0, 1, 2, 3} , {0, 2, 1, 3}, {0, 3, 2, 1} ,
											 {1, 0, 2, 3} , {1, 2, 0, 3}, {1, 3, 2, 0} ,
											 {2, 0, 1, 3} , {2, 1, 0, 3}, {2, 3, 0, 1} ,
											 {3, 0, 1, 2} , {3, 1, 0, 2}, {3, 2, 1, 0} };

struct Int3  {
	Int3() : i1(0) , i2(1) , i3(2) {}
	Int3(int _i1, int _i2, int _i3) : i1(_i1) , i2(_i2) , i3(_i3) {}
	int i1, i2, i3;
};

struct AABB  {
	AABB() : ll(1e10,1e10,1e10) , ur(-1e10,-1e10,-1e10) {}
	AABB(Vector3 _ll, Vector3 _ur) : ll(_ll) , ur(_ur) {}

	void expand_bound(Vector3 _pt)  {
		ll.expand_min_bound(_pt);
		ur.expand_max_bound(_pt);
	}

	bool is_overlapping(AABB _bounds)  {
		Vector3 other_ll = _bounds.ll, other_ur = _bounds.ur;
		if(other_ll.x < ll.x || other_ll.x > ur.x)
			return false;
		if(other_ll.y < ll.y || other_ll.y > ur.y)
			return false;
		if(other_ll.z < ll.z || other_ll.z > ur.z)
			return false;
		return true;
	}

	double diagonal()  { return ll.length(ur); }
	Vector3 center()  { return (ll+ur)*0.5; }

	Vector3 ll, ur;
};

struct Sphere  {
	Sphere()  {}
	Sphere(Vector3 _center, double _radius) : c(_center) , r(_radius)  {}

	bool is_overlapping(AABB _bounds)  {
		double cen_x = c.x, cen_y = c.y, cen_z = c.z;

		if(c.x < _bounds.ll.x)
			cen_x = _bounds.ll.x;
		else if(c.x > _bounds.ur.x)
			cen_x = _bounds.ur.x;
		if(c.y < _bounds.ll.y)
			cen_y = _bounds.ll.y;
		else if(c.y > _bounds.ur.y)
			cen_y = _bounds.ur.y;
		if(c.z < _bounds.ll.z)
			cen_z = _bounds.ll.z;
		else if(c.z > _bounds.ur.z)
			cen_z = _bounds.ur.z;

		Vector3 closest_c(cen_x, cen_y, cen_z);
		return closest_c.sqdDist(c) < (r*r);
	}

	Vector3 c;
	double r;
};

struct LinearFunction  {
	LinearFunction()  {}
	LinearFunction(Vector3 _n, double _d) : n(_n) , d(_d)  {}

	double signed_dist(Vector3 _pt)  {
		return n.dotProduct(_pt) + d;
	}

	Vector3 n;
	double d;
};

struct MPUTri  {
	MPUTri()  {}
	MPUTri(Vector3 _v1, Vector3 _v2, Vector3 _v3, Vector3 _n1, Vector3 _n2, Vector3 _n3)  {
		v1 = _v1;
		v2 = _v2;
		v3 = _v3;
		n1 = _n1;
		n2 = _n2;
		n3 = _n3;

		bounds.expand_bound(v1);
		bounds.expand_bound(v2);
		bounds.expand_bound(v3);

		normal = (v2-v1).cross(v3-v1);
		area = 0.5*normal.length();
		normal.normalize();
	}

	Vector3 barycenter()  { return (v1+v2+v3)*(1.0/3.0); }

	Vector3 v1, v2, v3;
	Vector3 n1, n2, n3;
	AABB bounds;
	Vector3 normal;
	double area;
};

class MPUOctree  {
	public:
		MPUOctree(MPUTri** _tris, ANNkd_tree* _tree, int _numTris, int _minSamples, double _fitEps, double _covering);
		MPUOctree(MPUOctree* _parentNode, AABB _bounds);
		MPUOctree(FILE* _file, int _minSamples, double _fitEps, double _covering, AABB _bounds);
		~MPUOctree();

		void build_tree();

		bool within_covering(Vector3 _pt);
		bool is_leaf()  { return isLeaf; }
		bool is_root()  { return depth == 0; }

		void add_tri(int _ind)  { tri_inds.push_back(_ind); }
		unsigned size()  { return tri_inds.size(); }
		int get_tri_ind(int _t)  { return tri_inds[_t]; }
		MPUTri* get_tri(int _t)  { return tris[ tri_inds[_t] ]; }

		double interpolatory_weight(Vector3 _pt);
		double bspline2(Vector3 _pt);
		Vector3 gradient_bspline2(Vector3 _pt);
		void mpu_eval(Vector3 _pt, double& func_eval, double& weight);

		double function_eval(Vector3 _pt);
		Vector3 gradient_eval(Vector3 _pt);
		int number_of_fits()  { return linear_fits.size(); }

		double boolean_eval(int _boolean, Vector3 _pt,
				LinearFunction _fit1, LinearFunction _fit2, LinearFunction _fit3, LinearFunction _fit4);
		double boolean_eval(int _boolean, Vector3 _pt, LinearFunction _fit1, LinearFunction _fit2, LinearFunction _fit3);
		double boolean_eval(int _boolean, Vector3 _pt, LinearFunction _fit1, LinearFunction _fit2);

		int boolean_eval_int(int _boolean, Vector3 _pt,
				LinearFunction _fit1, LinearFunction _fit2, LinearFunction _fit3, LinearFunction _fit4);
		int boolean_eval_int(int _boolean, Vector3 _pt, LinearFunction _fit1, LinearFunction _fit2, LinearFunction _fit3);
		int boolean_eval_int(int _boolean, Vector3 _pt, LinearFunction _fit1, LinearFunction _fit2);

		void write_to_file(FILE* _file);

		MPUOctree** children;
		AABB bounds;
		Sphere support;

		// shared across all octree nodes
		MPUTri** tris;
		ANNkd_tree* kd_tree;
		int min_samples;
		double fit_eps, covering_lambda, smooth_eps;
		int num_steps;
		double* nonuniform_steps;

	private:
		MPUOctree* parent_node;
		vector<int> tri_inds;

		vector<LinearFunction> linear_fits;
		unsigned boolean_flags;
		int ord_ind;

		int depth;
		bool isLeaf;
		bool is_convex_shape;

		void gather_points();
		int partition_points(int* point_ids);
		void determine_boolean();
		void fit();
};

class MPUSoup : public ImplicitFunction  {
	public:
		MPUSoup(BasicTriMesh* _polygonSoup, int _minSamples, double _fitEps, double _covering);
		MPUSoup(string _filename);
		~MPUSoup();

		void write_to_file(string _filename);

		virtual double function(Vector3 _pt);
		virtual Vector3 gradient(Vector3 _pt);
		virtual Vector3 normal(Vector3 _pt);
		virtual DenseMatrix hessian(Vector3 _pt);

		void octree_func(Vector3 _pt, MPUOctree* _node, double& func_numer, double& func_denom);
		void octree_func(Vector3 _pt, MPUOctree* _node, vector<MPUOctree*>* all_nodes);
		void octree_gradient(Vector3 _pt, MPUOctree* _node, double& numer1, double& numer2,
																			 Vector3& grad1, Vector3& grad2);

		virtual void bounds(Vector3& ll_bound, Vector3& ur_bound);

	private:
		BasicTriMesh* polygon_soup;
		double eps, area_eps, my_pi;
		int num_interior, num_leaf;

		MPUTri** all_tris;
		int num_soup;
		MPUOctree* mpu;

		int min_samples;
		double fit_eps, covering;

		bool debug;

		vector<double> nonuniform_steps;

		ANNpointArray point_array;
		ANNkd_tree* kd_tree;

		void init();
};

class InvMPUSoup : public ImplicitFunction  {
	public:
		InvMPUSoup(BasicTriMesh* _polygonSoup, int _minSamples, double _fitEps, double _covering);
		InvMPUSoup(string _filename);
		~InvMPUSoup();

		void write_to_file(string _filename);

		virtual double function(Vector3 _pt);
		virtual Vector3 gradient(Vector3 _pt);
		virtual Vector3 normal(Vector3 _pt);
		virtual DenseMatrix hessian(Vector3 _pt);

		void octree_func(Vector3 _pt, MPUOctree* _node, double& func_numer, double& func_denom);
		void octree_gradient(Vector3 _pt, MPUOctree* _node, double& numer1, double& numer2,
																			 Vector3& grad1, Vector3& grad2);
		virtual void bounds(Vector3& ll_bound, Vector3& ur_bound);

	private:
		BasicTriMesh* polygon_soup;
		double eps, area_eps, my_pi;
		int num_interior, num_leaf;

		MPUTri** all_tris;
		int num_soup;
		MPUOctree* mpu;

		int min_samples;
		double fit_eps, covering;

		bool debug;

		vector<double> nonuniform_steps;

		ANNpointArray point_array;
		ANNkd_tree* kd_tree;

		void init();
};

#endif
