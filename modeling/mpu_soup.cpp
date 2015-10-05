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

#include "mpu_soup.h"

// ------ MPUOctree ------ //

MPUOctree::MPUOctree(MPUTri** _tris, ANNkd_tree* _tree, int _numTris, int _minSamples, double _fitEps, double _covering)  {
	tris = _tris;
	kd_tree = _tree;

	for(int t = 0; t < _numTris; t++)  {
		tri_inds.push_back(t);
		MPUTri* tri = tris[t];
		bounds.expand_bound(tri->bounds.ll);
		bounds.expand_bound(tri->bounds.ur);
	}

	Vector3 diagonal = bounds.ur-bounds.ll;
	bounds.ll = bounds.ll - diagonal*0.01;
	bounds.ur = bounds.ur + diagonal*0.01;

	fit_eps = _fitEps*diagonal.length();
	smooth_eps = 0;

	num_steps = 50;
	nonuniform_steps = new double[num_steps];
	int step_ind = 0;
	for(int n = num_steps-1; n >= 1; n--)
		nonuniform_steps[step_ind++] = pow((4.0/5.0), n);
	nonuniform_steps[step_ind++] = 1-1e-3;

	min_samples = _minSamples;
	covering_lambda = _covering;
	depth = 0;

	isLeaf = false;
	support = Sphere(bounds.center(), covering_lambda*bounds.diagonal());
	fit();
}

MPUOctree::MPUOctree(MPUOctree* _parentNode, AABB _bounds)  {
	parent_node = _parentNode;
	bounds = _bounds;

	tris = parent_node->tris;
	kd_tree = parent_node->kd_tree;
	fit_eps = parent_node->fit_eps;
	smooth_eps = parent_node->smooth_eps;
	num_steps = parent_node->num_steps;
	nonuniform_steps = parent_node->nonuniform_steps;

	min_samples = parent_node->min_samples;
	covering_lambda = parent_node->covering_lambda;
	depth = parent_node->depth+1;

	isLeaf = false;
	support = Sphere(bounds.center(), covering_lambda*bounds.diagonal());
	gather_points();
	fit();
}

MPUOctree::MPUOctree(FILE* _file, int _minSamples, double _fitEps, double _covering, AABB _bounds)  {
	min_samples = _minSamples;
	fit_eps = _fitEps;
	covering_lambda = _covering;
	bounds = _bounds;
	support = Sphere(bounds.center(), covering_lambda*bounds.diagonal());

	size_t num_read;
	int leaf_int;
	num_read = fread(&leaf_int, sizeof(int), 1, _file);
	int num_fits;
	num_read = fread(&num_fits, sizeof(int), 1, _file);

	isLeaf = leaf_int == 1;

	for(int f = 0; f < num_fits; f++)  {
		double n_fit[3];
		num_read = fread(n_fit, sizeof(double), 3, _file);
		double d_fit;
		num_read = fread(&d_fit, sizeof(double), 1, _file);
		linear_fits.push_back(LinearFunction(Vector3(n_fit[0],n_fit[1],n_fit[2]), d_fit));
	}
	if(num_fits > 1)
		num_read = fread(&boolean_flags, sizeof(unsigned), 1, _file);
	if(num_fits > 2)
		num_read = fread(&ord_ind, sizeof(int), 1, _file);

	if(isLeaf)
		return;

	double x_bounds[] = { bounds.ll.x , 0.5*(bounds.ll.x+bounds.ur.x) , bounds.ur.x };
	double y_bounds[] = { bounds.ll.y , 0.5*(bounds.ll.y+bounds.ur.y) , bounds.ur.y };
	double z_bounds[] = { bounds.ll.z , 0.5*(bounds.ll.z+bounds.ur.z) , bounds.ur.z };

	children = new MPUOctree*[8];

	int child_ind = 0;
	for(int x = 0; x < 2; x++)  {
		for(int y = 0; y < 2; y++)  {
			for(int z = 0; z < 2; z++)  {
				Vector3 child_ll(x_bounds[x], y_bounds[y], z_bounds[z]);
				Vector3 child_ur(x_bounds[x+1], y_bounds[y+1], z_bounds[z+1]);
				AABB child_bounds(child_ll, child_ur);
				children[child_ind++] = new MPUOctree(_file, min_samples, fit_eps, covering_lambda, child_bounds);
			}
		}
	}
}

MPUOctree::~MPUOctree()  {
	tri_inds.clear();
	if(!isLeaf)  {
		for(int i = 0; i < 8; i++)
			delete children[i];
		delete [] children;
	}
}

void MPUOctree::build_tree()  {
	if(isLeaf)
		return;

	double x_bounds[] = { bounds.ll.x , 0.5*(bounds.ll.x+bounds.ur.x) , bounds.ur.x };
	double y_bounds[] = { bounds.ll.y , 0.5*(bounds.ll.y+bounds.ur.y) , bounds.ur.y };
	double z_bounds[] = { bounds.ll.z , 0.5*(bounds.ll.z+bounds.ur.z) , bounds.ur.z };

	children = new MPUOctree*[8];

	int child_ind = 0;
	for(int x = 0; x < 2; x++)  {
		for(int y = 0; y < 2; y++)  {
			for(int z = 0; z < 2; z++)  {
				Vector3 child_ll(x_bounds[x], y_bounds[y], z_bounds[z]);
				Vector3 child_ur(x_bounds[x+1], y_bounds[y+1], z_bounds[z+1]);
				AABB child_bounds(child_ll, child_ur);
				children[child_ind++] = new MPUOctree(this, child_bounds);
			}
		}
	}

	for(int c = 0; c < 8; c++)
		children[c]->build_tree();
}

bool MPUOctree::within_covering(Vector3 _pt)  {
	return _pt.sqdDist(support.c) < (support.r*support.r);
}

void MPUOctree::gather_points()  {
	double tight_r = support.r / 1.1;
	int our_size = 0;
	int iteration = 1;

	BCollision::VECTOR DUMMY_Q;

	while(our_size < min_samples)  {
		tight_r *= 1.1;
		support.r = tight_r;
		tri_inds.clear();

		BCollision::BCollisionDetection collision;

		BCollision::BLOB sphere_blob;
		sphere_blob.Radius = support.r;
		sphere_blob.Centre[0] = support.c.x;
		sphere_blob.Centre[1] = support.c.y;
		sphere_blob.Centre[2] = support.c.z;

		int num_box_overlap = 0, num_tri_overlap = 0;

		for(int p = 0; p < parent_node->size(); p++)  {
			int tri_ind = parent_node->get_tri_ind(p);
			MPUTri* tri = tris[tri_ind];

			BCollision::TRIANGLE triangle;
			triangle.P1[0] = tri->v1.x; triangle.P1[1] = tri->v1.y; triangle.P1[2] = tri->v1.z;
			triangle.P2[0] = tri->v2.x; triangle.P2[1] = tri->v2.y; triangle.P2[2] = tri->v2.z;
			triangle.P3[0] = tri->v3.x; triangle.P3[1] = tri->v3.y; triangle.P3[2] = tri->v3.z;
			triangle.N[0] = tri->normal.x; triangle.N[1] = tri->normal.y; triangle.N[2] = tri->normal.z;

			bool is_tri_sphere_intersect = collision.SphereTriangleIntersectionTest(triangle, sphere_blob, DUMMY_Q) == 1;
			if(is_tri_sphere_intersect)
				tri_inds.push_back(tri_ind);
		}
		our_size = tri_inds.size();
		if(our_size == 0 && iteration == 1)
			isLeaf = true;
		if(our_size == parent_node->size() && our_size < min_samples)
			cout << "wahh? parent doesn't contain min samples?" << endl;
		iteration++;
	}

	if(depth == MAX_OCTREE_DEPTH)
		isLeaf = true;
}

int MPUOctree::partition_points(int* point_ids)  {
	// too many points to consider sharp feature fit
	//if(true)  {
	if(this->size() > 2*min_samples)  {
		for(int i = 0; i < this->size(); i++)
			point_ids[i] = 0;
		return 1;
	}

	// TODO: handle corners of degree 3 and 4
	int num_fits = 1;

	// edge detection
	double min_dot_1 = 1e10;
	for(int p = 0; p < this->size(); p++)  {
		MPUTri* p_tri = tris[ tri_inds[p] ];
		Vector3 p_normal = p_tri->normal;
		for(int q = (p+1); q < this->size(); q++)  {
			MPUTri* q_tri = tris[ tri_inds[q] ];
			Vector3 q_normal = q_tri->normal;
			double dot = p_normal.dotProduct(q_normal);
			min_dot_1 = dot < min_dot_1 ? dot : min_dot_1;
		}
	}

	// no sharp feature - one plane
	if(min_dot_1 >= 0.85)  {
		for(int i = 0; i < this->size(); i++)
			point_ids[i] = 0;
		return num_fits;
	}

	num_fits++;
	// otherwise, get the two normals who form a max angle deviation, and cluster normals into 2 groups
	double min_dot_2 = 1e10;
	Vector3 feature_n1(0,0,0), feature_n2;(0,0,0);
	for(int p = 0; p < this->size(); p++)  {
		MPUTri* p_tri = tris[ tri_inds[p] ];
		Vector3 p_normal = p_tri->normal;
		for(int q = (p+1); q < this->size(); q++)  {
			MPUTri* q_tri = tris[ tri_inds[q] ];
			Vector3 q_normal = q_tri->normal;
			double dot = p_normal.dotProduct(q_normal);
			if(dot < min_dot_2)  {
				feature_n1 = p_normal;
				feature_n2 = q_normal;
				min_dot_2 = dot;
			}
		}
	}

	// cluster points into 2
	for(int i = 0; i < this->size(); i++)  {
		MPUTri* tri = tris[ tri_inds[i] ];
		Vector3 normal = tri->normal;
		if(normal.dotProduct(feature_n1) > normal.dotProduct(feature_n2))
			point_ids[i] = 0;
		else
			point_ids[i] = 1;
	}

	// construct normal orthogonal to two planes - if there exists a normal which agrees with it (>0.65), then we may have a corner?
	Vector3 feature_n3 = feature_n1.cross(feature_n2);
	double max_dot_3 = -1e10;
	for(int p = 0; p < this->size(); p++)  {
		MPUTri* p_tri = tris[ tri_inds[p] ];
		Vector3 p_normal = p_tri->normal;
		double local_dot = fabs(p_normal.dotProduct(feature_n3));
		max_dot_3 = local_dot > max_dot_3 ? local_dot : max_dot_3;
	}

	// if the max dot is sufficiently small (i.e. points exist along edge), then we have an edge, so return as is
	if(max_dot_3 < 0.65)
		return num_fits;

	// otherwise, make a new cluster for points which agree better with n3
	int num_3corner = 0;
	for(int i = 0; i < this->size(); i++)  {
		MPUTri* tri = tris[ tri_inds[i] ];
		Vector3 normal = tri->normal;
		double dot1 = fabs(normal.dotProduct(feature_n1));
		double dot2 = fabs(normal.dotProduct(feature_n2));
		double dot3 = fabs(normal.dotProduct(feature_n3));
		if(dot3 > dot1 && dot3 > dot2)  {
			point_ids[i] = 2;
			num_3corner++;
		}
	}

	if(num_3corner > 0)
		num_fits++;

	//cout << "num 3 corner: " << num_3corner << endl;

	// FINALLY, check for a fucking corner of degree 4 motherfucker - same test as edge test, applied to all groups!
	bool is_corner_4 = false;
	int group_split = -1;
	for(int i = 0; i < 3; i++)  {
		double min_dot_4 = 1e10;
		for(int p = 0; p < this->size(); p++)  {
			if(point_ids[p] != i)
				continue;
			MPUTri* p_tri = tris[ tri_inds[p] ];
			Vector3 p_normal = p_tri->normal;
			for(int q = (p+1); q < this->size(); q++)  {
				if(point_ids[q] != i)
					continue;
				MPUTri* q_tri = tris[ tri_inds[q] ];
				Vector3 q_normal = q_tri->normal;
				double dot = p_normal.dotProduct(q_normal);
				min_dot_4 = dot < min_dot_4 ? dot : min_dot_4;
			}
		}

		if(min_dot_4 < 0.85)  {
			is_corner_4 = true;
			group_split = i;
			break;
		}
	}

	// if no corner 4, then return
	if(!is_corner_4)
		return num_fits;

	// otherwise, break up this group, return
	double min_dot_5 = 1e10;
	Vector3 feature_n1_prime(0,0,0), feature_n2_prime;(0,0,0);
	for(int p = 0; p < this->size(); p++)  {
		if(point_ids[p] != group_split)
			continue;
		MPUTri* p_tri = tris[ tri_inds[p] ];
		Vector3 p_normal = p_tri->normal;
		for(int q = (p+1); q < this->size(); q++)  {
			if(point_ids[q] != group_split)
				continue;
			MPUTri* q_tri = tris[ tri_inds[q] ];
			Vector3 q_normal = q_tri->normal;
			double dot = p_normal.dotProduct(q_normal);
			if(dot < min_dot_5)  {
				feature_n1_prime = p_normal;
				feature_n2_prime = q_normal;
				min_dot_5 = dot;
			}
		}
	}
	/*
	for(int i = 0; i < this->size(); i++)  {
		MPUTri* tri = tris[ tri_inds[i] ];
		Vector3 normal = tri->normal;
		cout << "normal["<<i<<"] : " << normal << endl;
	}
	cout << "feature n1: " << feature_n1 << " : feature n2: " << feature_n2 << " ; feature n3: " << feature_n3 << " ; num corner3 : " << num_3corner << endl;
	cout << "prime feature n1: " << feature_n1_prime << " : feature n2: " << feature_n2_prime << endl;
	*/

	// cluster points into 2
	//cout << "group split " << group_split << " : " << feature_n1_prime << " : " << feature_n2_prime << endl;
	for(int i = 0; i < this->size(); i++)  {
		if(point_ids[i] != group_split)
			continue;
		MPUTri* tri = tris[ tri_inds[i] ];
		Vector3 normal = tri->normal;
		if(normal.dotProduct(feature_n1_prime) > normal.dotProduct(feature_n2_prime))
			point_ids[i] = group_split;
		else
			point_ids[i] = num_fits;
	}

	num_fits++;
	return num_fits;
}

void MPUOctree::determine_boolean()  {
	switch(linear_fits.size())  {
		case 1:
			boolean_flags = 0;
			return;
		case 2:
			{
				double max_error = 1e10;
				for(int b = 0; b < 2; b++)  {
					double local_max = -1e10;

					for(int i = 0; i < this->size(); i++)  {
						MPUTri* tri = tris[ tri_inds[i] ];
						Vector3 normal = tri->normal;
						Vector3 bary = tri->barycenter();
						double error = this->boolean_eval(b, bary, linear_fits[0], linear_fits[1]);
						double abs_error = fabs(error);
						if(abs_error > local_max)
							local_max = abs_error;
					}
					if(local_max < max_error)  {
						max_error = local_max;
						boolean_flags = b;
					}
				}
			}
			return;

		case 3:
			{
				double max_error = 1e10;
				for(unsigned o = 0; o < 3; o++)  {
					int* ord = int_triples[o];
					for(int b1 = 0; b1 < 2; b1++)  {
						for(int b2 = 0; b2 < 2; b2++)  {
							int local_bool = b1 | (2*b2);

							double local_max = -1e10;
							for(int i = 0; i < this->size(); i++)  {
								MPUTri* tri = tris[ tri_inds[i] ];
								Vector3 normal = tri->normal;
								Vector3 bary = tri->barycenter();
								double error = this->boolean_eval(local_bool, bary, linear_fits[ ord[0] ], linear_fits[ ord[1] ], linear_fits[ ord[2] ]);
								double abs_error = fabs(error);
								if(abs_error > local_max)
									local_max = abs_error;
							}
							if(local_max < max_error)  {
								max_error = local_max;
								boolean_flags = local_bool;
								ord_ind = o;
							}
						}
					}
				}

				int* best_ord = int_triples[ord_ind];
				//cout << "best ordering[3] : " << best_ord[0] << " : " << best_ord[1] << " : " << best_ord[2] << " : " << max_error << endl;
			}
			return;

		case 4:
			{
				double max_error = 1e10;
				for(unsigned o = 0; o < 12; o++)  {
					int* ord = int_quadruples[o];
					for(int b1 = 0; b1 < 2; b1++)  {
						for(int b2 = 0; b2 < 2; b2++)  {
							for(int b3 = 0; b3 < 2; b3++)  {
								int local_bool = b1 | (2*b2) | (4*b3);

								double local_max = -1e10;
								for(int i = 0; i < this->size(); i++)  {
									MPUTri* tri = tris[ tri_inds[i] ];
									Vector3 normal = tri->normal;
									Vector3 bary = tri->barycenter();
									double error = this->boolean_eval(local_bool, bary, linear_fits[ ord[0] ], linear_fits[ ord[1] ],
																					linear_fits[ ord[2] ], linear_fits[ ord[3] ]);
									double abs_error = fabs(error);
									if(abs_error > local_max)
										local_max = abs_error;
								}
								if(local_max < max_error)  {
									max_error = local_max;
									boolean_flags = local_bool;
									ord_ind = o;
								}
							}
						}
					}
				}

				int* best_ord = int_quadruples[ord_ind];
				//cout << "best ordering[4] : " << best_ord[0] << " : " << best_ord[1] << " : " << best_ord[2] << " : " << best_ord[3] << " : " << max_error << endl;
			}
			return;
	}
}

double MPUOctree::boolean_eval(int _boolean, Vector3 _pt,
		LinearFunction _fit1, LinearFunction _fit2, LinearFunction _fit3, LinearFunction _fit4)  {
	double func1 = _fit1.signed_dist(_pt), func2 = this->boolean_eval(_boolean, _pt, _fit2, _fit3, _fit4);
	if((_boolean&4) == 0)
		return func1 < func2 ? func1 : func2;
	else if((_boolean&4) == 4)
		return func1 > func2 ? func1 : func2;

	// else
	cerr << "_booolean ? " << _boolean << endl;
	return 0;
}

double MPUOctree::boolean_eval(int _boolean, Vector3 _pt,
		LinearFunction _fit1, LinearFunction _fit2, LinearFunction _fit3)  {
	double func1 = _fit1.signed_dist(_pt), func2 = this->boolean_eval(_boolean, _pt, _fit2, _fit3);
	if((_boolean&2) == 0)
		return func1 < func2 ? func1 : func2;
	else if((_boolean&2) == 2)
		return func1 > func2 ? func1 : func2;

	// else
	cerr << "_booolean ? " << _boolean << endl;
	return 0;
}

double MPUOctree::boolean_eval(int _boolean, Vector3 _pt, LinearFunction _fit1, LinearFunction _fit2)  {
	double func1 = _fit1.signed_dist(_pt), func2 = _fit2.signed_dist(_pt);
	if((_boolean&1) == 0)
		return func1 < func2 ? func1 : func2;
	else if((_boolean&1) == 1)
		return func1 > func2 ? func1 : func2;

	// else
	cerr << "_booolean ? " << _boolean << endl;
	return 0;
}

int MPUOctree::boolean_eval_int(int _boolean, Vector3 _pt,
		LinearFunction _fit1, LinearFunction _fit2, LinearFunction _fit3, LinearFunction _fit4)  {
	double func1 = _fit1.signed_dist(_pt), func2 = this->boolean_eval(_boolean, _pt, _fit2, _fit3, _fit4);
	int local_ind = this->boolean_eval_int(_boolean, _pt, _fit2, _fit3, _fit4);
	if((_boolean&4) == 0)
		return func1 < func2 ? 0 : local_ind+1;
	else if((_boolean&4) == 4)
		return func1 > func2 ? 0 : local_ind+1;

	// else
	cerr << "_booolean ? " << _boolean << endl;
	return 0;
}

int MPUOctree::boolean_eval_int(int _boolean, Vector3 _pt, LinearFunction _fit1, LinearFunction _fit2, LinearFunction _fit3)  {
	double func1 = _fit1.signed_dist(_pt), func2 = this->boolean_eval(_boolean, _pt, _fit2, _fit3);
	int local_ind = this->boolean_eval_int(_boolean, _pt, _fit2, _fit3);
	if((_boolean&2) == 0)
		return func1 < func2 ? 0 : local_ind+1;
	else if((_boolean&2) == 2)
		return func1 > func2 ? 0 : local_ind+1;

	// else
	cerr << "_booolean ? " << _boolean << endl;
	return 0;
}

int MPUOctree::boolean_eval_int(int _boolean, Vector3 _pt, LinearFunction _fit1, LinearFunction _fit2)  {
	double func1 = _fit1.signed_dist(_pt), func2 = _fit2.signed_dist(_pt);
	if((_boolean&1) == 0)
		return func1 < func2 ? 0 : 1;
	else if((_boolean&1) == 1)
		return func1 > func2 ? 0 : 1;

	// else
	cerr << "_booolean ? " << _boolean << endl;
	return 0;
}

void MPUOctree::fit()  {
	int num_approximations = 0;

	int* point_ids = new int[this->size()];
	int num_fits = this->partition_points(point_ids);

	for(int i = 0; i < num_fits; i++)  {
		Vector3 normal_fit(0,0,0);
		Vector3 center_of_mass(0,0,0);
		double total_weight = 0;

		for(int t = 0; t < this->size(); t++)  {
			if(point_ids[t] != i)
				continue;

			MPUTri* tri = tris[ tri_inds[t] ];
			Vector3 bary = tri->barycenter();
			Vector3 normal = tri->normal;
			double area = tri->area;
			//if(num_fits == 3 && depth == (MAX_OCTREE_DEPTH-2))
			//	cout << "normal["<<i<<"] : " << normal << " ; bary["<<i<<"] : " << bary << endl;

			TriIntegrator integrals(support.c, tri->v1, tri->v2, tri->v3, smooth_eps, 1.75);
			integrals.set_steps(num_steps, nonuniform_steps);

			if(integrals.useApproximation())
				num_approximations++;

			double constant = integrals.constant_integration();

			TriFunction x_coord(tri->v1.x,tri->v2.x,tri->v3.x);
			TriFunction y_coord(tri->v1.y,tri->v2.y,tri->v3.y);
			TriFunction z_coord(tri->v1.z,tri->v2.z,tri->v3.z);
			double x_linear = integrals.linear_integration(x_coord);
			double y_linear = integrals.linear_integration(y_coord);
			double z_linear = integrals.linear_integration(z_coord);
			normal_fit = normal_fit + normal*constant;

			center_of_mass = center_of_mass + Vector3(x_linear,y_linear,z_linear);
			total_weight += constant;
		}
		if(total_weight == 0)
			cout << "no tris for " << i << " / " << num_fits << endl;

		normal_fit.normalize();
		center_of_mass.scale(1.0/total_weight);
		//if(num_fits == 3 && depth == (MAX_OCTREE_DEPTH-2))
		//	cout << "com["<<i<<"] : " << center_of_mass << endl;
		double d_fit = -normal_fit.dotProduct(center_of_mass);

		linear_fits.push_back(LinearFunction(normal_fit, d_fit));
	}
	this->determine_boolean();

	delete [] point_ids;

	// drop out if we are already known to be a leaf
	if(isLeaf)  {
		//tri_inds.clear();
		return;
	}

	double max_error = -1e10;
	for(int t = 0; t < this->size(); t++)  {
		MPUTri* tri = tris[ tri_inds[t] ];
		Vector3 bary = tri->barycenter();
		double func_approx = this->function_eval(bary);
		double abs_func = fabs(func_approx);
		max_error = abs_func > max_error ? abs_func : max_error;
	}

	isLeaf = (max_error < fit_eps);
	if(!isLeaf)  {
		if(this->size() == 1)
			cout << " !------> ["<<depth<<"] " << this->size() << " failed: " << max_error << " " << num_approximations << "/" << this->size() << endl;
		else if(depth == (MAX_OCTREE_DEPTH-1))
			cout << "["<<depth<<"] " << this->size() << " failed: " << max_error << " " << num_approximations << "/" << this->size() << " : " << linear_fits.size() << endl;
		return;
	}
	else  {
		/*
		if(depth == (MAX_OCTREE_DEPTH-1) && this->size() == 1)
			cout << "["<<depth<<"] " << this->size() << " good error: " << max_error << " " << num_approximations << "/" << this->size() << endl;
			*/
	}

	//tri_inds.clear();
}

double MPUOctree::function_eval(Vector3 _pt)  {
	double funky = 0;
	switch(linear_fits.size())  {
		case 1:
			funky = linear_fits[0].signed_dist(_pt);
			break;
		case 2:
			funky = this->boolean_eval(boolean_flags, _pt, linear_fits[0], linear_fits[1]);
			break;
		case 3:
			{
				int* ord = int_triples[ord_ind];
				funky = this->boolean_eval(boolean_flags, _pt, linear_fits[ ord[0] ], linear_fits[ ord[1] ], linear_fits[ ord[2] ]);
			}
			break;
		case 4:
			{
				int* ord = int_quadruples[ord_ind];
				funky = this->boolean_eval(boolean_flags, _pt, linear_fits[ ord[0] ], linear_fits[ ord[1] ],
																				linear_fits[ ord[2] ], linear_fits[ ord[3] ]);
			}
			break;
	}
	return funky;
}

double MPUOctree::interpolatory_weight(Vector3 _pt)  {
	double b_d = _pt.length(support.c);
	double b_r = support.r;
	if(b_d > b_r)
		return 0;
	double singular_weight = (b_r-b_d)/(b_r*b_d);
	return singular_weight*singular_weight;
}

double MPUOctree::bspline2(Vector3 _pt)  {
	double b_d = _pt.length(support.c);
	double b_r = support.r;
	if(b_d > b_r)
		return 0;
	else  {
		b_d = 1.5*(b_d/b_r);
		if(b_d < 0.5)
			return (-b_d*b_d + 0.75);
		else
			return (0.5*(1.5-b_d)*(1.5-b_d));
	}
}

Vector3 MPUOctree::gradient_bspline2(Vector3 _pt)  {
	double b_d = _pt.length(support.c);
	double b_r = support.r;
	if(b_d > b_r)
		return Vector3(0,0,0);
	else  {
		double ref_b_d = 1.5*(b_d/b_r);
		if(ref_b_d < 0.5)  {
			Vector3 gradient_vec = _pt-support.c;
			double gradient_weight = -(2.0*2.25)/(b_r*b_r);
			gradient_vec.scale(gradient_weight);
			return gradient_vec;
		}
		else  {
			Vector3 gradient_vec = _pt-support.c;
			double gradient_weight = (2.25/(b_r*b_r)) - (2.25/(b_r*b_d));
			gradient_vec.scale(gradient_weight);
			return gradient_vec;
		}
	}
}

void MPUOctree::mpu_eval(Vector3 _pt, double& func_eval, double& weight)  {
	weight = this->bspline2(_pt);
	func_eval = this->function_eval(_pt);
}

Vector3 MPUOctree::gradient_eval(Vector3 _pt)  {
	Vector3 funky;
	switch(linear_fits.size())  {
		case 1:
			funky = linear_fits[0].n;
			break;
		case 2:
			{
				int csg_ind = this->boolean_eval_int(boolean_flags, _pt, linear_fits[0], linear_fits[1]);
				funky = linear_fits[csg_ind].n;
			}
			break;
		case 3:
			{
				int* ord = int_triples[ord_ind];
				int csg_ind = this->boolean_eval_int(boolean_flags, _pt, linear_fits[ ord[0] ], linear_fits[ ord[1] ], linear_fits[ ord[2] ]);
				funky = linear_fits[ ord[csg_ind] ].n;
			}
			break;
		case 4:
			{
				int* ord = int_quadruples[ord_ind];
				int csg_ind = this->boolean_eval_int(boolean_flags, _pt, linear_fits[ ord[0] ], linear_fits[ ord[1] ],
																							linear_fits[ ord[2] ], linear_fits[ ord[3] ]);
				funky = linear_fits[ ord[csg_ind] ].n;
			}
			break;
	}
	return funky;
}

void MPUOctree::write_to_file(FILE* _file)  {
	// write out root node info
	if(depth == 0)  {
		double full_bounds[] = { bounds.ll.x, bounds.ll.y, bounds.ll.z, bounds.ur.x, bounds.ur.y, bounds.ur.z };
		fwrite(full_bounds, sizeof(double), 6, _file);
		fwrite(&min_samples, sizeof(int), 1, _file);
		fwrite(&fit_eps, sizeof(double), 1, _file);
		fwrite(&covering_lambda, sizeof(double), 1, _file);
	}

	// and then write out the rest
	int leaf_int = isLeaf ? 1 : 0;
	fwrite(&leaf_int, sizeof(int), 1, _file);
	int num_fits = linear_fits.size();
	fwrite(&num_fits, sizeof(int), 1, _file);
	for(int i = 0; i < num_fits; i++)  {
		Vector3 n_fit = linear_fits[i].n;
		double full_normal[] = { n_fit.x, n_fit.y, n_fit.z };
		fwrite(full_normal, sizeof(double), 3, _file);
		double d_fit = linear_fits[i].d;
		fwrite(&d_fit, sizeof(double), 1, _file);
	}
	if(num_fits > 1)
		fwrite(&boolean_flags, sizeof(unsigned), 1, _file);
	if(num_fits > 2)
		fwrite(&ord_ind, sizeof(int), 1, _file);

	// traverse the octree for writing, if we are not a leaf
	if(isLeaf)
		return;
	for(int i = 0; i < 8; i++)
		children[i]->write_to_file(_file);
}

// ------ MPUOctree ------ //

// ------ MPUSoup ------ //

MPUSoup::MPUSoup(BasicTriMesh* _polygonSoup, int _minSamples, double _fitEps, double _covering)  {
	polygon_soup = _polygonSoup;
	min_samples = _minSamples;
	fit_eps = _fitEps;
	covering = _covering;
	my_pi = acos(-1);
	init();
}

MPUSoup::MPUSoup(string _filename)  {
	FILE* file_in = fopen(_filename.c_str(), "rb");
	size_t num_read;

	double full_bounds[6];
	num_read = fread(full_bounds, sizeof(double), 6, file_in);
	AABB mpu_bounds(Vector3(full_bounds[0],full_bounds[1],full_bounds[2]), Vector3(full_bounds[3],full_bounds[4],full_bounds[5]));
	num_read = fread(&min_samples, sizeof(int), 1, file_in);
	num_read = fread(&fit_eps, sizeof(double), 1, file_in);
	num_read = fread(&covering, sizeof(double), 1, file_in);

	mpu = new MPUOctree(file_in, min_samples, fit_eps, covering, mpu_bounds);

	fclose(file_in);
}

MPUSoup::~MPUSoup()  {
}

void MPUSoup::write_to_file(string _filename)  {
	FILE* file_out = fopen(_filename.c_str(), "wb");
	mpu->write_to_file(file_out);
	fclose(file_out);
}

double MPUSoup::function(Vector3 _pt)  {
	num_interior = num_leaf = 0;
	debug = false;

	vector<MPUOctree*> octree_nodes;
	this->octree_func(_pt, mpu, &octree_nodes);

	double numer = 0, denom = 0;
	double weight, function_evaluation;
	for(unsigned i = 0; i < octree_nodes.size(); i++)  {
		MPUOctree* node = octree_nodes[i];
		node->mpu_eval(_pt, function_evaluation, weight);
		numer += weight*function_evaluation;
		denom += weight;
	}

	return numer / denom;
}

void MPUSoup::octree_func(Vector3 _pt, MPUOctree* _node, vector<MPUOctree*>* all_nodes)  {
	if(!_node->within_covering(_pt))
		return;

	// if we hit a leaf, then add its contribution
	if(_node->is_leaf())  {
		all_nodes->push_back(_node);
		return;
	}

	// otherwise, we are an interior node and within covering -> continue traversal
	for(int c = 0; c < 8; c++)
		this->octree_func(_pt, _node->children[c], all_nodes);
}

/*
double MPUSoup::function(Vector3 _pt)  {
	num_interior = num_leaf = 0;
	debug = false;

	double numer, denom;
	this->octree_func(_pt, mpu, numer, denom);

	if(denom == 0)  {
		cout << "zero denom! " << _pt << " bounds : " << mpu->bounds.ll << " , " << mpu->bounds.ur << endl;
		debug = true;
		this->octree_func(_pt, mpu, numer, denom);
	}
	else  {
		//cout << "func: " << (numer/denom) << endl;
	}

	return numer / denom;
}
*/

void MPUSoup::octree_func(Vector3 _pt, MPUOctree* _node, double& func_numer, double& func_denom)  {
	func_numer = func_denom = 0;
	if(!_node->within_covering(_pt))
		return;

	// if we hit a leaf, then add its contribution
	if(_node->is_leaf())  {
		double weight, function_evaluation;
		_node->mpu_eval(_pt, function_evaluation, weight);
		func_numer = weight*function_evaluation;
		func_denom = weight;
		if(debug)
			cout << "function : " << function_evaluation << " weight : " << weight << endl;
		return;
	}

	// otherwise, we are an interior node and within covering -> continue traversal
	for(int c = 0; c < 8; c++)  {
		double child_numer = 0, child_denom = 0;
		this->octree_func(_pt, _node->children[c], child_numer, child_denom);
		func_numer += child_numer;
		func_denom += child_denom;
	}
}

Vector3 MPUSoup::gradient(Vector3 _pt)  {
	double numer1 = 0, numer2 = 0;
	Vector3 grad1(0,0,0), grad2(0,0,0);
	this->octree_gradient(_pt, mpu, numer1, numer2, grad1, grad2);

	grad1.scale(numer1);
	grad2.scale(-numer2);
	double weight_scale = 1.0/(numer1*numer1);

	Vector3 final_gradient = grad1+grad2;
	final_gradient.scale(weight_scale);

	return final_gradient;
}

Vector3 MPUSoup::normal(Vector3 _pt)  {
	Vector3 normal_vec = this->gradient(_pt);
	normal_vec.normalize();
	return normal_vec;
}

void MPUSoup::octree_gradient(Vector3 _pt, MPUOctree* _node, double& numer1, double& numer2,
																				 Vector3& grad1, Vector3& grad2)  {
	numer1 = numer2 = 0;
	grad1 = Vector3(0,0,0);
	grad2 = Vector3(0,0,0);

	if(!_node->within_covering(_pt))
		return;

	// if we hit a leaf, then add its contribution
	if(_node->is_leaf())  {
		double linear_f = _node->function_eval(_pt);
		double weight = _node->bspline2(_pt);
		Vector3 linear_n = _node->gradient_eval(_pt);
		Vector3 weight_gradient = _node->gradient_bspline2(_pt);

		numer1 = weight;
		numer2 = weight*linear_f;

		grad2 = weight_gradient;
		weight_gradient.scale(linear_f);
		linear_n.scale(weight);
		grad1 = weight_gradient+linear_n;

		return;
	}

	// otherwise, we are an interior node and within covering -> continue traversal
	for(int c = 0; c < 8; c++)  {
		double child_numer1 = 0, child_numer2 = 0;
		Vector3 child_grad1(0,0,0), child_grad2(0,0,0);
		this->octree_gradient(_pt, _node->children[c], child_numer1, child_numer2, child_grad1, child_grad2);

		numer1 += child_numer1;
		numer2 += child_numer2;
		grad1 = grad1 + child_grad1;
		grad2 = grad2 + child_grad2;
	}
}

DenseMatrix MPUSoup::hessian(Vector3 _pt)  {
	return DenseMatrix();
}

void MPUSoup::bounds(Vector3& ll_bound, Vector3& ur_bound)  {
	ll_bound = mpu->bounds.ll;
	ur_bound = mpu->bounds.ur;
}

void MPUSoup::init()  {
	Vector3* vertex_normals = new Vector3[polygon_soup->n_vertices()];
	for(int v = 0; v < polygon_soup->n_vertices(); v++)
		vertex_normals[v] = Vector3(0,0,0);

	num_soup = 0;
	for(int f = 0; f < polygon_soup->n_faces(); f++)  {
		BasicTriMesh::FaceHandle face(f);
		BasicTriMesh::FaceVertexIter fv_it = polygon_soup->fv_iter(face);
		OpenMesh::VertexHandle id0 = fv_it.handle();
		OpenMesh::VertexHandle id1 = (++fv_it).handle();
		OpenMesh::VertexHandle id2 = (++fv_it).handle();
		BasicTriMesh::Point p0 = polygon_soup->point(id0), p1 = polygon_soup->point(id1), p2 = polygon_soup->point(id2);
		Vector3 v0(p0[0],p0[1],p0[2]), v1(p1[0],p1[1],p1[2]), v2(p2[0],p2[1],p2[2]);

		Vector3 normal = (v1-v0).cross(v2-v0);
		double area = 0.5*normal.length();
		normal.normalize();
		normal.scale(area);

		if(area < 1e-11)
			continue;
		num_soup++;

		vertex_normals[id0.idx()] = vertex_normals[id0.idx()] + normal;
		vertex_normals[id1.idx()] = vertex_normals[id1.idx()] + normal;
		vertex_normals[id2.idx()] = vertex_normals[id2.idx()] + normal;
	}

	for(int v = 0; v < polygon_soup->n_vertices(); v++)
		vertex_normals[v].normalize();

	all_tris = new MPUTri*[num_soup];
	int tri_ind = 0;
	for(int f = 0; f < polygon_soup->n_faces(); f++)  {
		BasicTriMesh::FaceHandle face(f);
		BasicTriMesh::FaceVertexIter fv_it = polygon_soup->fv_iter(face);
		OpenMesh::VertexHandle id0 = fv_it.handle();
		OpenMesh::VertexHandle id1 = (++fv_it).handle();
		OpenMesh::VertexHandle id2 = (++fv_it).handle();
		BasicTriMesh::Point p0 = polygon_soup->point(id0), p1 = polygon_soup->point(id1), p2 = polygon_soup->point(id2);
		Vector3 v0(p0[0],p0[1],p0[2]), v1(p1[0],p1[1],p1[2]), v2(p2[0],p2[1],p2[2]);

		Vector3 normal = (v1-v0).cross(v2-v0);
		double area = 0.5*normal.length();
		if(area < 1e-11)
			continue;

		Vector3 n0 = vertex_normals[id0.idx()], n1 = vertex_normals[id1.idx()], n2 = vertex_normals[id2.idx()];
		all_tris[tri_ind++] = new MPUTri(v0, v1, v2, n0, n1, n2);
	}

	point_array = annAllocPts(num_soup, 3);
	for(int f = 0; f < num_soup; f++)  {
		MPUTri* tri = all_tris[f];
		Vector3 point = (tri->v1+tri->v2+tri->v3)*(1.0/3.0);
		point_array[f] = annAllocPt(3);
		point_array[f][0] = point.x;
		point_array[f][1] = point.y;
		point_array[f][2] = point.z;
	}
	kd_tree = new ANNkd_tree(point_array, num_soup, 3);

	mpu = new MPUOctree(all_tris, kd_tree, num_soup, min_samples, fit_eps, covering);
	mpu->build_tree();
}

// ------ MPUSoup ------ //

// ------ InvMPUSoup ------ //

InvMPUSoup::InvMPUSoup(BasicTriMesh* _polygonSoup, int _minSamples, double _fitEps, double _covering)  {
	polygon_soup = _polygonSoup;
	min_samples = _minSamples;
	fit_eps = _fitEps;
	covering = _covering;
	my_pi = acos(-1);
	init();
}

InvMPUSoup::InvMPUSoup(string _filename)  {
	FILE* file_in = fopen(_filename.c_str(), "rb");
	size_t num_read;

	double full_bounds[6];
	num_read = fread(full_bounds, sizeof(double), 6, file_in);
	AABB mpu_bounds(Vector3(full_bounds[0],full_bounds[1],full_bounds[2]), Vector3(full_bounds[3],full_bounds[4],full_bounds[5]));
	num_read = fread(&min_samples, sizeof(int), 1, file_in);
	num_read = fread(&fit_eps, sizeof(double), 1, file_in);
	num_read = fread(&covering, sizeof(double), 1, file_in);

	mpu = new MPUOctree(file_in, min_samples, fit_eps, covering, mpu_bounds);

	fclose(file_in);
}

InvMPUSoup::~InvMPUSoup()  {
}

void InvMPUSoup::write_to_file(string _filename)  {
	FILE* file_out = fopen(_filename.c_str(), "wb");
	mpu->write_to_file(file_out);
	fclose(file_out);
}

double InvMPUSoup::function(Vector3 _pt)  {
	num_interior = num_leaf = 0;
	debug = false;

	double numer, denom;
	this->octree_func(_pt, mpu, numer, denom);

	if(denom == 0)  {
		cout << "zero denom! " << _pt << " bounds : " << mpu->bounds.ll << " , " << mpu->bounds.ur << endl;
		debug = true;
		this->octree_func(_pt, mpu, numer, denom);
	}
	else  {
		//cout << "func: " << (numer/denom) << endl;
	}

	return -(numer / denom);
}

void InvMPUSoup::octree_func(Vector3 _pt, MPUOctree* _node, double& func_numer, double& func_denom)  {
	func_numer = func_denom = 0;
	if(!_node->within_covering(_pt))
		return;

	// if we hit a leaf, then add its contribution
	if(_node->is_leaf())  {
		double weight, function_evaluation;
		_node->mpu_eval(_pt, function_evaluation, weight);
		func_numer = weight*function_evaluation;
		func_denom = weight;
		if(debug)
			cout << "function : " << function_evaluation << " weight : " << weight << endl;
		return;
	}

	// otherwise, we are an interior node and within covering -> continue traversal
	for(int c = 0; c < 8; c++)  {
		double child_numer = 0, child_denom = 0;
		this->octree_func(_pt, _node->children[c], child_numer, child_denom);
		func_numer += child_numer;
		func_denom += child_denom;
	}
}

Vector3 InvMPUSoup::gradient(Vector3 _pt)  {
	double numer1 = 0, numer2 = 0;
	Vector3 grad1(0,0,0), grad2(0,0,0);
	this->octree_gradient(_pt, mpu, numer1, numer2, grad1, grad2);

	grad1.scale(numer1);
	grad2.scale(-numer2);
	double weight_scale = 1.0/(numer1*numer1);

	Vector3 final_gradient = grad1+grad2;
	final_gradient.scale(-1.0*weight_scale);

	return final_gradient;
}

Vector3 InvMPUSoup::normal(Vector3 _pt)  {
	Vector3 normal_vec = this->gradient(_pt);
	normal_vec.normalize();
	return normal_vec;
}

void InvMPUSoup::octree_gradient(Vector3 _pt, MPUOctree* _node, double& numer1, double& numer2,
																				 Vector3& grad1, Vector3& grad2)  {
	numer1 = numer2 = 0;
	grad1 = Vector3(0,0,0);
	grad2 = Vector3(0,0,0);

	if(!_node->within_covering(_pt))
		return;

	// if we hit a leaf, then add its contribution
	if(_node->is_leaf())  {
		double linear_f = _node->function_eval(_pt);
		double weight = _node->bspline2(_pt);
		Vector3 linear_n = _node->gradient_eval(_pt);
		Vector3 weight_gradient = _node->gradient_bspline2(_pt);

		numer1 = weight;
		numer2 = weight*linear_f;

		grad2 = weight_gradient;
		weight_gradient.scale(linear_f);
		linear_n.scale(weight);
		grad1 = weight_gradient+linear_n;

		return;
	}

	// otherwise, we are an interior node and within covering -> continue traversal
	for(int c = 0; c < 8; c++)  {
		double child_numer1 = 0, child_numer2 = 0;
		Vector3 child_grad1(0,0,0), child_grad2(0,0,0);
		this->octree_gradient(_pt, _node->children[c], child_numer1, child_numer2, child_grad1, child_grad2);

		numer1 += child_numer1;
		numer2 += child_numer2;
		grad1 = grad1 + child_grad1;
		grad2 = grad2 + child_grad2;
	}
}

DenseMatrix InvMPUSoup::hessian(Vector3 _pt)  {
	return DenseMatrix();
}

void InvMPUSoup::bounds(Vector3& ll_bound, Vector3& ur_bound)  {
	ll_bound = mpu->bounds.ll;
	ur_bound = mpu->bounds.ur;
}

void InvMPUSoup::init()  {
	Vector3* vertex_normals = new Vector3[polygon_soup->n_vertices()];
	for(int v = 0; v < polygon_soup->n_vertices(); v++)
		vertex_normals[v] = Vector3(0,0,0);

	num_soup = 0;
	for(int f = 0; f < polygon_soup->n_faces(); f++)  {
		BasicTriMesh::FaceHandle face(f);
		BasicTriMesh::FaceVertexIter fv_it = polygon_soup->fv_iter(face);
		OpenMesh::VertexHandle id0 = fv_it.handle();
		OpenMesh::VertexHandle id1 = (++fv_it).handle();
		OpenMesh::VertexHandle id2 = (++fv_it).handle();
		BasicTriMesh::Point p0 = polygon_soup->point(id0), p1 = polygon_soup->point(id1), p2 = polygon_soup->point(id2);
		Vector3 v0(p0[0],p0[1],p0[2]), v1(p1[0],p1[1],p1[2]), v2(p2[0],p2[1],p2[2]);

		Vector3 normal = (v1-v0).cross(v2-v0);
		double area = 0.5*normal.length();
		normal.normalize();
		normal.scale(area);

		if(area < 1e-11)
			continue;
		num_soup++;

		vertex_normals[id0.idx()] = vertex_normals[id0.idx()] + normal;
		vertex_normals[id1.idx()] = vertex_normals[id1.idx()] + normal;
		vertex_normals[id2.idx()] = vertex_normals[id2.idx()] + normal;
	}

	for(int v = 0; v < polygon_soup->n_vertices(); v++)
		vertex_normals[v].normalize();

	all_tris = new MPUTri*[num_soup];
	int tri_ind = 0;
	for(int f = 0; f < polygon_soup->n_faces(); f++)  {
		BasicTriMesh::FaceHandle face(f);
		BasicTriMesh::FaceVertexIter fv_it = polygon_soup->fv_iter(face);
		OpenMesh::VertexHandle id0 = fv_it.handle();
		OpenMesh::VertexHandle id1 = (++fv_it).handle();
		OpenMesh::VertexHandle id2 = (++fv_it).handle();
		BasicTriMesh::Point p0 = polygon_soup->point(id0), p1 = polygon_soup->point(id1), p2 = polygon_soup->point(id2);
		Vector3 v0(p0[0],p0[1],p0[2]), v1(p1[0],p1[1],p1[2]), v2(p2[0],p2[1],p2[2]);

		Vector3 normal = (v1-v0).cross(v2-v0);
		double area = 0.5*normal.length();
		if(area < 1e-11)
			continue;

		Vector3 n0 = vertex_normals[id0.idx()], n1 = vertex_normals[id1.idx()], n2 = vertex_normals[id2.idx()];
		all_tris[tri_ind++] = new MPUTri(v0, v1, v2, n0, n1, n2);
	}

	point_array = annAllocPts(num_soup, 3);
	for(int f = 0; f < num_soup; f++)  {
		MPUTri* tri = all_tris[f];
		Vector3 point = (tri->v1+tri->v2+tri->v3)*(1.0/3.0);
		point_array[f] = annAllocPt(3);
		point_array[f][0] = point.x;
		point_array[f][1] = point.y;
		point_array[f][2] = point.z;
	}
	kd_tree = new ANNkd_tree(point_array, num_soup, 3);

	mpu = new MPUOctree(all_tris, kd_tree, num_soup, min_samples, fit_eps, covering);
	mpu->build_tree();
}

// ------ InvMPUSoup ------ //

