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

#include "tri_integrator.h"

TriIntegrator::TriIntegrator(Vector3 _point, Vector3 _v1, Vector3 _v2, Vector3 _v3, double _eps, double _lambda)  {
	v1 = _v1;
	v2 = _v2;
	v3 = _v3;
	point = _point;
	eps = _eps;

	normal = (v2-v1).cross(v3-v1);
	area = 0.5*normal.length();
	normal.normalize();

	verify = false;

	// setup integration, determine if approximation is ok
	this->tri_integrals();
	double approx_lambda = _lambda;
	double shortest_length = point.length(closest_point);
	double d1 = v1.length(point), d2 = v2.length(point), d3 = v3.length(point);
	double inv_shortest_length = 1.0 / shortest_length;

	use_approximation = (inv_shortest_length*d1) < approx_lambda && (inv_shortest_length*d2) < approx_lambda &&
							  (inv_shortest_length*d3) < approx_lambda;
}

TriIntegrator::TriIntegrator(Vector3 _point, Vector3 _v1, Vector3 _v2, Vector3 _v3, double _eps)  {
	v1 = _v1;
	v2 = _v2;
	v3 = _v3;
	point = _point;
	eps = _eps;

	normal = (v2-v1).cross(v3-v1);
	area = 0.5*normal.length();
	normal.normalize();

	verify = false;
	use_approximation = false;

	this->tri_integrals();
}

// assumption: ONLY LINE INTEGRALS!
TriIntegrator::TriIntegrator(Vector3 _pt, double _eps)  {
	point = _pt;
	eps = _eps;
}

TriIntegrator::~TriIntegrator()  {
}

double TriIntegrator::weight_eval(Vector3 _p1, Vector3 _p2)  {
	return 1.0 / (_p1.sqdDist(_p2) + (eps*eps));
}

double TriIntegrator::weight_eval_sqd(Vector3 _p1, Vector3 _p2)  {
	double main_weight = _p1.sqdDist(_p2) + (eps*eps);
	return 1.0 / (main_weight*main_weight);
}

double TriIntegrator::weight_eval_cbd(Vector3 _p1, Vector3 _p2)  {
	double main_weight = _p1.sqdDist(_p2) + (eps*eps);
	return 1.0 / (main_weight*main_weight*main_weight);
}


double TriIntegrator::constant_integration()  {
	if(use_approximation)  {
		Vector3 bary = (v1+v2+v3)*(1.0/3.0);
		double w = this->weight_eval_sqd(bary, point);
		return area*w;
	}
	// else
	return this->numerical_constant_integration();
}

double TriIntegrator::numerical_constant_integration()  {
	double full_integral = 0, full_area = 0;
	double* inner_integration = new double[num_steps];

	for(int t = 0; t < num_subtris; t++)  {
		SubTri sub_tri = subdivided_tris[t];
		Vector3 weighty_pt, lesser_pt;
		if(this->weight_eval(sub_tri.v2, point) > this->weight_eval(sub_tri.v3, point))  {
			weighty_pt = sub_tri.v2;
			lesser_pt = sub_tri.v3;
		}
		else  {
			weighty_pt = sub_tri.v3;
			lesser_pt = sub_tri.v2;
		}

		Vector3 edge_perp = (lesser_pt-weighty_pt).cross(normal);
		edge_perp.normalize();
		Vector3 edge_dir = weighty_pt-sub_tri.o;
		edge_dir.normalize();
		Vector3 integral_dir = edge_dir.cross(normal);

		// inner analytical integration: solve at each step point
		for(unsigned i = 0; i < num_steps; i++)  {
			double alpha = nonuniform_steps[i];
			Vector3 e1 = sub_tri.o*(1.0-alpha) + lesser_pt*alpha;
			Vector3 e2 = ray_plane_intersection(weighty_pt, edge_perp, e1, edge_dir);
			if(e1 == e2)
				cout << "positions equal?? " << alpha << endl;
			double analytical_integral = this->constant_line_integration(e1, e2);
			//double analytical_integral = this->weight_eval_sqd((e1+e2)*0.5, point);
			inner_integration[i] = analytical_integral;
		}

		// outer numerical integration: trapezoidal rule
		double sub_integral = 0;
		for(unsigned i = 1; i < num_steps; i++)  {
			double alpha_0 = nonuniform_steps[i-1];
			Vector3 e1_0 = sub_tri.o*(1.0-alpha_0) + lesser_pt*alpha_0;
			double alpha_1 = nonuniform_steps[i];
			Vector3 e1_1 = sub_tri.o*(1.0-alpha_1) + lesser_pt*alpha_1;

			Vector3 proj_pt = orthogonal_projection(e1_0, e1_1, integral_dir);
			double step_size = e1_0.length(proj_pt);

			sub_integral += 0.5*(inner_integration[i-1]+inner_integration[i])*step_size;
		}
		full_integral += sub_integral;
	}
	delete [] inner_integration;

	// VERIFICATION: monte-carlo integration
	if(verify)  {
		int num_samples = 900000;
		double mc_integral = 0;
		for(int i = 0; i < num_samples; i++)  {
			double u = (double)rand() / (double)RAND_MAX;
			double v = (double)rand() / (double)RAND_MAX;
			double tmp = sqrt(u);
			double b1 = 1.0-tmp;
			double b2 = v*tmp;
			double b3 = (1.0-b1-b2);

			Vector3 rand_pt = v1*b1 + v2*b2 + v3*b3;
			mc_integral += this->weight_eval_sqd(rand_pt, point);
		}
		mc_integral *= (area / (double)num_samples);
		cout << "constant quadrature integral: " << full_integral << " ; mc integral: " << mc_integral << endl;
	}

	return full_integral;
}

double TriIntegrator::linear_integration(TriFunction _f)  {
	if(use_approximation)  {
		Vector3 bary = (v1+v2+v3)*(1.0/3.0);
		double w = this->weight_eval_sqd(bary, point);
		double bary_f = (_f.f1+_f.f2+_f.f3)/3.0;
		return bary_f*area*w;
	}
	// else
	return this->numerical_linear_integration(_f);
}

double TriIntegrator::numerical_linear_integration(TriFunction _f)  {
	double full_integral = 0, full_area = 0;
	double* inner_integration = new double[num_steps];

	for(int t = 0; t < num_subtris; t++)  {
		SubTri sub_tri = subdivided_tris[t];
		Vector3 weighty_pt, lesser_pt;
		if(this->weight_eval(sub_tri.v2, point) > this->weight_eval(sub_tri.v3, point))  {
			weighty_pt = sub_tri.v2;
			lesser_pt = sub_tri.v3;
		}
		else  {
			weighty_pt = sub_tri.v3;
			lesser_pt = sub_tri.v2;
		}

		Vector3 edge_perp = (lesser_pt-weighty_pt).cross(normal);
		edge_perp.normalize();
		Vector3 edge_dir = weighty_pt-sub_tri.o;
		edge_dir.normalize();
		Vector3 integral_dir = edge_dir.cross(normal);

		double b0_1, b1_1, b2_1;
		double b0_2, b1_2, b2_2;

		// inner analytical integration: solve at each step point
		for(unsigned i = 0; i < num_steps; i++)  {
			double alpha = nonuniform_steps[i];
			Vector3 e1 = sub_tri.o*(1.0-alpha) + lesser_pt*alpha;
			this->barycentric_coordinate(e1, b0_1, b1_1, b2_1);
			Vector3 e2 = ray_plane_intersection(weighty_pt, edge_perp, e1, edge_dir);
			this->barycentric_coordinate(e2, b0_2, b1_2, b2_2);

			double f1 = _f.f1*b0_1 + _f.f2*b1_1 + _f.f3*b2_1;
			double f2 = _f.f1*b0_2 + _f.f2*b1_2 + _f.f3*b2_2;
			double analytical_integral = this->linear_line_integration(e1, e2, f1, f2);
			inner_integration[i] = analytical_integral;
		}

		// outer numerical integration: trapezoidal rule
		double sub_integral = 0;
		for(unsigned i = 1; i < num_steps; i++)  {
			double alpha_0 = nonuniform_steps[i-1];
			Vector3 e1_0 = sub_tri.o*(1.0-alpha_0) + lesser_pt*alpha_0;
			double alpha_1 = nonuniform_steps[i];
			Vector3 e1_1 = sub_tri.o*(1.0-alpha_1) + lesser_pt*alpha_1;

			Vector3 proj_pt = orthogonal_projection(e1_0, e1_1, integral_dir);
			double step_size = e1_0.length(proj_pt);

			sub_integral += 0.5*(inner_integration[i-1]+inner_integration[i])*step_size;
		}
		full_integral += sub_integral;
	}
	delete [] inner_integration;

	// VERIFICATION: monte-carlo integration
	if(verify)  {
		int num_samples = 900000;
		double mc_integral = 0;
		for(int i = 0; i < num_samples; i++)  {
			double u = (double)rand() / (double)RAND_MAX;
			double v = (double)rand() / (double)RAND_MAX;
			double tmp = sqrt(u);
			double b1 = 1.0-tmp;
			double b2 = v*tmp;
			double b3 = (1.0-b1-b2);

			Vector3 rand_pt = v1*b1 + v2*b2 + v3*b3;
			double rand_f = _f.f1*b1 + _f.f2*b2 + _f.f3*b3;
			mc_integral += rand_f*this->weight_eval_sqd(rand_pt, point);
		}
		mc_integral *= (area / (double)num_samples);
		cout << "linear quadrature integral["<<num_subtris<<"] : " << full_integral << " ; mc integral: " << mc_integral << endl;
	}

	return full_integral;
}

double TriIntegrator::product_integration(TriFunction _f, TriFunction _g)  {
	if(use_approximation)  {
		Vector3 bary = (v1+v2+v3)*(1.0/3.0);
		double w = this->weight_eval_sqd(bary, point);
		double bary_f = (_f.f1+_f.f2+_f.f3)/3.0;
		double bary_g = (_g.f1+_g.f2+_g.f3)/3.0;
		return bary_f*bary_g*area*w;
	}
	// else
	return this->numerical_product_integration(_f, _g);
}

double TriIntegrator::numerical_product_integration(TriFunction _f, TriFunction _g)  {
	double full_integral = 0, full_area = 0;
	for(int t = 0; t < num_subtris; t++)  {
		SubTri sub_tri = subdivided_tris[t];
		Vector3 weighty_pt, lesser_pt;
		if(this->weight_eval(sub_tri.v2, point) > this->weight_eval(sub_tri.v3, point))  {
			weighty_pt = sub_tri.v2;
			lesser_pt = sub_tri.v3;
		}
		else  {
			weighty_pt = sub_tri.v3;
			lesser_pt = sub_tri.v2;
		}

		Vector3 edge_perp = (lesser_pt-weighty_pt).cross(normal);
		edge_perp.normalize();
		Vector3 edge_dir = weighty_pt-sub_tri.o;
		edge_dir.normalize();
		Vector3 integral_dir = edge_dir.cross(normal);

		double b0_1, b1_1, b2_1;
		double b0_2, b1_2, b2_2;

		// inner analytical integration: solve at each step point
		vector<double> inner_integration;
		for(unsigned i = 0; i < num_steps; i++)  {
			double alpha = nonuniform_steps[i];
			Vector3 e1 = sub_tri.o*(1.0-alpha) + lesser_pt*alpha;
			this->barycentric_coordinate(e1, b0_1, b1_1, b2_1);
			Vector3 e2 = ray_plane_intersection(weighty_pt, edge_perp, e1, edge_dir);
			this->barycentric_coordinate(e2, b0_2, b1_2, b2_2);

			double f1 = _f.f1*b0_1 + _f.f2*b1_1 + _f.f3*b2_1;
			double g1 = _g.f1*b0_1 + _g.f2*b1_1 + _g.f3*b2_1;
			double f2 = _f.f1*b0_2 + _f.f2*b1_2 + _f.f3*b2_2;
			double g2 = _g.f1*b0_2 + _g.f2*b1_2 + _g.f3*b2_2;
			double analytical_integral = this->quadratic_line_integration(e1, e2, f1, f2, g1, g2);
			inner_integration.push_back(analytical_integral);
		}

		// outer numerical integration: trapezoidal rule
		double sub_integral = 0;
		for(unsigned i = 1; i < num_steps; i++)  {
			double alpha_0 = nonuniform_steps[i-1];
			Vector3 e1_0 = sub_tri.o*(1.0-alpha_0) + lesser_pt*alpha_0;
			double alpha_1 = nonuniform_steps[i];
			Vector3 e1_1 = sub_tri.o*(1.0-alpha_1) + lesser_pt*alpha_1;

			Vector3 proj_pt = orthogonal_projection(e1_0, e1_1, integral_dir);
			double step_size = e1_0.length(proj_pt);

			sub_integral += 0.5*(inner_integration[i-1]+inner_integration[i])*step_size;
		}
		full_integral += sub_integral;
	}

	// VERIFICATION: monte-carlo integration
	if(verify)  {
		int num_samples = 900000;
		double mc_integral = 0;
		for(int i = 0; i < num_samples; i++)  {
			double u = (double)rand() / (double)RAND_MAX;
			double v = (double)rand() / (double)RAND_MAX;
			double tmp = sqrt(u);
			double b1 = 1.0-tmp;
			double b2 = v*tmp;
			double b3 = (1.0-b1-b2);

			Vector3 rand_pt = v1*b1 + v2*b2 + v3*b3;
			double rand_f = _f.f1*b1 + _f.f2*b2 + _f.f3*b3;
			double rand_g = _g.f1*b1 + _g.f2*b2 + _g.f3*b3;
			mc_integral += rand_f*rand_g*this->weight_eval_sqd(rand_pt, point);
		}
		mc_integral *= (area / (double)num_samples);
		cout << "quadratic quadrature integral: " << full_integral << " ; mc integral: " << mc_integral << endl;
	}

	return full_integral;
}

Vector3 TriIntegrator::gradient_integration()  {
	if(use_approximation)  {
		Vector3 bary = (v1+v2+v3)*(1.0/3.0);
		double w = this->weight_eval_cbd(bary, point);
		return (bary-point)*4.0*area*w;
	}
	// else
	return this->numerical_gradient_integration();
}

Vector3 TriIntegrator::numerical_gradient_integration()  {
	Vector3 full_integral(0,0,0);
	Vector3* inner_integration = new Vector3[num_steps];

	for(int t = 0; t < num_subtris; t++)  {
		SubTri sub_tri = subdivided_tris[t];
		Vector3 weighty_pt, lesser_pt;
		if(this->weight_eval(sub_tri.v2, point) > this->weight_eval(sub_tri.v3, point))  {
			weighty_pt = sub_tri.v2;
			lesser_pt = sub_tri.v3;
		}
		else  {
			weighty_pt = sub_tri.v3;
			lesser_pt = sub_tri.v2;
		}

		Vector3 edge_perp = (lesser_pt-weighty_pt).cross(normal);
		edge_perp.normalize();
		Vector3 edge_dir = weighty_pt-sub_tri.o;
		edge_dir.normalize();
		Vector3 integral_dir = edge_dir.cross(normal);

		// inner analytical integration: solve at each step point
		for(unsigned i = 0; i < num_steps; i++)  {
			double alpha = nonuniform_steps[i];
			Vector3 e1 = sub_tri.o*(1.0-alpha) + lesser_pt*alpha;
			Vector3 e2 = ray_plane_intersection(weighty_pt, edge_perp, e1, edge_dir);
			inner_integration[i] = this->gradient_line_integration(e1, e2);
		}

		// outer numerical integration: trapezoidal rule
		Vector3 sub_integral(0,0,0);
		for(unsigned i = 1; i < num_steps; i++)  {
			double alpha_0 = nonuniform_steps[i-1];
			Vector3 e1_0 = sub_tri.o*(1.0-alpha_0) + lesser_pt*alpha_0;
			double alpha_1 = nonuniform_steps[i];
			Vector3 e1_1 = sub_tri.o*(1.0-alpha_1) + lesser_pt*alpha_1;

			Vector3 proj_pt = orthogonal_projection(e1_0, e1_1, integral_dir);
			double step_size = e1_0.length(proj_pt);

			sub_integral = sub_integral + (inner_integration[i-1]+inner_integration[i])*0.5*step_size;
		}
		full_integral = full_integral + sub_integral;
	}
	delete [] inner_integration;

	// VERIFICATION: monte-carlo integration
	if(verify)  {
		int num_samples = 900000;
		Vector3 mc_integral(0,0,0);
		for(int i = 0; i < num_samples; i++)  {
			double u = (double)rand() / (double)RAND_MAX;
			double v = (double)rand() / (double)RAND_MAX;
			double tmp = sqrt(u);
			double b1 = 1.0-tmp;
			double b2 = v*tmp;
			double b3 = (1.0-b1-b2);

			Vector3 rand_pt = v1*b1 + v2*b2 + v3*b3;
			mc_integral = mc_integral + (rand_pt-point)*4.0*this->weight_eval_cbd(rand_pt,point);
		}
		mc_integral = mc_integral*(area / (double)num_samples);
		cout << "gradient quadrature integral: " << full_integral << " ; mc integral: " << mc_integral << endl;
	}

	return full_integral;
}


double TriIntegrator::inv_quadratic_integral(double _a, double _b, double _c)  {
	double c1 = sqrt(4.0*_a*_c-_b*_b);
	return (2.0/c1)*(atan( (2.0*_a+_b)/(c1) ) - atan( (_b)/(c1) ));
	//return 2.0/c1;
}

double TriIntegrator::inv_sqd_quadratic_integral(double _a, double _b, double _c)  {
	double c1_sqd = 4.0*_a*_c-_b*_b;
	double factor1 = ( ((2.0*_a+_b) / (_a+_b+_c)) - (_b/_c) ) / c1_sqd;
	double factor2 = (2.0*_a) / c1_sqd;
	return factor1 + factor2*this->inv_quadratic_integral(_a,_b,_c);
}

double TriIntegrator::inv_cubic_integral(double _a, double _b, double _c)  {
	double descr = 4.0*_a*_c-_b*_b;
	double sqd_sum = (_a+_b+_c)*(_a+_b+_c);
	double factor1 = ( ((2.0*_a+_b)/sqd_sum) - (_b/(_c*_c)) ) / (2.0*descr);
	double factor2 = (3.0*_a) / descr;
	return factor1 + factor2*this->inv_sqd_quadratic_integral(_a,_b,_c);
}

double TriIntegrator::inv_linear_sqd_quadratic_integral(double _a, double _b, double _c)  {
	double c1_sqd = 4.0*_a*_c-_b*_b;
	double factor1 = -( ((_b+2.0*_c) / (_a+_b+_c)) - 2.0 ) / c1_sqd;
	double factor2 = -_b / c1_sqd;
	return factor1 + factor2*this->inv_quadratic_integral(_a,_b,_c);
}

double TriIntegrator::inv_linear_cubic_integral(double _a, double _b, double _c)  {
	double descr = 4.0*_a*_c-_b*_b;
	double sqd_sum = (_a+_b+_c)*(_a+_b+_c);
	double factor1 = -( ((_b+2.0*_c)/sqd_sum) - (2.0/_c) ) / (2.0*descr);
	double factor2 = -(3.0*_b)/(2.0*descr);
	return factor1 + factor2*this->inv_sqd_quadratic_integral(_a,_b,_c);
}

double TriIntegrator::constant_line_integration(Vector3 _p1, Vector3 _p2)  {
	Vector3 v1 = point-_p1, v2 = _p1-_p2;
	double a = v2.dotProduct(v2), b = 2.0*v1.dotProduct(v2), c = v1.dotProduct(v1) + eps*eps;
	double descr = 4.0*a*c-b*b;
	if(descr == 0)  {
		cout << "line " << _p1 << " " << _p2 << " : " << point << endl;
	}
	double line_length = _p1.length(_p2);
	double full_integral = line_length*this->inv_sqd_quadratic_integral(a,b,c);
	return full_integral;
}

double TriIntegrator::linear_line_integration(Vector3 _p1, Vector3 _p2, double _f1, double _f2)  {
	Vector3 v1 = point-_p1, v2 = _p1-_p2;
	double a = v2.dotProduct(v2), b = 2.0*v1.dotProduct(v2), c = v1.dotProduct(v1) + eps*eps;
	double line_length = _p1.length(_p2);
	double full_integral = line_length*(_f1*this->inv_sqd_quadratic_integral(a,b,c) + (_f2-_f1)*this->inv_linear_sqd_quadratic_integral(a,b,c));
	return full_integral;
}

double TriIntegrator::quadratic_line_integration(Vector3 _p1, Vector3 _p2, double _f1, double _f2, double _g1, double _g2)  {
	Vector3 v1 = point-_p1, v2 = _p1-_p2;
	double a = v2.dotProduct(v2), b = 2.0*v1.dotProduct(v2), c = v1.dotProduct(v1) + eps*eps;
	double line_length = _p1.length(_p2);
	double alpha_1 = (_f2-_f1)*(_g2-_g1);
	double alpha_2 = (_f2-_f1)*_g1 + (_g2-_g1)*_f1;
	double alpha_3 = _f1*_g1;

	double simple_factor = (alpha_1/a)*this->inv_quadratic_integral(a,b,c);
	double const_factor = (alpha_3 - (alpha_1*c)/a)*this->inv_sqd_quadratic_integral(a,b,c);
	double linear_factor = (alpha_2 - (alpha_1*b)/a)*this->inv_linear_sqd_quadratic_integral(a,b,c);

	double full_integral = line_length*(simple_factor+const_factor+linear_factor);
	return full_integral;
}

Vector3 TriIntegrator::gradient_line_integration(Vector3 _p1, Vector3 _p2)  {
	Vector3 v1 = point-_p1, v2 = _p1-_p2;
	double a = v2.dotProduct(v2), b = 2.0*v1.dotProduct(v2), c = v1.dotProduct(v1) + eps*eps;
	double line_length = _p1.length(_p2);

	Vector3 factor1 = point;
	Vector3 factor2_1 = _p1, factor2_2 = _p2-_p1;

	double integral1 = line_length*this->inv_cubic_integral(a,b,c);
	double integral2 = line_length*this->inv_linear_cubic_integral(a,b,c);

	return ((factor2_1*integral1 + factor2_2*integral2) - factor1*integral1)*4.0;
}

void TriIntegrator::tri_integrals()  {
	if(point == v1)  {
		num_subtris = 1;
		subdivided_tris[0] = SubTri(v1, v2, v3);
		closest_point = v1;
		return;
	}
	else if(point == v2)  {
		num_subtris = 1;
		subdivided_tris[0] = SubTri(v2, v3, v1);
		closest_point = v2;
		return;
	}
	else if(point == v3)  {
		num_subtris = 1;
		subdivided_tris[0] = SubTri(v3, v1, v2);
		closest_point = v3;
		return;
	}

	// --- first step: project point onto plane, and subdivide triangle based on projected point --- //
	Vector3 orthog_pt = point + normal*(v1-point).dotProduct(normal);
	Vector3 area_normal1 = (v2-orthog_pt).cross(v3-orthog_pt);
	Vector3 area_normal2 = (v3-orthog_pt).cross(v1-orthog_pt);
	Vector3 area_normal3 = (v1-orthog_pt).cross(v2-orthog_pt);
	bool is_area1_positive = area_normal1.dotProduct(normal) > 0;
	bool is_area2_positive = area_normal2.dotProduct(normal) > 0;
	bool is_area3_positive = area_normal3.dotProduct(normal) > 0;

	// is our projection in the triangle?
	if(is_area1_positive && is_area2_positive && is_area3_positive)  {
		num_subtris = 3;
		subdivided_tris[0] = SubTri(orthog_pt, v1, v2);
		subdivided_tris[1] = SubTri(orthog_pt, v2, v3);
		subdivided_tris[2] = SubTri(orthog_pt, v3, v1);
		closest_point = orthog_pt;
		return;
	}

	// otherwise, find the point on the border of the triangle closest to orthog_pt
	Vector3 tri_verts[3];
	tri_verts[0] = v1; tri_verts[1] = v2; tri_verts[2] = v3;

	// first try the edges ...
	bool is_contained = false;
	for(int i = 0; i < 3; i++)  {
		Vector3 q1 = tri_verts[i], q2 = tri_verts[(i+1)%3], q3 = tri_verts[(i+2)%3];
		Vector3 e1 = q2-q1, e2 = q1-q2;
		Vector3 edge_perp = e2.cross(normal);

		is_contained = signed_dist(orthog_pt, q1, edge_perp) > 0 && signed_dist(orthog_pt, q1, e1) > 0 && signed_dist(orthog_pt, q2, e2) > 0;
		if(is_contained)  {
			Vector3 edge_point = this->edge_projection(q1, q2, normal, orthog_pt);
			num_subtris = 2;
			subdivided_tris[0] = SubTri(edge_point, q2, q3);
			subdivided_tris[1] = SubTri(edge_point, q3, q1);
			closest_point = edge_point;
			return;
		}
	}

	// ... and last, find closest vertex
	int min_ind = 2;
	double sqd_dist1 = v1.sqdDist(orthog_pt), sqd_dist2 = v2.sqdDist(orthog_pt), sqd_dist3 = v3.sqdDist(orthog_pt);
	if(sqd_dist1 < sqd_dist2 && sqd_dist1 < sqd_dist3)
		min_ind = 0;
	else if(sqd_dist2 < sqd_dist3)
		min_ind = 1;

	Vector3 q1 = tri_verts[min_ind], q2 = tri_verts[(min_ind+1)%3], q3 = tri_verts[(min_ind+2)%3];
	num_subtris = 1;
	subdivided_tris[0] = SubTri(q1,q2,q3);
	closest_point = q1;
}

Vector3 TriIntegrator::ray_plane_intersection(Vector3 _p, Vector3 _n, Vector3 _m, Vector3 _d)  {
	double t = (_p-_m).dotProduct(_n) / _d.dotProduct(_n);
	return _m + _d*t;
}
Vector3 TriIntegrator::orthogonal_projection(Vector3 _p, Vector3 _x, Vector3 _normal)  {
	return _p + _normal*((_p-_x).dotProduct(_normal)/_normal.dotProduct(_normal));
}
double TriIntegrator::signed_dist(Vector3 _p, Vector3 _x, Vector3 _normal)  {
	return ((_p-_x).dotProduct(_normal)) / _normal.dotProduct(_normal);
}
Vector3 TriIntegrator::edge_projection(Vector3 _v1, Vector3 _v2, Vector3 _normal, Vector3 _point)  {
	// get normal of edge -> just cross of edge with normal
	Vector3 edge_normal = (_v2-_v1).cross(_normal);
	edge_normal.normalize();
	// projection of _point onto plane spanned by edge normal
	return _point + edge_normal*(_v1-_point).dotProduct(edge_normal);
}

void TriIntegrator::barycentric_coordinate(Vector3 _point, double& b1, double& b2, double& b3)  {
	double area1 = 0.5*((_point-v2).cross(_point-v3)).length();
	double area2 = 0.5*((_point-v3).cross(_point-v1)).length();
	double area3 = 0.5*((_point-v1).cross(_point-v2)).length();
	double inv_area = 1.0 / area;
	b1 = area1*inv_area;
	b2 = area2*inv_area;
	b3 = area3*inv_area;
}
