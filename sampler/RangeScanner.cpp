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

#include "RangeScanner.h"

// ------ RangeSubset ------ //

RangeSubset::RangeSubset(int _resx, int _resy, double _gaussianSigma)  {
	res_x = _resx;
	res_y = _resy;
	gauss_sigma = _gaussianSigma;

	double falloff = 0.005;
	filter_half = ceil( sqrt(-0.5*log(falloff)*gauss_sigma*gauss_sigma) );
	filter_half = filter_half < 2 ? 2 : filter_half;
	filter_size = 2*filter_half+1;

	gaussian_filter = new double[filter_size];
	double sig_sqd2 = gauss_sigma*gauss_sigma;
	for(int x = -filter_half; x <= filter_half; x++)
		gaussian_filter[x+filter_half] = 1.0/sqrt(2.0*acos(-1)*sig_sqd2) * exp(-(2.0*x*x)/(sig_sqd2));

	start_x = 1e8;
	end_x = -1e8;
}

RangeSubset::~RangeSubset()  {
	delete [] gaussian_filter;

	if(contains_nothing)
		return;

	for(int y = 0; y < res_y; y++)
		delete [] intensities[y];
	delete [] intensities;
}

void RangeSubset::expand_range(int _x)  {
	start_x = _x < start_x ? _x : start_x;
	end_x = _x > end_x ? _x : end_x;
}

void RangeSubset::finalize_range()  {
	if(start_x == 1e8)  {
		contains_nothing = true;
		return;
	}
	contains_nothing = false;

	start_x -= (filter_half+2);
	if(start_x < 0)
		start_x = 0;
	end_x += (filter_half+2);
	if(end_x >= res_x)
		end_x = res_x-1;
	image_width = (end_x-start_x)+1;

	intensities = new double*[res_y];
	for(int y = 0; y < res_y; y++)  {
		intensities[y] = new double[image_width];
		for(int x = 0; x < image_width; x++)
			intensities[y][x] = 0;
	}
}

void RangeSubset::contribute_point(int _x, int _y, double _intensity)  {
	int true_x = _x-start_x;
	intensities[_y][true_x] = _intensity;
}

void RangeSubset::smooth()  {
	if(contains_nothing)
		return;

	double** smoothed_signal = new double*[res_y];
	for(int y = 0; y < res_y; y++)  {
		smoothed_signal[y] = new double[image_width];
		for(int x = 0; x < image_width; x++)  {
			double smoothed_val = 0, total_weight = 0;

			for(int nx = -filter_half; nx <= filter_half; nx++)  {
				int next_x = nx+x;
				if(next_x < 0)
					next_x = 0;
				if(next_x >= image_width)
					next_x = image_width-1;

				for(int ny = -filter_half; ny <= filter_half; ny++)  {
					int next_y = ny+y;
					if(next_y < 0)
						next_y = 0;
					if(next_y >= res_y)
						next_y = res_y-1;

					double weight = gaussian_filter[nx+filter_half]*gaussian_filter[ny+filter_half];
					smoothed_val += weight*intensities[next_y][next_x];
					total_weight += weight;
				}

			}

			smoothed_signal[y][x] = (total_weight == 0) ? 0 : smoothed_val / total_weight;
		}
	}

	for(int y = 0; y < res_y; y++)  {
		for(int x = 0; x < image_width; x++)
			intensities[y][x] = smoothed_signal[y][x];
	}

	for(int y = 0; y < res_y; y++)
		delete [] smoothed_signal[y];
	delete [] smoothed_signal;
}

void RangeSubset::write_to_file(string _filename, double** _r, double** _g, double** _b)  {
	if(contains_nothing)
		return;

	PNGImage png_image(_filename, res_x, res_y);
	for(int y = 0; y < res_y; y++)  {
		for(int x = 0; x < res_x; x++)  {
			double r = _r[y][x], g = _g[y][x], b = _b[y][x];
			if(x >= start_x && x <= end_x)  {
				double intensity = intensities[y][x-start_x];
				r = 1.0*intensity + r*(1.0-intensity);
				g = 0.0*intensity + g*(1.0-intensity);
				b = 0.0*intensity + b*(1.0-intensity);
			}
			png_image.set_png_pixel(x,(res_y-y)-1, r, g, b);
		}
	}
	png_image.write_to_file();
}

// ------ RangeOctree ------ //

RangeOctree::RangeOctree(SparseRangeScan* _scan)  {
	scan = _scan;
	ll = Vector3(1e10,1e10,1e10);
	ur = Vector3(-1e10,-1e10,-1e10);
	depth = 0;

	int num_points = scan->size();
	for(int p = 0; p < num_points; p++)  {
		point_inds.push_back(p);
		Vector3 scan_pt = scan->getPt(p);
		ll.expand_min_bound(scan_pt);
		ur.expand_max_bound(scan_pt);
	}
	this->finalize();
}

RangeOctree::RangeOctree(SparseRangeScan* _scan, int _depth)  {
	scan = _scan;
	depth = _depth+1;
	ll = Vector3(1e10,1e10,1e10);
	ur = Vector3(-1e10,-1e10,-1e10);
}

RangeOctree::~RangeOctree()  {
	if(!isLeaf)  {
		for(int i = 0; i < 8; i++)
			delete children[i];
		delete [] children;
	}
}

void RangeOctree::construct_tree()  {
	isLeaf = this->size() < 2;
	if(this->is_leaf())
		return;

	children = new RangeOctree*[8];
	int child_ind = 0;
	for(int x = 0; x < 2; x++)  {
		for(int y = 0; y < 2; y++)  {
			for(int z = 0; z < 2; z++)
				children[child_ind++] = new RangeOctree(scan, depth+1);
		}
	}

	double x_bounds[] = { ll.x , 0.5*(ll.x+ur.x) , ur.x };
	double y_bounds[] = { ll.y , 0.5*(ll.y+ur.y) , ur.y };
	double z_bounds[] = { ll.z , 0.5*(ll.z+ur.z) , ur.z };

	if(fabs(x_bounds[0]-x_bounds[2]) < 1e-5)  {
		x_bounds[0] = ll.x-0.001;
		x_bounds[2] = ur.x+0.001;
		x_bounds[1] = 0.5*(x_bounds[0]+x_bounds[2]);
	}
	if(fabs(y_bounds[0]-y_bounds[2]) < 1e-5)  {
		y_bounds[0] = ll.y-0.001;
		y_bounds[2] = ur.y+0.001;
		y_bounds[1] = 0.5*(y_bounds[0]+y_bounds[2]);
	}
	if(fabs(z_bounds[0]-z_bounds[2]) < 1e-5)  {
		z_bounds[0] = ll.z-0.001;
		z_bounds[2] = ur.z+0.001;
		z_bounds[1] = 0.5*(z_bounds[0]+z_bounds[2]);
	}

	for(int p = 0; p < this->size(); p++)  {
		int next_pt = point_inds[p];
		Vector3 scan_pt = scan->getPt( next_pt );
		bool found_cell = false;

		int child_ind = 0;
		for(int x = 0; x < 2; x++)  {
			for(int y = 0; y < 2; y++)  {
				for(int z = 0; z < 2; z++)  {
					Vector3 child_ll(x_bounds[x], y_bounds[y], z_bounds[z]);
					Vector3 child_ur(x_bounds[x+1], y_bounds[y+1], z_bounds[z+1]);
					if( (scan_pt.x >= child_ll.x && scan_pt.x < child_ur.x) &&
						 (scan_pt.y >= child_ll.y && scan_pt.y < child_ur.y) &&
						 (scan_pt.z >= child_ll.z && scan_pt.z < child_ur.z) )  {
						found_cell = true;
						children[child_ind]->add_point(next_pt);
					}	
					child_ind++;

					if(found_cell) break;
				}
				if(found_cell) break;
			}
			if(found_cell) break;
		}

		if(!found_cell)  {
			cerr << "["<<depth<<"] pt " << scan_pt << " doesn't belong to any child cell???" << endl;
			cerr << ll << " : " << ur << endl;
		}
	}

	point_inds.clear();

	for(int i = 0; i < 8; i++)
		children[i]->finalize();
	for(int i = 0; i < 8; i++)
		children[i]->construct_tree();
}

void RangeOctree::finalize()  {
	Vector3 diag = ur-ll;
	ll = ll - (diag*0.001);
	ur = ur + (diag*0.001);
}

void RangeOctree::add_point(int _p)  {
	point_inds.push_back(_p);
	Vector3 scan_pt = scan->getPt(_p);
	ll.expand_min_bound(scan_pt);
	ur.expand_max_bound(scan_pt);
}

void RangeOctree::intersect_stripe(LaserStripe _stripe, vector<int>* overlapping_inds)  {
	// no intersection against intermediate node -> return
	if(!isLeaf && !this->stripe_aabb_intersection(_stripe))  {
		//cout << "failed intersection at depth " << depth << endl;
		return;
	}

	// otherwise, if we are a parent, then keep going ...
	if(!isLeaf)  {
		for(int i = 0; i < 8; i++)
			children[i]->intersect_stripe(_stripe, overlapping_inds);
		return;
	}

	// if we are a leaf, then add point(s) to point inds, if they are inside
	for(unsigned i = 0; i < point_inds.size(); i++)  {
		Vector3 pt = scan->getPt( point_inds[i] );
		if(_stripe.minus_signed_distance(pt) > 0 && _stripe.plus_signed_distance(pt) > 0)
			overlapping_inds->push_back(point_inds[i]);
	}
}

bool RangeOctree::stripe_aabb_intersection(LaserStripe _stripe)  {
	// first check to see if the aabb center is inside the stripe
	Vector3 center = (ll+ur)*0.5;
	if(_stripe.minus_signed_distance(center) > 0 && _stripe.plus_signed_distance(center) > 0)
		return true;

	// now check to see if our aabb intersects the stripe - test against plus and minus planes
	if(this->aabb_plane_intersection(_stripe.n_minus, _stripe.d_minus) || this->aabb_plane_intersection(_stripe.n_plus, _stripe.d_plus))
		return true;

	// otherwise, outside
	return false;
}

bool RangeOctree::aabb_plane_intersection(Vector3 _n, double _d)  {
	Vector3 p_vertex, n_vertex;

	if(_n.x > 0)  {
		p_vertex.x = ur.x;
		n_vertex.x = ll.x;
	}
	else  {
		n_vertex.x = ur.x;
		p_vertex.x = ll.x;
	}

	if(_n.y > 0)  {
		p_vertex.y = ur.y;
		n_vertex.y = ll.y;
	}
	else  {
		n_vertex.y = ur.y;
		p_vertex.y = ll.y;
	}

	if(_n.z > 0)  {
		p_vertex.z = ur.z;
		n_vertex.z = ll.z;
	}
	else  {
		n_vertex.z = ur.z;
		p_vertex.z = ll.z;
	}

	if((p_vertex.dotProduct(_n) + _d) < 0)
		return false;
	else if((n_vertex.dotProduct(_n) + _d) > 0)
		return false;
	return true;
}

// ------ RangeScanner ------ //

RangeScanner::RangeScanner(ImplicitFunction* _shape, RayCaster _sensor, Vector3 _laserPos)  {
	shape = _shape;
	sensor = _sensor;
	laser_pos = _laserPos;
	lipschitz = 0;
	cout << "laser pos: " << laser_pos << endl;

	min_range = 127;
	max_range = 230;
	num_stripes = sensor.resx();

	double default_thickness = 10, default_fov = 45;
	double fov_thickness = (default_fov / (double)num_stripes) * default_thickness;
	laser_fov = (fov_thickness/180.0)*acos(-1);

	additive_noise = 0.0;
	laser_smoother = 0.2;

	peak_threshold = 0.05;
	std_threshold = 1.5;

	dump_stripes = false;
	stripes_processed = 0;
}

RangeScanner::RangeScanner(ImplicitFunction* _shape, RayCaster _sensor, Vector3 _laserPos, double _lipschitz)  {
	shape = _shape;
	sensor = _sensor;
	laser_pos = _laserPos;
	lipschitz = _lipschitz;
	cout << "laser pos: " << laser_pos << endl;

	min_range = 127;
	max_range = 230;
	num_stripes = sensor.resx();

	double default_thickness = 10, default_fov = 45;
	double fov_thickness = (default_fov / (double)num_stripes) * default_thickness;
	laser_fov = (fov_thickness/180.0)*acos(-1);

	//additive_noise = 0.35;
	//laser_smoother = 5.0;
	additive_noise = 0.0;
	laser_smoother = 0.2;

	peak_threshold = 0.05;
	std_threshold = 1.5;

	dump_stripes = false;
	stripes_processed = 0;
}

RangeScanner::~RangeScanner()  {
}

SparseRangeScan* RangeScanner::optical_triangulation()  {
	// first construct image...
	this->image_sensor();

	RangeOctree range_octree(&sparse_scan);
	range_octree.construct_tree();

	// first determine starting and ending points
	int mid_y = sensor.resy()/2, start_x = sensor.resx(), end_x = 0;
	int mid_x = sensor.resx()/2, start_y = sensor.resy(), end_y = 0;
	Vector3 start_pos, dir_start, end_pos, dir_end;
	sensor.generatePerspectiveRay(start_x, mid_y, start_pos, dir_start);
	sensor.generatePerspectiveRay(end_x, mid_y, end_pos, dir_end);
	//sensor.generatePerspectiveRay(mid_x, start_y, start_pos, dir_start);
	//sensor.generatePerspectiveRay(mid_x, end_y, end_pos, dir_end);
	Vector3 start_pt = start_pos + dir_start*max_range;
	Vector3 end_pt = end_pos + dir_end*min_range;

	// the start and end points determine a frustum for the laser - determine frustum angle
	Vector3 laser_dir_start = start_pt-laser_pos;
	laser_dir_start.normalize();
	Vector3 laser_dir_end = end_pt-laser_pos;
	laser_dir_end.normalize();
	double laser_frustum = acos(laser_dir_start.dotProduct(laser_dir_end));
	cout << "laser frustum: " << ((laser_frustum/acos(-1))*180.0) << endl;

	// from the laser frustum and number of stripes, determine angle increment
	double angle_increment = laser_frustum / (double)num_stripes;
	double laser_shift = angle_increment/2.0;

	stripes_processed = 0;

	// populate true image
	double** r_scan = new double*[sensor.resy()];
	double** g_scan = new double*[sensor.resy()];
	double** b_scan = new double*[sensor.resy()];
	for(int y = 0; y < sensor.resy(); y++)  {
		r_scan[y] = new double[sensor.resx()];
		g_scan[y] = new double[sensor.resx()];
		b_scan[y] = new double[sensor.resx()];
		for(int x = 0; x < sensor.resx(); x++)  {
			r_scan[y][x] = 1;
			g_scan[y][x] = 1;
			b_scan[y][x] = 1;
		}
	}

	// compute analytical normals, while also populating original range scan image
	Vector3* pt_normals = new Vector3[sparse_scan.size()];
	for(int p = 0; p < sparse_scan.size(); p++)  {
		Vector3 scan_pt = sparse_scan.getPt(p);
		Vector3 shape_normal = shape->normal(scan_pt);
		pt_normals[p] = shape_normal;
		PixelCoord pix = sparse_scan.getPixel(p);
		Vector3 light_vec = sensor.camera_pos()-scan_pt;
		light_vec.normalize();

		double dot = fabs(light_vec.dotProduct(shape_normal));
		r_scan[pix.y][pix.x] = (236.0/255.0)*dot;
		g_scan[pix.y][pix.x] = (211.0/255.0)*dot;
		b_scan[pix.y][pix.x] = (143.0/255.0)*dot;
	}

	// ... add points which sensor saw, but laser didn't
	for(int p = 0; p < image_scan.size(); p++)  {
		Vector3 scan_pt = image_scan.getPt(p);
		Vector3 shape_normal = shape->normal(scan_pt);
		PixelCoord pix = image_scan.getPixel(p);
		Vector3 light_vec = sensor.camera_pos()-scan_pt;
		light_vec.normalize();

		double dot = fabs(light_vec.dotProduct(shape_normal));
		r_scan[pix.y][pix.x] = (236.0/255.0)*dot;
		g_scan[pix.y][pix.x] = (211.0/255.0)*dot;
		b_scan[pix.y][pix.x] = (143.0/255.0)*dot;
	}

	// initialize gaussian fits
	GaussianFit*** gaussian_fits = new GaussianFit**[sensor.resy()];
	for(int y = 0; y < sensor.resy(); y++)  {
		gaussian_fits[y] = new GaussianFit*[sensor.resx()];
		for(int x = 0; x < sensor.resx(); x++)
			gaussian_fits[y][x] = new GaussianFit();
	}

	int fraction = num_stripes/5;

	// now go through each laser shift, and perform optical triangulation
	for(int s = 0; s < num_stripes; s++)  {
		if((s % fraction) == 0)
			cout << ((double)s / (double)num_stripes) << " done with stripe... " << endl;

		laser_shift += angle_increment;
		LaserStripe stripe;
		this->construct_laser_stripe(start_pt, laser_shift, laser_fov, stripe);

		// gather all points which the laser stripe intersects
		vector<int> intersected_geometry;
		range_octree.intersect_stripe(stripe, &intersected_geometry);

		// construct image subset, first determine area coverage
		RangeSubset range_subset(sensor.resx(), sensor.resy(), laser_smoother);
		for(unsigned g = 0; g < intersected_geometry.size(); g++)  {
			int pt_ind = intersected_geometry[g];
			PixelCoord pix = sparse_scan.getPixel(pt_ind);
			range_subset.expand_range(pix.x);
		}
		range_subset.finalize_range();

		// populate image with laser
		for(unsigned g = 0; g < intersected_geometry.size(); g++)  {
			int pt_ind = intersected_geometry[g];
			PixelCoord pix = sparse_scan.getPixel(pt_ind);
			Vector3 pt = sparse_scan.getPt(pt_ind);

			double light_sigma = pt.length(laser_pos)*tan(0.5*laser_fov);
			double pulse_dist = stripe.pulse_abs_distance(pt);
			double luminance = this->gaussian_beam(pulse_dist, 0.5*light_sigma);
			luminance = luminance < 0 ? 0 : luminance;

			double noise_scale = (1.0-luminance);
			double noise_magnitude = luminance < additive_noise ? luminance : additive_noise;
			double gaussian_noise = this->sampleNormalDistribution(noise_scale, noise_magnitude);

			Vector3 laser_ray = laser_pos-pt;
			laser_ray.normalize();
			Vector3 pt_normal = pt_normals[pt_ind];
			double cos_falloff = laser_ray.dotProduct(pt_normal);

			double radiance = luminance*cos_falloff+gaussian_noise;
			radiance = radiance > 1 ? 1 : radiance;
			radiance = radiance < 0 ? 0 : radiance;
			int quantized_radiance = 255*radiance;
			radiance = quantized_radiance / 255.0;
			range_subset.contribute_point(pix.x, pix.y, radiance);
		}

		if(range_subset.containsNothing())
			continue;

		// smooth image
		range_subset.smooth();

		if(dump_stripes)  {
			char* file_char = new char[300];
			int s_f = stripe_start+stripes_processed+1;
			if(s_f < 10)
				sprintf(file_char, "%s0000%u.png", base_stripe.c_str(), s_f);
			else if(s_f < 100)
				sprintf(file_char, "%s000%u.png", base_stripe.c_str(), s_f);
			else if(s_f < 1000)
				sprintf(file_char, "%s00%u.png", base_stripe.c_str(), s_f);
			else if(s_f < 10000)
				sprintf(file_char, "%s0%u.png", base_stripe.c_str(), s_f);
			else
				sprintf(file_char, "%s%u.png", base_stripe.c_str(), s_f);
			string filename = file_char;
			range_subset.write_to_file(filename, r_scan, g_scan, b_scan);
			delete [] file_char;
		}
		stripes_processed++;

		// gather center-of-gravity times
		for(int y = 0; y < sensor.resy(); y++)  {
			for(int x = range_subset.startX(); x <= range_subset.endX(); x++)  {
				double intensity = range_subset.get_intensity(x,y);
				if(intensity > 1e-6)
					gaussian_fits[y][x]->add_point(laser_shift, intensity);
			}
		}
	}

	for(int y = 0; y < sensor.resy(); y++)  {
		delete [] r_scan[y];
		delete [] g_scan[y];
		delete [] b_scan[y];
	}
	delete [] r_scan;
	delete [] g_scan;
	delete [] b_scan;

	double outlier_threshold = 4.0;
	SparseRangeScan* depth_scan = new SparseRangeScan();
	int num_retained = 0, num_peak_rejections = 0, num_std_rejections = 0;
	double ave_weight = 0;
	for(int p = 0; p < sparse_scan.size(); p++)  {
		Vector3 scan_pt = sparse_scan.getPt(p);
		PixelCoord pix = sparse_scan.getPixel(p);

		if(gaussian_fits[pix.y][pix.x]->do_fit())  {
			double peak = gaussian_fits[pix.y][pix.x]->get_peak();
			double std = gaussian_fits[pix.y][pix.x]->get_std();
			if(peak < peak_threshold)  {
				//cout << "rejected " << pix.x << " " << pix.y << " due to peak threshold..." << endl;
				num_peak_rejections++;
				continue;
			}
			if(std > (0.5*std_threshold*laser_fov))  {
				//cout << "rejected " << pix.x << " " << pix.y << " due to std threshold..." << endl;
				num_std_rejections++;
				continue;
			}

			double proper_time = gaussian_fits[pix.y][pix.x]->get_mean();
			num_retained++;
			Vector3 ray_pos, ray_dir;
			sensor.generatePerspectiveRay(pix.x, pix.y, ray_pos, ray_dir);
			double reference_depth = ray_pos.length(scan_pt);

			LaserStripe stripe;
			this->construct_laser_stripe(start_pt, proper_time, laser_fov, stripe);
			double local_depth = -(ray_pos.dotProduct(stripe.n)+stripe.d) / ray_dir.dotProduct(stripe.n);

			if(fabs(local_depth-reference_depth) > outlier_threshold)
				continue;

			Vector3 depth_pt = ray_pos + ray_dir*local_depth;
			depth_scan->addPoint(pix, depth_pt);
		}
	}
	cout << "retained " << num_retained << " / " << sparse_scan.size() << " peak rejections: " <<
		num_peak_rejections << " ; std rejections: " << num_std_rejections << endl;

	for(int y = 0; y < sensor.resy(); y++)  {
		for(int x = 0; x < sensor.resx(); x++)
			delete gaussian_fits[y][x];
	}
	for(int y = 0; y < sensor.resy(); y++)
		delete [] gaussian_fits[y];
	delete [] gaussian_fits;
	delete [] pt_normals;

	return depth_scan;
}

void RangeScanner::interpolate_stripe(LaserStripe _stripe, Vector3 _pt, Vector3& n, double& d, double& w)  {
	double plus_dist = _stripe.plus_signed_distance(_pt);
	double minus_dist = _stripe.minus_signed_distance(_pt);
	double total_dist = plus_dist+minus_dist;

	if(plus_dist < 0 || minus_dist < 0)
		cout << "BADNESS IN INTERSECTION TESTING!!!" << endl;

	double alpha_minus = minus_dist / total_dist;
	Vector3 n_plus = _stripe.n_plus, n_minus = _stripe.n_minus*-1.0;
	double d_plus = _stripe.d_plus, d_minus = -_stripe.d_minus;
	n = n_minus*(1.0-alpha_minus) + n_plus*alpha_minus;
	d = d_minus*(1.0-alpha_minus) + d_plus*alpha_minus;

	w = 1.0 - 2.0*fabs(alpha_minus-0.5);
}

void RangeScanner::interpolate_stripe(Vector3 _startPt, double _shift, double _laserFOV, Vector3 _pt, Vector3& n, double& d, double& w)  {
	Vector3 laser_start = _startPt-laser_pos;
	laser_start.normalize();
	Vector3 baseline = sensor.camera_pos() - laser_pos;
	baseline.normalize();

	Vector3 laser_pt = _pt-laser_pos;
	laser_pt.normalize();
	double pt_angle = acos(laser_pt.dotProduct(laser_start));

	//cout << "min shift: " << (_shift-0.5*_laserFOV) << " : " << pt_angle << " max shift: " << (_shift+0.5*_laserFOV) << endl;

	Vector3 rotation_axis = laser_start.cross(baseline);
	rotation_axis.normalize();
	RotationMatrix rotation(rotation_axis, pt_angle);
	Vector3 rotated_laser_centered = rotation.rotate(laser_start);

	n = rotation_axis.cross(rotated_laser_centered);
	d = -(n.dotProduct(laser_pos));
	w = 1;
}

void RangeScanner::construct_laser_stripe(Vector3 _startPt, double _shift, double _laserFOV, LaserStripe& stripe)  {
	Vector3 laser_start = _startPt-laser_pos;
	laser_start.normalize();
	Vector3 baseline = sensor.camera_pos() - laser_pos;
	baseline.normalize();

	Vector3 rotation_axis = laser_start.cross(baseline);
	rotation_axis.normalize();

	// rotate laser start about rotation axis, by _shift amount and width _laserFOV
	RotationMatrix rotation(rotation_axis, _shift);
	Vector3 rotated_laser_centered = rotation.rotate(laser_start);
	RotationMatrix minus_rotation(rotation_axis, _shift-0.5*_laserFOV);
	Vector3 rotated_laser_minus = minus_rotation.rotate(laser_start);
	RotationMatrix plus_rotation(rotation_axis, _shift+0.5*_laserFOV);
	Vector3 rotated_laser_plus = plus_rotation.rotate(laser_start);

	//cout << "laser stripe: " << rotated_laser_minus << " : " << rotated_laser_centered << " : " << rotated_laser_plus << endl;

	// obtain laser stripe planes
	//Vector3 minus_n = rotated_laser_minus.cross(rotation_axis);
	Vector3 minus_n = rotation_axis.cross(rotated_laser_minus);
	minus_n.normalize();
	double minus_d = -(minus_n.dotProduct(laser_pos));

	//Vector3 plus_n = rotation_axis.cross(rotated_laser_plus);
	Vector3 plus_n = rotated_laser_plus.cross(rotation_axis);
	plus_n.normalize();
	double plus_d = -(plus_n.dotProduct(laser_pos));

	Vector3 centered_n = rotated_laser_centered.cross(rotation_axis);
	centered_n.normalize();
	double centered_d = -(centered_n.dotProduct(laser_pos));

	stripe = LaserStripe(minus_n, minus_d, centered_n, centered_d, plus_n, plus_d);
}

void RangeScanner::image_sensor()  {
	if(false)  {
		sparse_scan = SparseRangeScan("dense_scan.scan");
		return;
	}

	ImplicitNewton* implicit_solver = (lipschitz == 0) ? new ImplicitNewton(shape, 1.0, 1e-7) : new ImplicitNewton(shape, 1.0, 1e-7, lipschitz);
	double my_pi = acos(-1);

	int total_res = sensor.resx()*sensor.resy();
	int fraction = total_res/10;
	int ind = 0;
	int num_hits = 0, num_valid_hits = 0;

	Vector3 ray_pos, ray_dir;
	double t_hit;

	for(int y = 0; y < sensor.resy(); y++)  {
		for(int x = 0; x < sensor.resx(); x++)  {
			ind++;

			if((ind % fraction) == 0)
				cout << ((double)ind / (double)total_res) << " done with range image... " << num_valid_hits << " / " << num_hits << endl;

			sensor.generatePerspectiveRay(x, y, ray_pos, ray_dir);
			if(implicit_solver->intersection(ray_pos, ray_dir, t_hit))  {
				num_hits++;
				Vector3 range_pt = ray_pos + ray_dir*t_hit;

				Vector3 laser_dir = range_pt-laser_pos;
				laser_dir.normalize();

				double laser_t;
				if(!implicit_solver->intersection(laser_pos, laser_dir, laser_t))  {
					image_scan.addPoint(PixelCoord(x,y), range_pt);
					continue;
				}
				Vector3 laser_strike = laser_pos + laser_dir*laser_t;
				double laser_diff = laser_strike.length(range_pt);
				if(laser_diff > 1e-2)  {
					image_scan.addPoint(PixelCoord(x,y), range_pt);
					continue;
				}

				num_valid_hits++;
				sparse_scan.addPoint(PixelCoord(x,y), range_pt);
			}
		}
	}

	delete implicit_solver;
	//sparse_scan.write_to_file("dense_scan.scan");
}

void RangeScanner::ray_tracer(string _filename)  {
	ImplicitNewton* implicit_solver = (lipschitz == 0) ? new ImplicitNewton(shape, 1.0, 1e-7) : new ImplicitNewton(shape, 1.0, 1e-7, lipschitz);
	double my_pi = acos(-1);

	int total_res = sensor.resx()*sensor.resy();
	int fraction = total_res/10;
	int ind = 0;
	int num_hits = 0, num_valid_hits = 0;

	Vector3 ray_pos, ray_dir;
	double t_hit;

	double** r_scan = new double*[sensor.resy()];
	double** g_scan = new double*[sensor.resy()];
	double** b_scan = new double*[sensor.resy()];
	for(int y = 0; y < sensor.resy(); y++)  {
		r_scan[y] = new double[sensor.resx()];
		g_scan[y] = new double[sensor.resx()];
		b_scan[y] = new double[sensor.resx()];
		for(int x = 0; x < sensor.resx(); x++)  {
			r_scan[y][x] = 1;
			g_scan[y][x] = 1;
			b_scan[y][x] = 1;
		}
	}

	for(int y = 0; y < sensor.resy(); y++)  {
		for(int x = 0; x < sensor.resx(); x++)  {
			ind++;

			if((ind % fraction) == 0)
				cout << ((double)ind / (double)total_res) << " done with range image... " << num_valid_hits << " / " << num_hits << endl;

			sensor.generatePerspectiveRay(x, y, ray_pos, ray_dir);
			if(implicit_solver->intersection(ray_pos, ray_dir, t_hit))  {
				num_hits++;
				Vector3 range_pt = ray_pos + ray_dir*t_hit;
				Vector3 range_normal = shape->normal(range_pt);

				num_valid_hits++;

				Vector3 light_vec = sensor.camera_pos()-range_pt;
				light_vec.normalize();
				double dot = fabs(light_vec.dotProduct(range_normal));
				r_scan[y][x] = (236.0/255.0)*dot;
				g_scan[y][x] = (211.0/255.0)*dot;
				b_scan[y][x] = (143.0/255.0)*dot;
			}
		}
	}

	PNGImage png_image(_filename, sensor.resx(), sensor.resy());
	for(int y = 0; y < sensor.resy(); y++)  {
		for(int x = 0; x < sensor.resx(); x++)  {
			double r = r_scan[y][x], g = g_scan[y][x], b = b_scan[y][x];
			png_image.set_png_pixel(x,(sensor.resy()-y)-1, r, g, b);
		}
	}
	png_image.write_to_file();

	for(int y = 0; y < sensor.resy(); y++)  {
		delete [] r_scan[y];
		delete [] g_scan[y];
		delete [] b_scan[y];
	}
	delete [] r_scan;
	delete [] g_scan;
	delete [] b_scan;

	delete implicit_solver;
}

void RangeScanner::depth_image(SparseRangeScan* _scan, string _filename)  {
	double** r_scan = new double*[sensor.resy()];
	double** g_scan = new double*[sensor.resy()];
	double** b_scan = new double*[sensor.resy()];
	for(int y = 0; y < sensor.resy(); y++)  {
		r_scan[y] = new double[sensor.resx()];
		g_scan[y] = new double[sensor.resx()];
		b_scan[y] = new double[sensor.resx()];
		for(int x = 0; x < sensor.resx(); x++)  {
			r_scan[y][x] = 0;
			g_scan[y][x] = 0;
			b_scan[y][x] = 0;
		}
	}

	double min_depth = 1e10, max_depth = -1e10;
	for(int p = 0; p < _scan->size(); p++)  {
		PixelCoord pix = _scan->getPixel(p);
		Vector3 pt = _scan->getPt(p);
		double depth = pt.length(sensor.camera_pos());
		min_depth = depth < min_depth ? depth : min_depth;
		max_depth = depth > max_depth ? depth : max_depth;
	}
	for(int p = 0; p < _scan->size(); p++)  {
		PixelCoord pix = _scan->getPixel(p);
		Vector3 pt = _scan->getPt(p);
		double depth = pt.length(sensor.camera_pos());
		double normalized_depth = 1.0 - ((depth-min_depth) / (max_depth-min_depth));
		r_scan[pix.y][pix.x] = normalized_depth;
		g_scan[pix.y][pix.x] = normalized_depth;
		b_scan[pix.y][pix.x] = normalized_depth;
	}

	PNGImage png_image(_filename, sensor.resx(), sensor.resy());
	for(int y = 0; y < sensor.resy(); y++)  {
		for(int x = 0; x < sensor.resx(); x++)  {
			double r = r_scan[y][x], g = g_scan[y][x], b = b_scan[y][x];
			png_image.set_png_pixel(x,(sensor.resy()-y)-1, r, g, b);
		}
	}
	png_image.write_to_file();

	for(int y = 0; y < sensor.resy(); y++)  {
		delete [] r_scan[y];
		delete [] g_scan[y];
		delete [] b_scan[y];
	}
	delete [] r_scan;
	delete [] g_scan;
	delete [] b_scan;
}

double RangeScanner::sampleNormalDistribution(double _sigma, double _magnitude)  {
	/*
	double rand_val = _sigma*((double)rand() / (double)RAND_MAX);
	double gaussian_val = _magnitude*exp(-(2.0*rand_val*rand_val)/(_sigma*_sigma));
	double neg_val = (double)rand() / (double)RAND_MAX;
	return neg_val < 0.5 ? -gaussian_val : gaussian_val;
	*/
	double rand_val = 0;
	for(int i = 0; i < 12; i++)
		rand_val += ((double)rand() / (double)RAND_MAX);
	rand_val -= 6;

	double normal_x = rand_val*_sigma;
	double gaussian_val = _magnitude*exp(-(2.0*normal_x*normal_x)/(_sigma*_sigma));
	double neg_val = (double)rand() / (double)RAND_MAX;
	//return neg_val < 0.5 ? -gaussian_val : gaussian_val;
	return -gaussian_val;
}

double RangeScanner::gaussian_beam(double _x, double _sigma)  {
	return exp((-2.0*(_x*_x)) / (_sigma*_sigma));
}
