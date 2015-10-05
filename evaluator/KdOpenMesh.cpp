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

#include "KdOpenMesh.h"

KdOpenMesh::KdOpenMesh(BasicTriMesh* _mesh)  {
	mesh = _mesh;
	this->construct_tree();
}

KdOpenMesh::~KdOpenMesh()  {
	delete kd_tree;
}

bool KdOpenMesh::intersect(Vector3 _origin, Vector3 _dir, Vector3& tri_pt, Vector3& tri_normal)  {
	Intersection intersection;
	Ray ray( Point(_origin.x,_origin.y,_origin.z), Vector(_dir.x,_dir.y,_dir.z), 1e-10 );
	if(kd_tree->Intersect(ray, &intersection))  {
		DifferentialGeometry diff_geom = intersection.dg;
		Point geom_p = diff_geom.p;
		tri_pt = Vector3(geom_p.x,geom_p.y,geom_p.z);
		Normal geom_normal = diff_geom.nn;
		tri_normal = Vector3(geom_normal.x,geom_normal.y,geom_normal.z);
		return true;
	}

	// else no intersection
	return false;
}

bool KdOpenMesh::intersect(Vector3 _origin, Vector3 _dir, Vector3& tri_pt, int& tri_index)  {
	return false;
}

void KdOpenMesh::construct_tree()  {
	Transform identity;
	bool reverse = false;
	int n_tris = mesh->n_faces();
	int n_verts = mesh->n_vertices();
	int* tri_list = new int[3*n_tris];
	Point* point_list = new Point[n_verts];

	for(int p = 0; p < n_verts; p++)  {
		BasicTriMesh::Point next_pt = mesh->point(OpenMesh::VertexHandle(p));
		Point new_pt = Point(next_pt[0], next_pt[1], next_pt[2]);
		point_list[p].x = next_pt[0];
		point_list[p].y = next_pt[1];
		point_list[p].z = next_pt[2];
	}
	for(int t = 0; t < n_tris; t++)  {
		BasicTriMesh::FaceVertexIter fv_it = mesh->fv_iter(OpenMesh::FaceHandle(t));
		BasicTriMesh::VertexHandle v0 = fv_it.handle();
		BasicTriMesh::VertexHandle v1 = (++fv_it).handle();
		BasicTriMesh::VertexHandle v2 = (++fv_it).handle();

		tri_list[3*t] = v0.idx();
		tri_list[3*t+1] = v1.idx();
		tri_list[3*t+2] = v2.idx();
	}

	// shape
	ParamSet shape_params;
	shape_params.AddPoint("P", point_list, n_verts);
	shape_params.AddInt("nv", &n_verts, 1);
	shape_params.AddInt("nt", &n_tris, 1);
	shape_params.AddInt("indices", tri_list, 3*n_tris);
	shape = CreateTriangleMeshShape(&identity, &identity, reverse, shape_params);

	// dummy material
	ParamSet geomp, matp;
	map<string, Reference<Texture<float> > > ft;
	map<string, Reference<Texture<Spectrum> > > st;
	TextureParams texture_params(geomp, matp, ft, st);
	dummy_material = CreateMatteMaterial(identity, texture_params);
	if(dummy_material == NULL)
		std::cout << "bad material?" << std::endl;

	// primitive
	primitive = new GeometricPrimitive(shape, dummy_material, NULL);

	// and finally, the kd-tree
	ParamSet kd_params;
	int max_prims = 2, max_depth = 22;
	kd_params.AddInt("maxprims", &max_prims, 1);
	kd_params.AddInt("maxdepth", &max_depth, 1);
	vector<Reference<Primitive> > single_primitive;
	Reference<Primitive> prim_ref = primitive;
	single_primitive.push_back(prim_ref);
	std::cout << "creating kd tree..." << std::flush;
	kd_tree = CreateKdTreeAccelerator(single_primitive, kd_params);
	std::cout << "... done creating kd tree..." << std::endl;
}

void KdOpenMesh::ray_tracer(string _filename, int _view, double _samplingRadius, int _res, string _algorithm)  {
	// algorithm colors
	Vector3 lgreen_color = this->hash_to_rgb("B6DFA7");
	Vector3 green_color = this->hash_to_rgb("7FBF67");
	Vector3 lpink_color = this->hash_to_rgb("E9AEB7");
	Vector3 pink_color = this->hash_to_rgb("D37281");
	Vector3 lblue_color = this->hash_to_rgb("96B8C6");
	Vector3 blue_color = this->hash_to_rgb("4E7B8C");
	Vector3 lyellow_color = this->hash_to_rgb("EFE3B3");
	Vector3 yellow_color = this->hash_to_rgb("E0CA79");
	Vector3 lpurple_color = this->hash_to_rgb("A6A0CD");
	Vector3 purple_color = this->hash_to_rgb("655D98");
	Vector3 gold_color = Vector3((236.0/255.0), (211.0/255.0), (143.0/255.0));

	// determine algorithm color
	Vector3 main_color = gold_color;
	if(_algorithm.compare("apss") == 0)
		main_color = green_color;
	else if(_algorithm.compare("fourier") == 0)
		main_color = pink_color;
	else if(_algorithm.compare("imls") == 0)
		main_color = blue_color;
	else if(_algorithm.compare("mpu") == 0)
		main_color = yellow_color;
	else if(_algorithm.compare("mpusmooth") == 0)
		main_color = purple_color;
	else if(_algorithm.compare("poisson") == 0)
		main_color = lgreen_color;
	else if(_algorithm.compare("rbf") == 0)
		main_color = lpink_color;
	else if(_algorithm.compare("scattered") == 0)
		main_color = lblue_color;
	else if(_algorithm.compare("spss") == 0)
		main_color = lyellow_color;
	else if(_algorithm.compare("wavelet") == 0)
		main_color = lpurple_color;

	// figure out sensor settings
	double my_pi = acos(-1);
	string particle_file = "./particle_sampler/particles30.npts";
	OrientedPointCloud sphere_pc;
	NPTSReader npts_reader(particle_file);
	if(!npts_reader.read_pc(&sphere_pc))  {
		cerr << "unable to read pc " << particle_file << "!" << endl;
		return;
	}
	Vector3 center(0,0,0);

	Vector3 sphere_sample = sphere_pc.getPoint(_view);
	Vector3 look_from = center + sphere_sample*_samplingRadius;
	Vector3 look_at = center;
	Vector3 view_dir = look_at-look_from;
	view_dir.normalize();

	Vector3 up;
	if(view_dir.x == 0 && view_dir.y == 1 && view_dir.z == 0)
		up = Vector3(0,0,1);
	else  {
		Vector3 fake_up(0,1,0);
		Vector3 crossed = view_dir.cross(fake_up);
		crossed.normalize();
		up = crossed.cross(view_dir);
		up.normalize();
	}

	bool reverse_up = false;
	if(reverse_up)
		up.scale(-1.0);

	RayCaster sensor(look_from, view_dir, up, _res, _res, 50, my_pi/4.0);


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

			Vector3 hit_pt, hit_normal;
			sensor.generatePerspectiveRay(x, y, ray_pos, ray_dir);
			if(this->intersect(ray_pos, ray_dir, hit_pt, hit_normal))  {
				num_hits++;
				num_valid_hits++;

				Vector3 light_vec = sensor.camera_pos()-hit_pt;
				light_vec.normalize();
				double dot = fabs(light_vec.dotProduct(hit_normal));
				r_scan[y][x] = main_color.x*dot;
				g_scan[y][x] = main_color.y*dot;
				b_scan[y][x] = main_color.z*dot;
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
}

Vector3 KdOpenMesh::hash_to_rgb(string _hex)  {
	int ind = 0;

	Vector3 rgb_val;
	for(int i = 0; i < 3; i++)  {
		double next_num = 0;

		for(int a = 0; a < 2; a++)  {
			double base_val = 0;
			string hex_digit = _hex.substr(ind,1);

			if(hex_digit.compare("A") == 0)
				base_val = 10;
			else if(hex_digit.compare("B") == 0)
				base_val = 11;
			else if(hex_digit.compare("C") == 0)
				base_val = 12;
			else if(hex_digit.compare("D") == 0)
				base_val = 13;
			else if(hex_digit.compare("E") == 0)
				base_val = 14;
			else if(hex_digit.compare("F") == 0)
				base_val = 15;
			else
				base_val = atof(hex_digit.c_str());

			next_num = (a == 1) ? next_num+base_val : (next_num+16*base_val);
			ind++;
		}

		if(i == 0)
			rgb_val.x = next_num / 255.0;
		else if(i == 1)
			rgb_val.y = next_num / 255.0;
		else if(i == 2)
			rgb_val.z = next_num / 255.0;
	}

	return rgb_val;
}
