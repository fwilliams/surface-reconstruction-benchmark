/*
Copyright 2007 Szymon Rusinkiewicz
Princeton University

autorient.cc
Automatically orient meshes for VRIP
Fragile - pretty much works only for scanalyze output

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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "XForm.h"
#include "bsphere.h"
#include <iostream>
using namespace std;


// VRIP-style translation + quaternion
typedef Vec<7,double> VRIPxform;


// Convert VRIP bmesh line transform to an xform
xform VRIP2xform(const VRIPxform &v)
{
	// The VRIP convention for quaternions is that the real part is last
	// and THE @#$#@$!! THING IS $#@$@ INVERTED!!!   (Grumble, grumble)
	double s2 = sqrt(sqr(v[3]) + sqr(v[4]) + sqr(v[5]));
	double c2 = -v[6];  // Yes, this is what's needed...
	double angle = 2.0 * atan2(s2, c2);
	return xform::trans(v[0], v[1], v[2]) *
	       xform::rot(angle, v[3], v[4], v[5]);
}


// Convert xform to VRIP bmesh line
VRIPxform xform2VRIP(const xform &xf)
{
	VRIPxform v;
	v[0] = xf[12];
	v[1] = xf[13];
	v[2] = xf[14];
	v[3] = xf[6] - xf[9];
	v[4] = xf[8] - xf[2];
	v[5] = xf[1] - xf[4];
	v[6] = xf[0] + xf[5] + xf[10] + xf[15];
	double qscale = 1.0f / sqrt(sqr(v[3])+sqr(v[4])+sqr(v[5])+sqr(v[6]));
	v[6] *= qscale;
	qscale = -qscale; // Flipped quat convention 
	v[3] *= qscale;
	v[4] *= qscale;
	v[5] *= qscale;
	return v;
}


// Given a mesh filename, together with translation and rotation (as a
// quaternion), orient the mesh and write out a new line to vripconf
void autorient(const char *filename, const VRIPxform &v, FILE *vripconf)
{
	// Read mesh
	TriMesh *mesh = TriMesh::read(filename);
	if (!mesh) {
		fprintf(stderr, "Couldn't read %s\n", filename);
		return;
	}

	// Figure out the name of the .xf file
	const char *dot = rindex(filename, '.');
	if (!dot)
		dot = filename + strlen(filename);
	int l = dot - filename;
	char xfname[l+4];
	strncpy(xfname, filename, l);
	xfname[l] = '.';
	xfname[l+1] = 'x';
	xfname[l+2] = 'f';
	xfname[l+3] = '\0';


	// Transform the mesh to world coordinates
	xform xf = VRIP2xform(v);
	int nv = mesh->vertices.size();
	for (int i = 0; i < nv; i++)
		mesh->vertices[i] = xf * mesh->vertices[i];
	mesh->need_neighbors();
	mesh->need_adjacentfaces();

	if (mesh->normals.empty()) {
		mesh->need_normals();
		diffuse_normals(mesh, 0.5f * mesh->feature_size());
	} else {
		xform r = rot_only(xf);
		for (int i = 0; i < nv; i++)
			mesh->normals[i] = r * mesh->normals[i];
	}


	// Compute center of mass, avg normal, and bounding cone of normals
	vec avgnorm;
	point COM;
	Miniball<3,float> mb;
	int nvalid = 0;

	for (int i = 0; i < nv; i++) {
		// Ignore unconnected points and boundary vertices
		if (mesh->neighbors[i].empty())
			continue;
		if (mesh->neighbors[i].size() != mesh->adjacentfaces[i].size())
			continue;

		nvalid++;
		mb.check_in(mesh->normals[i]);
		avgnorm += mesh->normals[i];
		COM += mesh->vertices[i];
	}

	if (!nvalid) {
		printf("Bad scan... nvalid = 0\n");
		mesh->write(filename);
		xform new_xf = xform();
		new_xf.write(xfname);
		fprintf(vripconf, "bmesh %s 0 0 0 0 0 0 1\n", filename);
		fflush(vripconf);
		delete mesh;
		return;
	}

	mb.build();
	avgnorm /= nvalid;
	COM /= nvalid;

	vec centernorm;
	if (mb.squared_radius() < 0.99f) {
		printf("Using bcone\n");
		centernorm = mb.center();
	} else {
		printf("Using avg\n");
		centernorm = avgnorm;
	}
	normalize(centernorm);

	// Transform the mesh so that centernorm faces the Z axis
	vec rotaxis = centernorm CROSS vec(0,0,1);
	float rotamount = atan2(len(rotaxis), centernorm DOT vec(0,0,1));
	// If we want to flip all the way around, rotaxis will be unstable
	if (rotamount > 3.14f)
		rotaxis = vec(0,1,0);
	normalize(rotaxis);
	xform new_xf = xform::rot(rotamount, rotaxis) * xform::trans(-COM);
	for (int i = 0; i < nv; i++)
		mesh->vertices[i] = new_xf * mesh->vertices[i];
	mesh->write(filename);

	// Write out inverse transform and bmesh line
	invert(new_xf);
	new_xf.write(xfname);
	VRIPxform new_v = xform2VRIP(new_xf);
	fprintf(vripconf, "bmesh %s %f %f %f %f %f %f %f\n",
		filename, new_v[0], new_v[1], new_v[2],
		new_v[3], new_v[4], new_v[5], new_v[6]);
	fflush(vripconf);

	delete mesh;
}


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s vrip.conf\n", myname);
	exit(1);
}


int main(int argc, const char *argv[])
{
	if (argc < 2)
		usage(argv[0]);

	// Rename vrip.conf and open it
	const char *vripconfname = argv[1];
	char *vripconfbakname = new char[strlen(vripconfname)+5];
	sprintf(vripconfbakname, "%s.bak", vripconfname);

	if (rename(vripconfname, vripconfbakname) != 0)
		usage(argv[0]);

	FILE *oldvripconf = fopen(vripconfbakname, "r");
	FILE *newvripconf = fopen(vripconfname, "w");
	if (!oldvripconf || !newvripconf)
		usage(argv[0]);

	VRIPxform v;
	char mesh[1024];
	while (fscanf(oldvripconf, "bmesh %s %lf %lf %lf %lf %lf %lf %lf\n",
				   mesh, &v[0], &v[1], &v[2],
				   &v[3], &v[4], &v[5], &v[6]) == 8) {
		autorient(mesh, v, newvripconf);
	}

	fclose(oldvripconf);
	fclose(newvripconf);
}

