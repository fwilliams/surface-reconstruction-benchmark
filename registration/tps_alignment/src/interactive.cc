/*
Copyright 2007 Benedict Brown and Szymon Rusinkiewicz
Princeton University

interactive.cc
Interactive locally weighted ICP
Heavily based on plyv.cc by Szymon Rusinkiewicz

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
#include <assert.h>
#include "EdgeMesh.h"
#include "TriMesh_algo.h"
#include "XForm.h"
#include "KDtree.h"
#include "GLCamera.h"
#include "ICP.h"
// #include "tps_defs.h"
#include <GL/glut.h>
#include "glui.h"
#include <vector>
using namespace std;

// diffusion sigmas
static float d_geom = 0;
static float d_norm = 1;
static float d_curv = 1;

// maximum edge length
static float max_el_factor = -1;

// other parameters
static int num_points = 500;
static float max_kd_dist = 5;
static float gauss_dist = 20; // david; use 5 for penguin
#define MAX_COLOR_DIST 5

// mouse defines
#define LBUTTON 1<<0
#define MBUTTON 1<<1
#define RBUTTON 1<<2
#define ROTBUTTON LBUTTON
#define TRANSXYBUTTON MBUTTON
#define TRANSZBUTTON1 (LBUTTON | RBUTTON)
#define TRANSZBUTTON2 (LBUTTON | MBUTTON)
#define ICPBUTTON RBUTTON

// window variables
int main_win;
GLUI *glui_win;

// Global vars
static GLCamera camera;
static xform camxf;
static float feature_size;

// flags for what to draw
#define SRC_MESH 0
#define TGT_MESH 1
#define ORG_MESH 2
#define FLS_TYPE 0
#define HUE_TYPE 1
#define RGB_TYPE 2
#define ITEM_CURRENT 0
#define ITEM_ORIG    1
#define ITEM_LOCAL   2
#define ITEM_GLOBAL  3
#define ITEM_CURVE   4
#define ITEM_EDGES   5
#define ITEM_NONE    6

static int draw_meshes[3] = { true, true, false }; // src, tgt, orig
static int draw_icp_points = false;
static int use_indexes = 0;
static int do_affine = 0;
static int draw_type = 1; // 0 = false color, 1 = distance (hue), 2 = channels
static int hue_info = ITEM_CURRENT;
static int channel_info[3] = { ITEM_CURRENT, ITEM_ORIG, ITEM_NONE }; // what to draw in r, g, b channels
static bool draw_edges = false;
static bool draw_faces = true;

// structures with what to draw
static vector<EdgeMesh *> meshes;
typedef struct { float cur,  orig, stability, orig_stability, curvature, edge; } dist_t;
static vector<point> points;
static vector<dist_t> dists[3];
static vector< pair<int, int> > ipairs;
static vector<point> point_colors;
typedef struct { float r, g, b; } color_t;
static vector<color_t> draw_colors[3];

static bool icp_click = false; // used to distinguish clicking from motion

// feature position for weighted ICP, and direction of last ICP computation
// static bool backwards_icp = true;
static int icp_vertex_center = -1;

// Various trivial utility functions
static void cls()
{
	glDisable(GL_DITHER);
	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_NORMALIZE);
	glDisable(GL_LIGHTING);
	glDisable(GL_COLOR_MATERIAL);
	// glClearColor(0, 0, 0, 1);
	glClearColor(1, 1, 1, 1);
	glClearDepth(1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

static void need_redraw() { glutSetWindow(main_win); glutPostRedisplay(); }
// static int screenx() { return glutGet(GLenum(GLUT_WINDOW_WIDTH)); }
static int screeny() { return glutGet(GLenum(GLUT_WINDOW_HEIGHT)); }

// Find closest point to p on segment from v0 to v1
static point closest_on_segment(const point &v0, const point &v1, const point &p) {
  vec v01 = v1 - v0;
  float d = (p - v0) DOT v01;
  d /= len2(v01);
  if (d < 0.0f)
    d = 0.0f;
  else if (d > 1.0f)
    d = 1.0f;
  return v0 + d * v01;
}


// Find closest point to p on face i of mesh
static point closest_on_face(const TriMesh *mesh, int i, const point &p) {
  const TriMesh::Face &f = mesh->faces[i];
  const point &v0 = mesh->vertices[f[0]];
  const point &v1 = mesh->vertices[f[1]];
  const point &v2 = mesh->vertices[f[2]];
  vec a = v1 - v0, b = v2 - v0, p1 = p - v0, n = a CROSS b;

  float A[3][3] = { { a[0], b[0], n[0] },
                    { a[1], b[1], n[1] },
                    { a[2], b[2], n[2] } };
  float x[3] = { p1[0], p1[1], p1[2] };
  int indx[3];
  ludcmp(A, indx);
  lubksb(A, indx, x);

  if (x[0] >= 0.0f && x[1] >= 0.0f && x[0] + x[1] <= 1.0f)
    return v0 + x[0] * a + x[1] * b;

  point c01 = closest_on_segment(v0, v1, p);
  point c12 = closest_on_segment(v1, v2, p);
  point c20 = closest_on_segment(v2, v0, p);
  float d01 = dist2(c01, p);
  float d12 = dist2(c12, p);
  float d20 = dist2(c20, p);
  if (d01 < d12) {
    if (d01 < d20) return c01; else return c20;
  } else {
    if (d12 < d20) return c12; else return c20;
  }
}


// Find (good approximation to) closest point on mesh to p.
// Finds closest vertex, then checks all faces that touch it.
// If the mesh has no faces, returns the closest vertex
// Returns nearest vertex num, with nearest point in out
// Return < 0 for no match
static int closest_pt(const TriMesh *mesh, const KDtree *kd, float mdist, const point &p, point &out) {
  // float maxdist2 = sqr(dist(p, mesh->bsphere.center) + mesh->bsphere.r);
  const float *match = kd->closest_to_pt(p, mdist);
  if (!match) return -1;

  int ind = (match - (const float *) &(mesh->vertices[0][0])) / 3;
  if (ind < 0 || ind >= (int) mesh->vertices.size()) {
    fprintf(stderr, "This can't happen - bogus index\n");
    return -1;
  }

  if (mesh->faces.size() == 0)  {
    out = mesh->vertices[ind];
    return ind;
  } else assert(mesh->adjacentfaces.size());

  const vector<int> &a = mesh->adjacentfaces[ind];
  if (a.empty()) {
    fprintf(stderr, "Found match to unconnected point\n");
    out =  mesh->vertices[ind];
    return ind;
  }

  float closest_dist2 = 3.3e33;
  for (unsigned int i = 0; i < a.size(); i++) {
    point c1 = closest_on_face(mesh, a[i], p);
    float this_dist2 = dist2(c1, p);
    if (this_dist2 < closest_dist2) {
      closest_dist2 = this_dist2;
      out = c1;
    }
  }
  return ind;
}

static void update_draw_dists(void) {
  for (int mesh_num = 0; mesh_num < 3; mesh_num++) {
    if (draw_type == HUE_TYPE) {
      float maxval = 0.1;
      for (size_t i = 0; i < meshes[mesh_num]->vertices.size(); i++)
        maxval = max(maxval, ((float *) &dists[mesh_num][i])[hue_info]);
      for (size_t i = 0; i < meshes[mesh_num]->vertices.size(); i++) {
        // this is a special case to show the current alignment using nice hues
        float &r = draw_colors[mesh_num][i].r;
        float &g = draw_colors[mesh_num][i].g;
        float &b = draw_colors[mesh_num][i].b;
        if (((float *) &dists[mesh_num][i])[hue_info] < 0) {
          // nice peach color
          r = 244 / 255.0;
          g = 211 / 255.0;
          b = 153 / 255.0;
        } else {
          float a = ((float *) &dists[mesh_num][i])[hue_info] / maxval;
          float h = 4.0f * (1.0f - a);
          float s = 0.7f +  0.3f * a;
          float v = s;
          
          if (s < 0) {
            r = g = b = v;
          } else {
            h = fmodf(h, 2 * (float) M_PI);
            while (h < 0) h += 2 * (float) M_PI;
            h /= ((float) M_PI / 3.0f);
            int i = int(floorf(h));
            float f = h - i;
            float p = v * (1 - s);
            float q = v * (1 - s * f);
            float t = v * (1 - s * (1 - f));
            switch (i) {
            case 0 : r = v; g = t; b = p; break;
            case 1 : r = q; g = v; b = p; break;
            case 2 : r = p; g = v; b = t; break;
            case 3 : r = p; g = q; b = v; break;
            case 4 : r = t; g = p; b = v; break;
            default: r = v; g = p; b = q; break;
            }
          }
        }
      }
    } else {
      for (size_t i = 0; i < meshes[mesh_num]->vertices.size(); i++) {
        draw_colors[mesh_num][i].r = channel_info[0] == ITEM_NONE ? 0 :
          ((float *) &dists[mesh_num][i])[channel_info[0]];
        draw_colors[mesh_num][i].g = channel_info[1] == ITEM_NONE ? 0 :
          ((float *) &dists[mesh_num][i])[channel_info[1]];
        draw_colors[mesh_num][i].b = channel_info[2] == ITEM_NONE ? 0 :
          ((float *) &dists[mesh_num][i])[channel_info[2]];
      }
    }
  }
}

static void resetview() {
	camxf = xform::trans(0, 0, -4.0f * meshes[0]->bsphere.r) *
		xform::trans(-meshes[0]->bsphere.center);
	camera.stopspin();
}


// Draw triangle strips
static void draw_tstrips(const EdgeMesh *mesh) {
	const int *t = &mesh->tstrips[0];
	const int *end = t + mesh->tstrips.size();
	while (likely(t < end)) {
		int striplen = *t++;   
		glDrawElements(GL_TRIANGLE_STRIP, striplen, GL_UNSIGNED_INT, t);
		t += striplen;
	}
}

static void drawmesh() {
	for (size_t mesh_num = 0; mesh_num < meshes.size(); mesh_num++) {
		if (!draw_meshes[mesh_num])
			continue;

		EdgeMesh *themesh = meshes[mesh_num];

		// Set up vertex arrays
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, &themesh->vertices[0][0]);
		glEnableClientState(GL_NORMAL_ARRAY);
		glNormalPointer(GL_FLOAT, 0, &themesh->normals[0][0]);

		// Set up color
		if (draw_type != FLS_TYPE) {
			glDisable(GL_LIGHTING);
			glEnableClientState(GL_COLOR_ARRAY);
			glColorPointer(3, GL_FLOAT, 0, &draw_colors[mesh_num][0]);
		} else {
      glEnable(GL_LIGHTING);
      glDisableClientState(GL_COLOR_ARRAY);
			glColor3f(mesh_num == 0, mesh_num == 2, mesh_num == 1);
		}

		glPolygonOffset(1.0f, 1.0f);
		glEnable(GL_POLYGON_OFFSET_FILL);

    if (themesh->tstrips.empty()) {
      // No triangles - draw as points
      glPointSize(1);
      glDisable(GL_LIGHTING);
      glDrawArrays(GL_POINTS, 0, themesh->vertices.size());
    } else if (draw_faces) {
#if 1
			draw_tstrips(themesh);
      // glDrawElements(GL_TRIANGLES, themesh->faces.size(), GL_UNSIGNED_INT, &themesh->faces[0][0]);
#else
      // this code was written to work around a bug in the ATI Linux driver
      // which should be fixed now
      vector<vec> &n = themesh->normals;
      vector<point> &v = themesh->vertices;
      vector<EdgeMesh::Face> &f = themesh->faces;
      vector<color_t> &d = draw_colors[mesh_num];
      assert(n.size() == v.size());
      assert(d.size() == v.size());
      glBegin(GL_TRIANGLES);
      glColor3f(mesh_num == 0, mesh_num == 2, mesh_num == 1);
      for (unsigned int i = 0; i < f.size(); i++) {
        unsigned int a = f[i][0]; assert(a < v.size());
        unsigned int b = f[i][1]; assert(b < v.size());
        unsigned int c = f[i][2]; assert(c < v.size());
        if (draw_type != FLS_TYPE) glColor3f(d[a].r, d[a].g, d[a].b);
        glNormal3f(n[a][0], n[a][1], n[a][2]);
        glVertex3f(v[a][0], v[a][1], v[a][2]);
        if (draw_type != FLS_TYPE) glColor3f(d[b].r, d[b].g, d[b].b);
        glNormal3f(n[b][0], n[b][1], n[b][2]);
        glVertex3f(v[b][0], v[b][1], v[b][2]);
        if (draw_type != FLS_TYPE) glColor3f(d[c].r, d[c].g, d[c].b);
        glNormal3f(n[c][0], n[c][1], n[c][2]);
        glVertex3f(v[c][0], v[c][1], v[c][2]);
      }
      glEnd();
#endif
    }

		// glDisableClientState(GL_COLOR_ARRAY);
		glEnable(GL_LIGHTING);

		// Draw ICP center
		if (mesh_num == 0 && icp_vertex_center >= 0) {
			glPushMatrix();
			const point &c = meshes[0]->vertices[icp_vertex_center];
			glTranslatef(c[0], c[1], c[2]);
			glColor3f(1, 1, 1);
			glutSolidSphere(5.0 * feature_size, 10, 10);
			glPopMatrix();
		}

		// draw selected points for optimization
		if (mesh_num <= 1 && draw_icp_points) {
      while (point_colors.size() < ipairs.size())
        point_colors.push_back(point(drand48(), drand48(), drand48()));

			for (unsigned int i = 0; i < ipairs.size(); i++) {
				glPushMatrix();
        glColor3f(point_colors[i][0], point_colors[i][1], point_colors[i][2]);
        int vnum = (mesh_num == 1) ? ipairs[i].first : ipairs[i].second;
        if (vnum >= (int) meshes[mesh_num]->vertices.size()) {
          fprintf(stderr, "%d %u %d %d %u\n", vnum, (unsigned int) meshes[mesh_num]->vertices.size(),
                  (int) mesh_num, i, (unsigned int) ipairs.size());
          assert(0);
        }
        glTranslatef(meshes[mesh_num]->vertices[vnum][0],
                     meshes[mesh_num]->vertices[vnum][1],
                     meshes[mesh_num]->vertices[vnum][2]);
				glutSolidSphere(1.5, 10, 10);
				glPopMatrix();
			}
		}

    // draw edges
		if (draw_edges) {
			glDisable(GL_LIGHTING);
			glColor3f(mesh_num == 0, mesh_num == 2, mesh_num == 1);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			draw_tstrips(themesh);
			glPolygonMode(GL_FRONT, GL_FILL);
		}
	}
}

// Redraw function
static void redraw() {
	camera.setupGL(camxf * meshes[0]->bsphere.center, meshes[0]->bsphere.r);
	glLoadMatrixd(camxf);
	cls();

	GLfloat mat_diffuse[4] = { 0.9, 0.9, 0.9, 0.9 };
	GLfloat mat_specular[4] = { 0.25, 0.25, 0.25, 0.25 };
	GLfloat mat_shininess[] = { 127 };
	GLfloat global_ambient[] = { 0.05, 0.05, 0.05, 0.05 };
	GLfloat light0_ambient[] = { 0, 0, 0, 0 };
	GLfloat light0_diffuse[] = { 0.85, 0.85, 0.85, 0.85 };
	GLfloat light0_specular[] = { 0.85, 0.85, 0.85, 0.85 };
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);
	// glDisable(GL_DEPTH_TEST);

	glCullFace(GL_BACK);
	glEnable(GL_CULL_FACE);
	// glDisable(GL_CULL_FACE);

	// Draw
	drawmesh();

	glutSwapBuffers();
}


// do icp at clicked location
static void do_icp(int x, int y) {
	fprintf(stderr, "Doing ICP\n");

	// first render id buffer to figure out
  // which face was clicked on
	glDisable(GL_LIGHTING);
	glDisable(GL_CULL_FACE);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glLoadMatrixd(camxf);

	glClearColor(1, 1, 1, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glBegin(GL_TRIANGLES);
	for (size_t i = 0; i < meshes[0]->faces.size(); i++) {
		unsigned char r = (unsigned char)  i        & 255;
		unsigned char g = (unsigned char) (i >>  8) & 255;
		unsigned char b = (unsigned char) (i >> 16) & 255;
		glColor3ub(r, g, b);
		int v = meshes[0]->faces[i][0];
		glVertex3fv(meshes[0]->vertices[v]);
		v = meshes[0]->faces[i][1];
		glVertex3fv(meshes[0]->vertices[v]);
		v = meshes[0]->faces[i][2];
		glVertex3fv(meshes[0]->vertices[v]);
	}
	glEnd();

	// read back face_num
	unsigned char rgb[3];
	glReadBuffer(GL_BACK);
	glReadPixels(x, screeny() - y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, rgb);
	unsigned f = (((unsigned) rgb[2]) << 16) +
		     (((unsigned) rgb[1]) << 8) +
		       (unsigned) rgb[0];
	if (f >= meshes[0]->faces.size()) {
		fprintf(stderr, "Hit background\n");
		return;
	}

  // select the first vertex of this face as weighted ICP center
	int v = meshes[0]->faces[f][0];
	fprintf(stderr, "rgb: (%d, %d, %d), f: %d, v: %d\n",
		(int) rgb[0], (int) rgb[1], (int) rgb[2], f, v);

	icp_vertex_center = v;

	// do weighted icp alignment
	float maxdist = 50;
	int verbose = 2;

	EdgeMesh *tgt_mesh = meshes[1];
	KDtree *tgt_kd = new KDtree(tgt_mesh->vertices);
  ICP_struct tgt_icp(tgt_mesh, xform::identity(), tgt_kd, true);

	EdgeMesh *src_mesh = meshes[2];
	KDtree *src_kd = new KDtree(src_mesh->vertices);
  ICP_struct src_icp(src_mesh, xform::identity(), src_kd, true);

#define WEIGHT_SCALE 0.01f

  // compute stability weights
  compute_overlaps(tgt_icp, src_icp, maxdist, verbose);

  // find closest point to ICP center on other mesh, so we
  // can compute proper PDFs
  const point &p = meshes[0]->vertices[icp_vertex_center];
	point psrc = p;
	const float *tmp = src_kd->closest_to_pt(p, sqr(gauss_dist));
	if (tmp)
		psrc = point(tmp);
	point ptgt = p;
	tmp = tgt_kd->closest_to_pt(p, sqr(gauss_dist));
	if (tmp)
		ptgt = point(tmp);

  // weight probability of each point by it's distance to ICP center
	for (size_t j = 0; j < tgt_mesh->vertices.size(); j++) {
		if (!tgt_icp.weights[j])
			continue;
		const point &q = tgt_mesh->vertices[j];
		float sqdist = dist2(ptgt, q);
		tgt_icp.weights[j] = WEIGHT_SCALE / (sqr(gauss_dist) + sqdist);
	}

	for (size_t j = 0; j < src_mesh->vertices.size(); j++) {
		if (!src_icp.weights[j])
			continue;
		const point &q = src_mesh->vertices[j];
		float sqdist = dist2(psrc, q);
		src_icp.weights[j] = WEIGHT_SCALE / (sqr(gauss_dist) + sqdist);
	}

  float stability, max_stability;
  ipairs.clear();
  // backwards_icp = false;
  float icp_err = ICP(tgt_icp, src_icp, stability, max_stability, maxdist,
                      verbose, &points, ipairs, false, false);
	bool icp_stable = ((icp_err >= 0.0f) && (stability > 0.01f));

  cerr << src_icp.xf; // print out transform

  // apply transform
  vector<point> &v0 = meshes[0]->vertices;
  vector<point> &v2 = meshes[2]->vertices;
  for (size_t j = 0; j < v2.size(); j++)
    v0[j] = src_icp.xf * v2[j];
  assert(meshes[0]->normals.size() == meshes[2]->normals.size());
  xform rot_xf = rot_only(src_icp.xf);
  for (unsigned int j = 0; j < meshes[2]->normals.size(); j++)
    meshes[0]->normals[j] = rot_xf * meshes[2]->normals[j];
  fprintf(stderr, "ICP Stability = %f, Max Stability = %g ICP Err = %f\n",
          stability, max_stability, icp_err);

  if (!icp_stable)
		fprintf(stderr, "Unstable ICP\n");

  // update dists array
  KDtree *out_kd = new KDtree(meshes[0]->vertices);
  for (int mesh_num = 0; mesh_num < 3; mesh_num++) {
    ICP_struct &m = (mesh_num != 1) ? src_icp : tgt_icp;
    KDtree *kd = (mesh_num == 1) ? out_kd : tgt_kd;

    float max_stability = m.pdf.empty() ? 0 : *max_element(m.pdf.begin(), m.pdf.end());
    for (size_t i = 0; i < meshes[mesh_num]->vertices.size(); i++) {
      const point &p = meshes[mesh_num]->vertices[i];
      point q;
      int idx = closest_pt(mesh_num == 1 ? meshes[0] : meshes[1], kd, sqr(MAX_COLOR_DIST), p, q);
      dists[mesh_num][i].cur = (idx > 0) ? dist(p, q) / MAX_COLOR_DIST : -1;
      dists[mesh_num][i].stability = m.pdf.empty() ? 0 : m.pdf[i] / max_stability;
    }
  }

  delete out_kd;

	delete src_kd;
	delete tgt_kd;
	update_draw_dists();
	need_redraw();
}

// Mouse and keyboard GLUT callbacks
static unsigned buttonstate = 0;
static void mousemotionfunc(int x, int y) {
	icp_click = false;
	Mouse::button b = Mouse::NONE;
	if (buttonstate & (1 << 3))
		b = Mouse::WHEELUP;
	else if (buttonstate & (1 << 4))
		b = Mouse::WHEELDOWN;
	else if (buttonstate & (1 << 30))
		b = Mouse::LIGHT;
	else if (buttonstate == ROTBUTTON)
		b = Mouse::ROTATE;
	else if (buttonstate == TRANSXYBUTTON)
		b = Mouse::MOVEXY;
	else if (buttonstate == TRANSZBUTTON1)
		b = Mouse::MOVEZ;
	else if (buttonstate == TRANSZBUTTON2)
		b = Mouse::MOVEZ;

	glLoadIdentity();
	camera.mouse(x, y, b,
		     camxf * meshes[0]->bsphere.center, meshes[0]->bsphere.r,
		     camxf);
	if (b != Mouse::NONE)
		need_redraw();
}

static void mousebuttonfunc(int button, int state, int x, int y) {
	if (state == GLUT_DOWN)
		buttonstate |= (1 << button);
	else
		buttonstate &= ~(1 << button);
	if (glutGetModifiers() & GLUT_ACTIVE_SHIFT)
		buttonstate |= (1 << 30);
	else
		buttonstate &= ~(1 << 30);

	if (buttonstate == ICPBUTTON)
		icp_click = true;
	else if (buttonstate == 0 && icp_click)
		do_icp(x, y);
	else
		mousemotionfunc(x, y);
}


static void keyboardfunc(unsigned char key, int /* x */, int /* y */) {
	static int reps = 0;
	static bool write_outputs = false;
	switch (key) {
  case 'e':
    draw_edges = !draw_edges; need_redraw(); break;
  case 'f':
    draw_faces = !draw_faces; need_redraw(); break;
  case ' ':
    resetview(); need_redraw(); break;
  case 'S':
  case 's': {
    // save mesh to out.ply
    char fname[100];
    meshes[2]->write("mesha.ply");
    meshes[1]->write("meshb.ply");
    sprintf(fname, "out_%05d.ply", reps);
    meshes[0]->write(fname);
    break;
  }
  case 'W':
  case 'w':
    write_outputs = !write_outputs;
    printf("Write outputs: %d\n", write_outputs);
    break;
  case 'b': {
    // compute stable sampling of regions
    vector<int> ipoints;
    stable_sample(meshes[0], num_points, ipoints);
    points.resize(ipoints.size());
    for (unsigned int i = 0; i < ipoints.size(); i++)
      points[i] = meshes[0]->vertices[ipoints[i]];
    if (draw_icp_points) need_redraw();
  } break;
  case '\033': // Esc
  case '\021': // Ctrl-Q
  case 'Q':
  case 'q':
    exit(0);
	}
}

// needed to keep glui happy
static void reshape(int, int) {
  GLUI_Master.auto_set_viewport();
  need_redraw();
}

#if 0
// Idle callback, enable this for autospinning
static void idle() {
	if (camera.autospin(camxf))
		need_redraw();
	else
		usleep(10000);
}
#endif

static void glui_checkbox_callback(int /* id */) {
  need_redraw();
}

static void glui_draw_type_callback(int /* id */) {
  update_draw_dists();
  need_redraw();
}

static void glui_channel_callback(int id) {
  update_draw_dists();
  if (draw_type == RGB_TYPE && id < 3) need_redraw();
  else if (draw_type == HUE_TYPE && id == 3) need_redraw();
}

// prune edges and corresponding triangle longer than
// specified length (for e.g. sea floor data)
static void prune_edges(EdgeMesh *themesh, float max_length) {
  vector<point> v = themesh->vertices;
  vector<EdgeMesh::Face> f = themesh->faces;

  // keep faces with all edges short enough
  themesh->faces.clear();
  for (unsigned int i = 0; i < f.size(); i++) {
    assert(f[i][0] < (int) v.size());
    assert(f[i][1] < (int) v.size());
    assert(f[i][2] < (int) v.size());
    if (dist(v[f[i][0]], v[f[i][1]]) < max_length &&
        dist(v[f[i][1]], v[f[i][2]]) < max_length &&
        dist(v[f[i][2]], v[f[i][0]]) < max_length) {
      themesh->faces.push_back(f[i]);
    }
  }

  // now go through list of vertices on each kept faces,
  // assigning them new IDs starting at zero; that way
  // we clear out all the vertices we aren't using
  vector<int> vertex_map(v.size(), -1);
  themesh->vertices.clear();
  int maxv = -1;
  for (unsigned int i = 0; i < themesh->faces.size(); i++) {
    for (int j = 0; j < 3; j++) {
      if (vertex_map[themesh->faces[i][j]] > 0) {
        themesh->faces[i][j] = vertex_map[themesh->faces[i][j]];
      } else {
        maxv++;
        vertex_map[themesh->faces[i][j]] = maxv;
        themesh->vertices.push_back(v[themesh->faces[i][j]]);
        themesh->faces[i][j] = maxv;
      }
    }
  }

  assert(themesh->faces.size() > 0);
  fprintf(stderr, "Done pruning faces.\n");
}

// fetch the name of the mesh; either this should be a fully-qualified
// .ply file, with res == -1, or a .set file (minus the extension) with
// the resolution specified in res.
static void get_mesh_name(const char *plyfilename, char *fname, int res) {
	if (res == -1) {
		strcpy(fname, plyfilename);
		// strcat(fname, ".ply");
		return;
	}

	// extract path component of plyfilename
	puts(plyfilename);
	int path_comp = strlen(plyfilename);
	while (path_comp && plyfilename[path_comp - 1] != '/') path_comp--;
	printf("path_comp = %d\n", path_comp);

	char tmp[1024];
	sprintf(tmp, "%s.set", plyfilename);
	puts(tmp);
	FILE *f = fopen(tmp, "r");
	assert(f);

	// skip first two lines
	while (getc(f) != '\n');
	while (getc(f) != '\n');

	// now get to the right resolution
	for (int i = 0; i < res; i++) while (getc(f) != '\n');

	// read in the appropriate ply file name
	fscanf(f, "%*s %*s %s", tmp);
	strncpy(fname, plyfilename, path_comp);
	strcpy(fname + path_comp, tmp);

	fclose(f);
}


// apply rigid-body transformcontained in a .xf to a mesh
static void apply_xform(EdgeMesh *themesh, const char *plyfilename) {
	char xfname[1024];
	sprintf(xfname, "%s", plyfilename);
  if (strlen(xfname) > 4 && xfname[strlen(xfname) - 4] == '.')
    strcpy(xfname + strlen(xfname) - 3, "xf");
  else strcat(xfname, ".xf");

	xform xf;
  xf.read(xfname);

	vector<point> &v = themesh->vertices;
	for (size_t i = 0; i < v.size(); i++)
		v[i] = xf * v[i];
}


static void usage(const char *myname) {
	fprintf(stderr, "Usage: %s [ -res res1 ] mesh1 [ -res res2 ] mesh2 [ options ]\n"
          "Options: -kd <max_kd_dist>\n"
          "         -dg <diffuse_geometry_amount>\n"
          "         -dn <diffuse_normal_amount>\n"
          "         -dc <diffuse_curvature_amount>\n"
          "         -el <max_edge_length>\n"
          "         -ss (sort stability for sampling on initial alignment)\n"
          "         -af (do affine initial alignment)\n"
          "All option values are scaled by the mesh feature size.\n", myname);
	exit(1);
}

static void process_args(int argc, char **argv, char *fmesh[2], int res[2]) {
  // by default, the mesh is a mesh, not a .set
  res[0] = res[1] = -1;

  int cur_mesh = 0;

  for (int curarg = 1; curarg < argc; curarg++) {
    if (!strcmp(argv[curarg], "-res")) {
      // specify resolution of mesh
      curarg++;
      res[cur_mesh] = atoi(argv[curarg]);
      continue;
    } else if (!strcmp(argv[curarg], "-kd")) {
      // max distance for kd searches
      // mesh is not scaled if this is specified
      curarg++;
      max_kd_dist = atof(argv[curarg]);
      continue;
    } else if (!strcmp(argv[curarg], "-dg")) {
      // diffuse geometry
      curarg++;
      d_geom = atof(argv[curarg]);
      continue;
    } else if (!strcmp(argv[curarg], "-dn")) {
      // diffuse noramls
      curarg++;
      d_norm = atof(argv[curarg]);
      continue;
    } else if (!strcmp(argv[curarg], "-dc")) {
      // diffuse curvatures
      curarg++;
      d_curv = atof(argv[curarg]);
      continue;
    } else if (!strcmp(argv[curarg], "-el")) {
      // prune edges longer than this times that feature size
      curarg++;
      max_el_factor = atof(argv[curarg]);
      continue;
    } else if (!strcmp(argv[curarg], "-ss")) {
      // sort stabilities initially
      // afterwards this is controlled by a check box
      use_indexes = true;
    } else if (!strcmp(argv[curarg], "-af")) {
      // affine initial alignment
      // afterwards this is controlled by a check box
      do_affine = true;
    } else {
      // assume we want a file name
      fmesh[cur_mesh] = argv[curarg];
      cur_mesh++;
      assert(cur_mesh <= 2);
      continue;
    }
  }

  // make sure we specified two meshes
  assert(cur_mesh == 2);
}

static void read_meshes(char **fmesh, int *res) {
	for (int curarg = 0; curarg < 2; curarg++) {
    char fname[1024];

		get_mesh_name(fmesh[curarg], fname, res[curarg]);
		puts(fname);

		EdgeMesh *themesh = EdgeMesh::read_preprocessed(fname);
    if (!themesh) {
      themesh = EdgeMesh::read(fname);
      assert(!themesh->faces.empty());
      if (max_el_factor > 0) {
        themesh->normals.clear();
        prune_edges(themesh, max_el_factor * themesh->feature_size());
      }
      themesh->need_tstrips();
      themesh->need_neighbors();
      themesh->need_adjacentfaces();

      themesh->need_edges();
    }
		if (!themesh) {
			fprintf(stderr, "Couldn't open file %s\n", fname);
			exit(1);
		}

		apply_xform(themesh, fmesh[curarg]);

		meshes.push_back(themesh);
	}
}

static void setup_window(char *title) {
  // Create the window, register callbacks, go for it
  main_win = glutCreateWindow(title);

  glutDisplayFunc(redraw);
  GLUI_Master.set_glutMouseFunc(mousebuttonfunc);
  // GLUI_Master.set_glutMotionFunc(mousemotionfunc);
  glutMotionFunc(mousemotionfunc);
  GLUI_Master.set_glutKeyboardFunc(keyboardfunc);
  GLUI_Master.set_glutReshapeFunc(reshape);
  // GLUI_Master::set_glutIdleFunc(idle); // commented out to prevent spinning

  glui_win = GLUI_Master.create_glui_subwindow(main_win, GLUI_SUBWINDOW_LEFT);
  glui_win->set_main_gfx_window(main_win);

  glui_win->add_checkbox("Source",   &draw_meshes[0], 0, glui_checkbox_callback);
  glui_win->add_checkbox("Target",   &draw_meshes[1], 1, glui_checkbox_callback);
  glui_win->add_checkbox("Original", &draw_meshes[2], 2, glui_checkbox_callback);

  glui_win->add_separator();

  glui_win->add_checkbox("Draw ICP Points", &draw_icp_points, 3, glui_checkbox_callback);

  glui_win->add_separator();

  // this ensures the command-line options only affect the initial alignment
  use_indexes = false;
  do_affine = false;
  glui_win->add_checkbox("Affine ICP", &do_affine, 4, glui_checkbox_callback);
  glui_win->add_checkbox("Sort Stabilities", &use_indexes, 5, glui_checkbox_callback);

  glui_win->add_separator();

  GLUI_Listbox *rgbl[4];
  GLUI_RadioGroup *color_group = glui_win->add_radiogroup(&draw_type, 0, glui_draw_type_callback);
  glui_win->add_radiobutton_to_group(color_group, "False Color");
  glui_win->add_radiobutton_to_group(color_group, "Hue");
  rgbl[3] = glui_win->add_listbox("Hue:",   &hue_info, 3, glui_channel_callback);

  glui_win->add_radiobutton_to_group(color_group, "Channels");

  rgbl[0] = glui_win->add_listbox("Red:",   &channel_info[0], 0, glui_channel_callback);
  rgbl[1] = glui_win->add_listbox("Green:", &channel_info[1], 1, glui_channel_callback);
  rgbl[2] = glui_win->add_listbox("Blue:",  &channel_info[2], 2, glui_channel_callback);

  for (int i = 0; i < 4; i++) {
    rgbl[i]->add_item(ITEM_LOCAL,    "Local Stability");
    rgbl[i]->add_item(ITEM_GLOBAL,   "Global Stability");
    rgbl[i]->add_item(ITEM_CURRENT , "Current Distance");
    rgbl[i]->add_item(ITEM_ORIG,     "Original Distance");
    rgbl[i]->add_item(ITEM_CURVE,    "Curvature");
    rgbl[i]->add_item(ITEM_EDGES,    "Edges");
    rgbl[i]->add_item(ITEM_NONE,     "None");
  }

  glui_win->sync_live();

  resetview();
}

int main(int argc, char *argv[]) {

  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInit(&argc, argv);

  char *fmesh[2];
  int  res[2];
  process_args(argc, argv, fmesh, res);

	// Read in .ply file
	if (argc < 2)
		usage(argv[0]);

  EdgeMesh::set_verbose(0);
  read_meshes(fmesh, res);

	printf("num_points: %d\nmax_kd_dist: %g\n\n",	num_points, max_kd_dist);

	if (int(meshes[0]->vertices.size()) < num_points)
		num_points = int(meshes[0]->vertices.size());

  // EdgeMesh::BBox b = union_bbox();

	for (size_t mesh = 0; mesh < meshes.size(); mesh++) {
		EdgeMesh *themesh = meshes[mesh];

		// Do all the smoothing, then compute tstrips and bsphere, ...
    fprintf(stderr, "Mesh feature size: %f\n", themesh->feature_size());
		themesh->confidences.clear();
		themesh->colors.clear();

    if (d_geom > 0) {
      themesh->normals.clear();
      smooth_mesh(themesh, d_geom * themesh->feature_size());
    }

    if (themesh->normals.empty()) {
      themesh->need_normals();
      if (d_norm > 0) diffuse_normals(themesh, d_norm * themesh->feature_size());
      themesh->need_curvatures();
      if (d_curv > 0) diffuse_curv(themesh, d_curv * themesh->feature_size());
    }

		if (themesh->faces.size()) {
      themesh->tstrips.clear();
      themesh->need_tstrips();
    }
		themesh->need_bsphere();
	}
	feature_size = meshes[0]->feature_size();

  // compute initial (unweighted) alignemnt
  float stability, max_stability;
  vector<float> src_weights, tgt_weights, out_pdf_1, out_pdf_2;
  KDtree *src_kd = new KDtree(meshes[0]->vertices);
  KDtree *tgt_kd = new KDtree(meshes[1]->vertices);
  // backwards_icp = false;

  ICP_struct m0(meshes[0], xform::identity(), src_kd, true);
  ICP_struct m1(meshes[1], xform::identity(), tgt_kd, true);

  // use_indexes, do_affine
  // set true, false for FUR,
  // false, false for most other stuff
  float icp_err = ICP(m1, m0, stability, max_stability, 15, 2,
                      &points, ipairs, use_indexes, do_affine);
  // assert(out_pdf.size() == meshes[0]->vertices.size());
  if ((m0.pdf.size() != meshes[0]->vertices.size()) && (m1.pdf.size() > 0)) {
    fprintf(stderr, "%u %u %u\n", (unsigned int) m0.pdf.size(),
            (unsigned int) meshes[0]->vertices.size(), (unsigned int) meshes[1]->vertices.size());
    assert(0);
  }

  fprintf(stderr, "ICP Stability = %f, Max Stability = %g ICP Err = %f\n",
          stability, max_stability, icp_err);
	bool icp_stable = (icp_err >= 0.0f && stability > 0.01f);
  cerr << m0.xf;

	if (!icp_stable) {
		fprintf(stderr, "Unstable pair by correspond standards.\n");
	}	else {} /* stupid trick to always transform regardless */  {
		// Transform mesh1 according to xf1
		vector<point> &v0 = meshes[0]->vertices;
		for (size_t i = 0; i < v0.size(); i++)
			v0[i] = m0.xf * v0[i];
    delete src_kd;
    src_kd = new KDtree(meshes[0]->vertices);
	}

  // init dist arrays
  dists[0].resize(meshes[0]->vertices.size());
  draw_colors[0].resize(meshes[0]->vertices.size());
  dists[1].resize(meshes[1]->vertices.size());
  draw_colors[1].resize(meshes[1]->vertices.size());
  dists[2].resize(meshes[0]->vertices.size());
  draw_colors[2].resize(meshes[0]->vertices.size());

  for (unsigned int i = 0; i < meshes[0]->vertices.size(); i++) {
    point p;
    int idx = closest_pt(meshes[1], tgt_kd, sqr(MAX_COLOR_DIST), meshes[0]->vertices[i], p);
    dists[0][i].cur = dists[0][i].orig = dists[2][i].cur = dists[2][i].orig =
      idx > 0 ? dist(meshes[0]->vertices[i], p) / MAX_COLOR_DIST : -1;
  }

  for (unsigned int i = 0; i < meshes[1]->vertices.size(); i++) {
    point p;
    int idx = closest_pt(meshes[0], src_kd, sqr(MAX_COLOR_DIST), meshes[1]->vertices[i], p);
    dists[1][i].cur = dists[1][i].orig =
      idx > 0 ? dist(meshes[1]->vertices[i], point(p)) / MAX_COLOR_DIST : -1;
  }

  // update dists array with stability and curvature info
  if (m0.pdf.size() > 0) {
    float max_stability = *max_element(m0.pdf.begin(), m0.pdf.end());
    float max_curve = 0;
    for (size_t i = 0; i < meshes[0]->vertices.size(); i++) {
      float curve = meshes[0]->curv1[i] + meshes[0]->curv2[i];
      if (max_curve < curve) max_curve = curve;
    } max_curve = 1 / max_curve;
    for (size_t i = 0; i < meshes[0]->vertices.size(); i++) {
      dists[0][i].orig_stability = dists[0][i].stability =
        dists[2][i].orig_stability = dists[2][i].stability = m0.pdf[i] / max_stability;
      dists[0][i].curvature  = dists[2][i].curvature  = powf((meshes[0]->curv1[i] + meshes[0]->curv2[i]) * max_curve, 0.5f);
      dists[0][i].edge = dists[2][i].edge = meshes[0]->isedge(i);
    }
  } else {
    // for (size_t i = 0; i < meshes[0]->vertices.size(); i++)
      // dists[0][i].change = dists[2][i].change = 0;
  }

  if (m1.pdf.size() > 0) {
    float max_stability = *max_element(m1.pdf.begin(), m1.pdf.end());
    float max_curve = 0;
    for (size_t i = 0; i < meshes[1]->vertices.size(); i++) {
      float curve = meshes[1]->curv1[i] + meshes[1]->curv2[i];
      if (max_curve < curve) max_curve = curve;
    } max_curve = 1 / max_curve;
    for (size_t i = 0; i < meshes[1]->vertices.size(); i++) {
      dists[1][i].orig_stability = dists[1][i].stability = m1.pdf[i] / max_stability;
      dists[1][i].curvature  = powf((meshes[1]->curv1[i] + meshes[1]->curv2[i]) * max_curve, 0.5f);
      dists[1][i].edge = meshes[0]->isedge(i);
    }
  } else {
    // for (size_t i = 0; i < meshes[1]->vertices.size(); dists[2][i++].change = 0);
  }

  delete src_kd;
  delete tgt_kd;

	// Create new mesh: meshes[2] starts out as a copy of meshes[0]
	meshes.push_back(new EdgeMesh(*meshes[0]));

  update_draw_dists();

  setup_window(argv[0]);

  glutMainLoop();
}
