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

#include "Rigid.h"

Rigid::Rigid()  {
	rotation = new DenseMatrix(3);
	translation = Vector3(0,0,0);
}

Rigid::Rigid(string _filename)  {
	FILE* file_in = fopen(_filename.c_str(), "rb");
	size_t num_read;

	double flattened_matrix[9];
	double flattened_tran[3];
	num_read = fread(&flattened_matrix, sizeof(double), 9, file_in);
	num_read = fread(&flattened_tran, sizeof(double), 3, file_in);

	rotation = new DenseMatrix(3);
	rotation->setEntry(0,0, flattened_matrix[0]) ; rotation->setEntry(1,0, flattened_matrix[1]) ; rotation->setEntry(2,0, flattened_matrix[2]);
	rotation->setEntry(0,1, flattened_matrix[3]) ; rotation->setEntry(1,1, flattened_matrix[4]) ; rotation->setEntry(2,1, flattened_matrix[5]);
	rotation->setEntry(0,2, flattened_matrix[6]) ; rotation->setEntry(1,2, flattened_matrix[7]) ; rotation->setEntry(2,2, flattened_matrix[8]);
	translation = Vector3(flattened_tran[0], flattened_tran[1], flattened_tran[2]);

	fclose(file_in);
}

Rigid::Rigid(vector<Vector3> _srcPoints, vector<Vector3> _dstPoints)  {
	rotation = new DenseMatrix(3);
	this->solve_rigid(_srcPoints, _dstPoints, rotation, translation);
}

Rigid::~Rigid()  {
	delete rotation;
}

void Rigid::update_rigid(vector<Vector3> _srcPoints, vector<Vector3> _dstPoints)  {
	DenseMatrix new_rotation(3);
	Vector3 new_translation;

	this->solve_rigid(_srcPoints, _dstPoints, &new_rotation, new_translation);

	rotation = new_rotation.mult(rotation);
	translation = new_rotation.mult(translation) + new_translation;
}

Vector3 Rigid::apply_transformation(Vector3 _pt)  {
	return rotation->mult(_pt) + translation;
}

Vector3 Rigid::apply_rotation(Vector3 _vec)  {
	return rotation->mult(_vec);
}

void Rigid::solve_rigid(vector<Vector3> _srcPoints, vector<Vector3> _dstPoints, DenseMatrix* rot, Vector3& tran)  {
	// compute centroids
	unsigned size = _srcPoints.size();
	Vector3 src_centroid(0,0,0);
	Vector3 dst_centroid(0,0,0);
	for(unsigned p = 0; p < size; p++)  {
		src_centroid = src_centroid + _srcPoints[p];
		dst_centroid = dst_centroid + _dstPoints[p];
	}
	src_centroid.scale(1.0 / (double)size);
	dst_centroid.scale(1.0 / (double)size);

	// construct covariance matrix
	double c_xx, c_xy, c_xz, c_yx, c_yy, c_yz, c_zx, c_zy, c_zz;
	c_xx = c_xy = c_xz = c_yx = c_yy = c_yz = c_zx = c_zy = c_zz = 0;
	for(unsigned p = 0; p < size; p++)  {
		Vector3 src_pt = _srcPoints[p]-src_centroid;
		Vector3 dst_pt = _dstPoints[p]-dst_centroid;
		c_xx += src_pt.x*dst_pt.x; c_xy += src_pt.x*dst_pt.y; c_xz += src_pt.x*dst_pt.z;
		c_yx += src_pt.y*dst_pt.x; c_yy += src_pt.y*dst_pt.y; c_yz += src_pt.y*dst_pt.z;
		c_zx += src_pt.z*dst_pt.x; c_zy += src_pt.z*dst_pt.y; c_zz += src_pt.z*dst_pt.z;
	}
	DenseMatrix covariance(3);
	covariance.setEntry(0, 0, c_xx); covariance.setEntry(1, 0, c_xy); covariance.setEntry(2, 0, c_xz);
	covariance.setEntry(0, 1, c_yx); covariance.setEntry(1, 1, c_yy); covariance.setEntry(2, 1, c_yz);
	covariance.setEntry(0, 2, c_zx); covariance.setEntry(1, 2, c_zy); covariance.setEntry(2, 2, c_zz);

	DenseMatrix svd_u(3);
	DenseMatrix svd_vt(3);
	double singular_values[3];
	covariance.svd(&svd_u, singular_values, &svd_vt);

	DenseMatrix* svd_v = svd_vt.transpose();
	DenseMatrix* svd_ut = svd_u.transpose();

	DenseMatrix* vu_t = svd_v->mult(svd_ut);
	DenseEigenSystem eigen_vu(vu_t, 0, 2);
	eigen_vu.solve();
	double determinant = eigen_vu.getEigenvalue(0)*eigen_vu.getEigenvalue(1)*eigen_vu.getEigenvalue(2);
	double det_sign = determinant < 0 ? -1 : 1;
	svd_ut->setEntry(0, 2, det_sign*svd_ut->getEntry(0,2));
	svd_ut->setEntry(1, 2, det_sign*svd_ut->getEntry(1,2));
	svd_ut->setEntry(2, 2, det_sign*svd_ut->getEntry(2,2));

	DenseMatrix* new_rot = svd_v->mult(svd_ut);
	for(int i = 0; i < 3; i++)  {
		for(int j = 0; j < 3; j++)
			rot->setEntry(i,j,new_rot->getEntry(i,j));
	}
	tran = dst_centroid - rot->mult(src_centroid);

	delete new_rot;
	delete svd_v;
	delete svd_ut;
	delete vu_t;
}
