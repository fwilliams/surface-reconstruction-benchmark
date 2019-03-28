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

#include "RotationMatrix.h"

RotationMatrix::RotationMatrix()  {
	rotation_matrix = new DenseMatrix(3);
}

RotationMatrix::RotationMatrix(Vector3 _axis, double _angle)  {
	axis = _axis;
	angle = _angle;

	DenseMatrix* outer_product = new DenseMatrix(axis.x*axis.x, axis.x*axis.y, axis.x*axis.z, axis.y*axis.y,
			axis.y*axis.z, axis.z*axis.z);

	DenseMatrix* orthogonal_subspace = new DenseMatrix(3);
	orthogonal_subspace->sub(outer_product);
	orthogonal_subspace->scale(cos(angle));

	DenseMatrix* skew_symmetric = new DenseMatrix(3);
	skew_symmetric->setEntry(0,0, 0); skew_symmetric->setEntry(1,1, 0); skew_symmetric->setEntry(2,2, 0);
	skew_symmetric->setEntry(1,0, -axis.z); skew_symmetric->setEntry(0,1, axis.z);
	skew_symmetric->setEntry(2,0, axis.y); skew_symmetric->setEntry(0,2, -axis.y);
	skew_symmetric->setEntry(2,1, -axis.x); skew_symmetric->setEntry(1,2, axis.x);
	skew_symmetric->scale(sin(angle));

	rotation_matrix = new DenseMatrix(outer_product);
	rotation_matrix->add(orthogonal_subspace);
	rotation_matrix->add(skew_symmetric);

	delete outer_product;
	delete orthogonal_subspace;
	delete skew_symmetric;
}

RotationMatrix::~RotationMatrix()  {
	delete rotation_matrix;
}

Vector3 RotationMatrix::rotate(Vector3 _vector)  {
	return rotation_matrix->mult(_vector);
}

void RotationMatrix::printMatrix()  {
	rotation_matrix->printMatrix();
}
