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

#ifndef DENSEMATRIX_H
#define DENSEMATRIX_H

#include "Vector.h"
#include "Vector3.h"
#include "linalg_consts.h"

#include <string>

using namespace std;

class DenseMatrix  {

	public:
		DenseMatrix();
		DenseMatrix(int _square);
		DenseMatrix(int _rows, int _columns);

		DenseMatrix(double _m00, double _m01, double _m11);

		DenseMatrix(double _m00, double _m01, double _m02,
										 double _m11, double _m12,
														  double _m22);

		DenseMatrix(double _m00, double _m01, double _m02, double _m03,
										 double _m11, double _m12, double _m13,
														  double _m22, double _m23,
														               double _m33);

		DenseMatrix(int* _entries, int _rows, int _columns);
		DenseMatrix(double* _entries, int _rows, int _columns);
		DenseMatrix(BLinAlg::Vector** _rowVectors, int _numRows);
		DenseMatrix(DenseMatrix* _matrix);

		~DenseMatrix();

		int rows()  {
			return R;
		}

		int columns()  {
			return C;
		}

		// i: column, j: row
		int flatten(int _i, int _j)  {
			return _i + C*_j;
		}

		void setEntry(int _i, int _j, double _entry)  {
			int ind = flatten(_i, _j);
			entries[ind] = _entry;
		}

		double getEntry(int _i, int _j)  {
			return entries[flatten(_i, _j)];
		}

		void setDiagonal(double *_entries);

		BLinAlg::Vector* getColumnVector(int _i);

		BLinAlg::Vector* getRowVector(int _j);

		DenseMatrix* transpose();
		DenseMatrix* normalMatrix();

		void pivotRow(int _row1, int _row2);

		void add(DenseMatrix *_mat);
		void sub(DenseMatrix *_mat);
		void scale(double _scale);
		DenseMatrix* mult(DenseMatrix *_rMatrix);
		BLinAlg::Vector* mult(BLinAlg::Vector *_vec);
		Vector3 mult(Vector3 _vec);

		DenseMatrix* covarianceMatrix(BLinAlg::Vector *_expectedVector);
		BLinAlg::Vector* expectedVector();

		double* formLAPACKMatrix();

		void luFactorization(DenseMatrix *_lMatrix, DenseMatrix *_uMatrix);
		BLinAlg::Vector* linearSystem(BLinAlg::Vector* _rhs);
		BLinAlg::Vector* linearLeastSquares(BLinAlg::Vector* _rhs);
		void svd(DenseMatrix *_uMatrix, double *_singularValues, DenseMatrix *_vMatrix);

		DenseMatrix* inverse();

		void printMatrix();

	private:
		int R;
		int C;
		double* entries;

};

#endif
