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

#include "DenseMatrix.h"

extern "C"  {
	// LU factorization
	void FORTRANIZE(dgetrf)(
			int* rows, int* columns, double* matrix, int* columns2, int* pivot, int* info
	);

	// linear system
	void FORTRANIZE(dgesv)(
			int* size, int* singleRHS, double* system, int* size2, int* pivots, double* rhs, int* size3, int* info
	);

	// linear least squares by way of QR factorization
	void FORTRANIZE(dgels)(
			char* trans, int* rows, int* columns, int* singleRHS, double* system, int* rows2, double* rhs, int* rows3,
			double* workspace, int* workspaceDim, int* info
	);

	// linear least squares by way of svd
	void FORTRANIZE(dgelsd)(
			int* rows, int* columns, int* singleRHS, double* system, int* rows2, double* rhs, int* rows3,
			double* singularValues, double* cond, int* rank, double* workspace, int* workspaceDim, int* iwork, int* info
	);

	// SVD
	void FORTRANIZE(dgesvd)(
			char* U, char* V, int* rows, int* columns, double* matrix, int* rows2, double* singularValues, double* uMatrix, int* rows3, double* vMatrix, int* columns2,
			double* workspace, int* workspaceDim, int* info
	);
}

DenseMatrix::DenseMatrix()  {
	R = 1;
	C = 1;
	entries = new double[R*C];
	entries[0] = 1.0;
}

DenseMatrix::DenseMatrix(int _square)  {
	R = _square;
	C = _square;
	entries = new double[R*C];

	for(int i = 0; i < (R*C); i++)
		entries[i] = 0.0;

	double *identity = new double[_square];
	for(int i = 0; i < _square; i++)
		identity[i] = 1.0;
	setDiagonal(identity);
	delete [] identity;
}

DenseMatrix::DenseMatrix(int _rows, int _columns)  {
	R = _rows;
	C = _columns;
	entries = new double[R*C];
	for(int i = 0; i < (R*C); i++)
		entries[i] = 0.0;

	int dim = R < C ? R : C;
	double *identity = new double[dim];
	for(int i = 0; i < dim; i++)
		identity[i] = 1.0;
	setDiagonal(identity);
	delete [] identity;
}

DenseMatrix::DenseMatrix(int *_entries, int _rows, int _columns)  {
	R = _rows;
	C = _columns;
	entries = new double[R*C];

	for(int i = 0; i < (R*C); i++)
		entries[i] = _entries[i];
}

// shorthand for symmetric 2x2 matrix
DenseMatrix::DenseMatrix(double _m00, double _m01, double _m11)  {
	R = 2;
	C = 2;
	entries = new double[R*C];

	this->setEntry(0,0,_m00);
	this->setEntry(0,1,_m01);
	this->setEntry(1,0,_m01);
	this->setEntry(1,1,_m11);
}

// shorthand for symmetric 3x3 matrix
DenseMatrix::DenseMatrix(double _m00, double _m01, double _m02, double _m11, double _m12, double _m22)  {
	R = 3;
	C = 3;
	entries = new double[R*C];

	this->setEntry(0,0,_m00);
	this->setEntry(1,0,_m01); this->setEntry(0,1,_m01);
	this->setEntry(2,0,_m02); this->setEntry(0,2,_m02);
	this->setEntry(1,1,_m11);
	this->setEntry(2,1,_m12); this->setEntry(1,2,_m12);
	this->setEntry(2,2,_m22);
}

// shorthand for symmetric 4x4 matrix
DenseMatrix::DenseMatrix(double _m00, double _m01, double _m02, double _m03,
												  double _m11, double _m12, double _m13,
																   double _m22, double _m23,
																					 double _m33)  {
	R = 4;
	C = 4;
	entries = new double[R*C];

	this->setEntry(0,0,_m00);
	this->setEntry(1,0,_m01); this->setEntry(0,1,_m01);
	this->setEntry(2,0,_m02); this->setEntry(0,2,_m02);
	this->setEntry(3,0,_m03); this->setEntry(0,3,_m03);

	this->setEntry(1,1,_m11);
	this->setEntry(2,1,_m12); this->setEntry(1,2,_m12);
	this->setEntry(3,1,_m13); this->setEntry(1,3,_m13);

	this->setEntry(2,2,_m22);
	this->setEntry(3,2,_m23); this->setEntry(2,3,_m23);

	this->setEntry(3,3,_m33);
}

DenseMatrix::DenseMatrix(double *_entries, int _rows, int _columns)  {
	R = _rows;
	C = _columns;
	entries = new double[R*C];

	for(int i = 0; i < (R*C); i++)
		entries[i] = _entries[i];
}

DenseMatrix::DenseMatrix(BLinAlg::Vector **_rowVectors, int _rows)  {
	R = _rows;
	C = _rowVectors[R-1]->dim();
	entries = new double[R*C];

	for(int i = 0; i < R; i++)  {
		BLinAlg::Vector* rowVector = _rowVectors[i];
		for(int j = 0; j < C; j++)
			setEntry(j, i, rowVector->getEntry(j));
	}
}

DenseMatrix::DenseMatrix(DenseMatrix *_matrix)  {
	R = _matrix->rows();
	C = _matrix->columns();
	entries = new double[R*C];

	for(int j = 0; j < R; j++)  {
		for(int i = 0; i < C; i++)  {
			int nextInd = flatten(i, j);
			entries[nextInd] = _matrix->getEntry(i, j);
		}
	}
}

DenseMatrix::~DenseMatrix()  {
	delete [] entries;
}

void DenseMatrix::setDiagonal(double *_entries)  {
	int dim = R < C ? R : C;
	for(int i = 0; i < dim; i++)  {
		int ind = flatten(i, i);
		entries[ind] = _entries[i];
	}
}

BLinAlg::Vector* DenseMatrix::getColumnVector(int _i)  {
	BLinAlg::Vector* columnVector = new BLinAlg::Vector(R);
	for(int j = 0; j < R; j++)
		columnVector->setEntry(j, getEntry(_i, j));
	return columnVector;
}

BLinAlg::Vector* DenseMatrix::getRowVector(int _j)  {
	BLinAlg::Vector* rowVector = new BLinAlg::Vector(C);
	for(int i = 0; i < C; i++)
		rowVector->setEntry(i, getEntry(i, _j));
	return rowVector;
}

DenseMatrix* DenseMatrix::transpose()  {
	double *transpose_entries = new double[C*R];
	for(int j = 0; j < R; j++)  {
		for(int i = 0; i < C; i++)  {
			double entry = getEntry(i, j);
			int trans_ind = j + R*i;
			transpose_entries[trans_ind] = entry;
		}
	}

	DenseMatrix* transposeMatrix = new DenseMatrix(transpose_entries, C, R);
	delete [] transpose_entries;

	return transposeMatrix;
}

DenseMatrix* DenseMatrix::normalMatrix()  {
	DenseMatrix* transposed = this->transpose();
	DenseMatrix* normalMat = transposed->mult(this);

	delete transposed;

	return normalMat;
}

void DenseMatrix::pivotRow(int _row1, int _row2)  {
	BLinAlg::Vector *row1Vector = getRowVector(_row1);
	BLinAlg::Vector *row2Vector = getRowVector(_row2);

	for(int i = 0; i < C; i++)  {
		setEntry(i, _row1, row2Vector->getEntry(i));
		setEntry(i, _row2, row1Vector->getEntry(i));
	}

	delete row1Vector;
	delete row2Vector;
}

void DenseMatrix::add(DenseMatrix *_mat)  {
	if(R != _mat->rows() || C != _mat->columns())  {
		cerr << "unable to add matrix: not of same dimensions\n";
		return;
	}

	for(int j = 0; j < R; j++)  {
		for(int i = 0; i < C; i++)  {
			int nextInd = flatten(i, j);
			entries[nextInd] += _mat->getEntry(i, j);
		}
	}
}

void DenseMatrix::sub(DenseMatrix *_mat)  {
	if(R != _mat->rows() || C != _mat->columns())  {
		cerr << "unable to subtract matrix: not of same dimensions\n";
		return;
	}

	for(int j = 0; j < R; j++)  {
		for(int i = 0; i < C; i++)  {
			int nextInd = flatten(i, j);
			entries[nextInd] -= _mat->getEntry(i, j);
		}
	}
}

void DenseMatrix::scale(double _scale)  {
	for(int i = 0; i < (C*R); i++)
		entries[i] *= _scale;
}

DenseMatrix* DenseMatrix::mult(DenseMatrix *_rMatrix)  {
	if(C != _rMatrix->rows())  {
		cerr << "unable to perform matrix multiplication: mismatched rows/columns\n";
		return new DenseMatrix();
	}

	double *mult_entries = new double[R * _rMatrix->columns()];
	int ind = 0;

	for(int j = 0; j < R; j++)  {
		for(int i = 0; i < _rMatrix->columns(); i++)  {
			double innerProduct = 0.0;
			for(int k = 0; k < _rMatrix->rows(); k++)
				innerProduct += (getEntry(k, j) * _rMatrix->getEntry(i, k));
			mult_entries[ind] = innerProduct;
			ind++;
		}
	}

	DenseMatrix* result = new DenseMatrix(mult_entries, R, _rMatrix->columns());
	delete [] mult_entries;

	return result;
}

BLinAlg::Vector* DenseMatrix::mult(BLinAlg::Vector *_vec)  {
	if(_vec->dim() != C)  {
		cerr << "unable to perform matrix-vector multiplication: mismatched rows/columns\n";
		return new BLinAlg::Vector();
	}

	BLinAlg::Vector* dotResult = new BLinAlg::Vector(R);
	for(int j = 0; j < R; j++)  {
		double innerProduct = 0.0;
		for(int i = 0; i < C; i++)
			innerProduct += (getEntry(i, j) * _vec->getEntry(i));
		dotResult->setEntry(j, innerProduct);
	}

	return dotResult;
}

Vector3 DenseMatrix::mult(Vector3 _vec)  {
	if(C != 3 || R != 3)  {
		cerr << "unable to perform matrix-vector multiplication: mismatched rows/columns\n";
		return Vector3();
	}

	double arrayed[] = { _vec.x, _vec.y, _vec.z };
	double result[] = { 0,0,0 };

	for(int j = 0; j < 3; j++)  {
		double innerProduct = 0.0;
		for(int i = 0; i < 3; i++)
			innerProduct += (getEntry(i, j) * arrayed[i]);
		result[j] = innerProduct;
	}

	return Vector3(result);
}

DenseMatrix* DenseMatrix::covarianceMatrix(BLinAlg::Vector *_expectedVector)  {
	DenseMatrix *expectedMatrix = new DenseMatrix(R, C);
	for(int j = 0; j < R; j++)  {
		for(int i = 0; i < C; i++)
			expectedMatrix->setEntry(i, j, _expectedVector->getEntry(j));
	}

	DenseMatrix *centeredMatrix = new DenseMatrix(entries, R, C);
	centeredMatrix->sub(expectedMatrix);
	DenseMatrix *covarianceMatrix = new DenseMatrix(R, R);
	for(int i = 0; i < R; i++)  {
		for(int j = i; j < R; j++)  {
			double innerProduct = 0.0;
			for(int k = 0; k < C; k++)
				innerProduct += (centeredMatrix->getEntry(k, i) * centeredMatrix->getEntry(k, j));
			covarianceMatrix->setEntry(j, i, innerProduct);
		}
	}
	for(int i = 0; i < R; i++)  {
		for(int j = 0; j < i; j++)
			covarianceMatrix->setEntry(j, i, covarianceMatrix->getEntry(i, j));
	}
	covarianceMatrix->scale(1.0 / (double)C);

	delete expectedMatrix;
	delete centeredMatrix;

	return covarianceMatrix;
}

BLinAlg::Vector* DenseMatrix::expectedVector()  {
	BLinAlg::Vector *expected = new BLinAlg::Vector(R);
	for(int i = 0; i < C; i++)  {
		BLinAlg::Vector *columnVector = getColumnVector(i);
		expected->add(columnVector);
		delete columnVector;
	}

	double invScale = 1.0 / (double)C;
	expected->scale(invScale);
	return expected;
}

double* DenseMatrix::formLAPACKMatrix()  {
	double *lapack_matrix = new double[R*C];
	for(int j = 0; j < R; j++)  {
		for(int i = 0; i < C; i++)  {
			double entry = getEntry(i, j);
			int trans_ind = j + R*i;
			lapack_matrix[trans_ind] = entry;
		}
	}
	return lapack_matrix;
}

void DenseMatrix::luFactorization(DenseMatrix *_lMatrix, DenseMatrix *_uMatrix)  {
	if(R != C)  {
		cerr << "unable to perform factorization on a non-square matrix\n";
		return;
	}

	double *lu_matrix = formLAPACKMatrix();
	int *pivot = new int[C], info = 10;

	FORTRANIZE(dgetrf)(&R, &C, lu_matrix, &C, pivot, &info);
	DenseMatrix *lapack_matrix = new DenseMatrix(lu_matrix, R, C);
	DenseMatrix *full_matrix = lapack_matrix->transpose();

	for(int j = 0; j < R; j++)  {
		for(int i = 0; i < C; i++)  {
			double nextEntry = full_matrix->getEntry(i, j);

			if(i == j)  {
				_lMatrix->setEntry(i, j, 1.0);
				_uMatrix->setEntry(i, j, nextEntry);
			}
			else if(i < j)  {
				_lMatrix->setEntry(i, j, nextEntry);
				_uMatrix->setEntry(i, j, 0.0);
			}
			else  {
				_lMatrix->setEntry(i, j, 0.0);
				_uMatrix->setEntry(i, j, nextEntry);
			}
		}
	}

	for(int i = 0; i < R; i++)
		_lMatrix->pivotRow(i, pivot[i]-1);

	delete lapack_matrix;
	delete full_matrix;
	delete [] lu_matrix;
	delete [] pivot;
}

BLinAlg::Vector* DenseMatrix::linearSystem(BLinAlg::Vector* _rhs)  {
	if(R != C)  {
		cerr << "Linear system solve assumes a square matrix!\n";
		return new BLinAlg::Vector(R);
	}

	int singleRHS = 1;
	double* lapack_matrix = formLAPACKMatrix();
	int* pivotArray = new int[R];
	double* rhs_data = _rhs->dataCopy();
	int info;

	FORTRANIZE(dgesv)(&R, &singleRHS, lapack_matrix, &R, pivotArray, rhs_data, &R, &info);

	if(info != 0)
		cerr << "Error on linear solver." << endl;

	BLinAlg::Vector* solution = new BLinAlg::Vector(R);
	for(int i = 0; i < R; i++)
		solution->setEntry(i, rhs_data[i]);

	delete [] lapack_matrix;
	delete [] pivotArray;
	delete [] rhs_data;

	return solution;
}

BLinAlg::Vector* DenseMatrix::linearLeastSquares(BLinAlg::Vector* _rhs)  {
	char trans = 'N';
	int singleRHS = 1;
	double* lapack_matrix = formLAPACKMatrix();

	double* rhs_data = _rhs->dataCopy();
	double* singularValues = new double[C];
	double cond = 1e-10;
	int rank = 1;

	int workspace_size = -1;
	double* workspace = new double[1];
	int* iwork = new int[3*C*20 + 11*C];
	int info = 10;

	// first query for workspace size
	FORTRANIZE(dgelsd)(&R, &C, &singleRHS, lapack_matrix, &R, rhs_data, &R,
			singularValues, &cond, &rank, workspace, &workspace_size, iwork, &info);

	workspace_size = (int)workspace[0];
	delete [] workspace;
	workspace = new double[workspace_size];

	// then solve it
	FORTRANIZE(dgelsd)(&R, &C, &singleRHS, lapack_matrix, &R, rhs_data, &R,
			singularValues, &cond, &rank, workspace, &workspace_size, iwork, &info);

	if(info != 0)
		cerr << "Error on linear least squares solver." << endl;

	BLinAlg::Vector* solution = new BLinAlg::Vector(C);
	for(int i = 0; i < C; i++)
		solution->setEntry(i, rhs_data[i]);

	delete [] lapack_matrix;
	delete [] rhs_data;
	delete [] workspace;

	return solution;
}

void DenseMatrix::svd(DenseMatrix *_uMatrix, double *_singularValues, DenseMatrix *_vMatrix)  {
	char returnU = 'A', returnV = 'A';
	double *lapack_matrix = formLAPACKMatrix();
	double *lapack_u = new double[R*R], *lapack_v = new double[C*C];
	double *dummy_workspace = new double[1];
	int workspace_dim = -1, info = 10;

	// format: (A, A, R, C, matrix, R, singular, uMatrix, R, vMatrix, C, workspace, workspace_dim, info)

	// perform svd
	workspace_dim = 7*C;
	double *workspace = new double[workspace_dim];
	FORTRANIZE(dgesvd)(&returnU, &returnV, &R, &C, lapack_matrix, &R, _singularValues, lapack_u, &R, lapack_v, &C, workspace, &workspace_dim, &info);

	// set the matrices - right singular vectors, then left singular vectors
	int ind = 0;
	for(int j = 0; j < R; j++)  {
		for(int i = 0; i < R; i++)  {
			_uMatrix->setEntry(j, i, lapack_u[ind]);
			ind++;
		}
	}

	ind = 0;
	for(int j = 0; j < C; j++)  {
		for(int i = 0; i < C; i++)  {
			_vMatrix->setEntry(j, i, lapack_v[ind]);
			ind++;
		}
	}

	delete [] lapack_matrix;
	delete [] lapack_u;
	delete [] lapack_v;
	delete [] workspace;
	delete [] dummy_workspace;
}

DenseMatrix* DenseMatrix::inverse()  {
	if(R != C)  {
		cerr << "can only invert a square matrix\n";
		return new DenseMatrix(R,R);
	}

	DenseMatrix* uMatrix = new DenseMatrix(R,R);
	DenseMatrix* vMatrixT = new DenseMatrix(R,R);
	double* singularValues = new double[R];

	this->svd(uMatrix, singularValues, vMatrixT);

	DenseMatrix* uMatrixT = uMatrix->transpose();
	DenseMatrix* vMatrix = vMatrixT->transpose();

	DenseMatrix* singularMatrix = new DenseMatrix(R,R);
	double* invSingularValues = new double[R];
	for(int i = 0; i < R; i++)  {
		//cout << "singular value[" << i << "] : " << singularValues[i] << endl;
		invSingularValues[i] = 1.0 / singularValues[i];
	}
	singularMatrix->setDiagonal(invSingularValues);

	DenseMatrix* partialInverted = vMatrix->mult(singularMatrix);
	DenseMatrix* invertedMatrix = partialInverted->mult(uMatrixT);

	delete uMatrix;
	delete vMatrix;
	delete uMatrixT;
	delete vMatrixT;
	delete singularMatrix;
	delete partialInverted;

	delete [] singularValues;
	delete [] invSingularValues;

	return invertedMatrix;
}

void DenseMatrix::printMatrix()  {
	for(int j = 0; j < R; j++)  {
		cout << "| ";
		for(int i = 0; i < C; i++)  {
			float entry = getEntry(i, j);
			cout << entry << " ";
		}
		cout << "|\n";
	}
}
