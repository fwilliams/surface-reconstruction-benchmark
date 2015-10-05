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

#include "DenseEigenSystem.h"

extern "C"  {
	void FORTRANIZE(dspevx)(
			char* all, char* range, char* upper, int* dim, double* matrix, double* dummy, double* dummy1, int* dummy2, int* dummy3,
			double* tolerance, int* numEigen, double* eigvalues, double* eigvectors, int* mat_dim, double* workspace, int* workspace1,
			int* fail, int* info
	);

	double FORTRANIZE(dlamch)(
			char* prec
	);
}

DenseEigenSystem::DenseEigenSystem(DenseMatrix *_matrix, int _il, int _iu)  {
	matrix = _matrix;
	il = _il+1;
	iu = _iu+1;
}

DenseEigenSystem::~DenseEigenSystem()  {
	delete [] eigenvalues;
	for(int i = 0; i < systemSize; i++)
		delete eigenvectors[i];
	delete [] eigenvectors;
}

void DenseEigenSystem::solve()  {
	// dspevx(V, A, U, N, matrix, vl, vu, il, iu, tolerance, numEigen, values, vectors, N, work_space, iwork_space, fail, info)

	// dlamch('S') - accuracy
	char epsilon = 'S';
	double tolerance = 2.0*FORTRANIZE(dlamch)(&epsilon);

	// setup upper-triangular matrix from symmetric matrix
	int mat_dim = matrix->rows();
	int triangularStorage = (mat_dim * (mat_dim+1)) / 2;
	double *upperTriangular = new double[triangularStorage];
	int ind = 0;
	for(int i = 0; i < mat_dim; i++)  {
		for(int j = 0; j <= i; j++)  {
			upperTriangular[ind] = matrix->getEntry(i, j);
			ind++;
		}
	}

	// setup for the call!

	// inputs first ...
	int dummy_int0 = 0, dummy_int1 = 0;
	double dummy_double0 = 0, dummy_double1 = 0;
	char all = 'V', range = 'I', upper = 'U';

	// ... now outputs
	int numEigen = 0;
	double *flattened_vectors = new double[(iu-il+1)*mat_dim];
	double *flattened_values = new double[mat_dim];
	double *work_space = new double[10*mat_dim];
	int *iwork_space = new int[6*mat_dim], info, *fail = new int[mat_dim];

	// compute eigenvectors & eigenvalues
	FORTRANIZE(dspevx)(&all, &range, &upper, &mat_dim, upperTriangular, &dummy_double0, &dummy_double1, &il, &iu,
			 &tolerance, &numEigen, flattened_values, flattened_vectors, &mat_dim, work_space, iwork_space, fail, &info);

	if(info != 0)  {
		cerr << "error! something went wrong in the eigensystem decomposition...\n";
		return;
	}
	if(numEigen != (iu-il+1))
		cerr << "warning! degenerate eigensystem, but we are going to plow on through ...\n";

	systemSize = numEigen;
	eigenvectors = new BLinAlg::Vector*[systemSize];
	eigenvalues = new double[systemSize];
	for(int i = 0; i < systemSize; i++)  {
		double eigenColumn[mat_dim];
		int start_dim = mat_dim*i;
		for(int j = 0; j < mat_dim; j++)
			eigenColumn[j] = flattened_vectors[start_dim+j];
		eigenvectors[i] = new BLinAlg::Vector(eigenColumn, mat_dim);
		eigenvalues[i] = flattened_values[i];
	}

	if(CHECK_EIGENSYSTEM)  {
		for(int i = 0; i < systemSize; i++)  {
			BLinAlg::Vector *vector1 = new BLinAlg::Vector(eigenvectors[i]);
			BLinAlg::Vector *vector2 = matrix->mult(vector1);
			vector1->scale(eigenvalues[i]);

			cerr << "Eigenvalue[" << i << "] : " << eigenvalues[i] << endl;
			cerr << "Vector 1[" << i << "] : " << vector1 << endl;
			cerr << "Vector 2[" << i << "] : " << vector2 << endl;

			delete vector1;
			delete vector2;
		}
	}

	delete[] upperTriangular;
	delete[] flattened_values;
	delete[] flattened_vectors;
	delete[] work_space;
	delete[] iwork_space;
	delete[] fail;
}
