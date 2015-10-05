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

#include "Vector.h"

namespace BLinAlg  {

Vector::Vector()  {
	n = 0;
}

Vector::Vector(int _n)  {
	n = _n;
	d = new double[_n];
	for(int i = 0; i < n; i++)
		d[i] = 0.0;
}

Vector::Vector(double *_d, int _n)  {
	n = _n;
	d = new double[_n];
	for(int i = 0; i < n; i++)
		d[i] = _d[i];
}

Vector::Vector(Vector* _vector)  {
	n = _vector->dim();
	d = new double[n];
	for(int i = 0; i < n; i++)
		d[i] = _vector->getEntry(i);
}

Vector::~Vector()  {
	delete [] d;
}

double Vector::l2Norm()  {
	double ssd = 0.0;
	for(int i = 0; i < n; i++)
		ssd += d[i]*d[i];
	return sqrt(ssd);
}

void Vector::normalize()  {
	double inv_length = 1.0 / l2Norm();
	for(int i = 0; i < n; i++)
		d[i] *= inv_length;
}

double Vector::dot(Vector* _vec)  {
	if(_vec->dim() != n)  {
		cerr << "[VECTOR] unable to perform dot product: mismatched dimensionality." << endl;
		return 0.0;
	}

	double dotted = 0.0;
	for(int i = 0; i < n; i++)
		dotted += d[i]*_vec->getEntry(i);
	return dotted;
}

void Vector::add(Vector* _vec)  {
	if(_vec->dim() != n)  {
		cerr << "[VECTOR] unable to perform addition: mismatched dimensionality." << endl;
		return;
	}

	for(int i = 0; i < n; i++)
		d[i] += _vec->getEntry(i);
}

void Vector::sub(Vector* _vec)  {
	if(_vec->dim() != n)  {
		cerr << "[VECTOR] unable to perform subtraction: mismatched dimensionality." << endl;
		return;
	}

	for(int i = 0; i < n; i++)
		d[i] -= _vec->getEntry(i);
}

void Vector::scale(double _scale)  {
	for(int i = 0; i < n; i++)
		d[i] *= _scale;
}

void Vector::minAndMax(int* min_index, double* minimum, int* max_index, double* maximum)  {
	double the_min = 1e10, the_max = -1e10;
	int min_ind = -1, max_ind = -1;
	for(int i = 0; i < n; i++)  {
		if(d[i] < the_min)  {
			the_min = d[i];
			min_ind = i;
		}
		if(d[i] > the_max)  {
			the_max = d[i];
			max_ind = i;
		}
	}

	*min_index = min_ind;
	*minimum = the_min;
	*max_index = max_ind;
	*maximum = the_max;
}

void Vector::min(int* index, double* minimum)  {
	double the_min = 1e10;
	int the_index = -1;
	for(int i = 0; i < n; i++)  {
		if(d[i] < the_min)  {
			the_min = d[i];
			the_index = i;
		}
	}
	*index = the_index;
	*minimum = the_min;
}

void Vector::max(int* index, double* maximum)  {
	double the_max = -1e10;
	int the_index = -1;
	for(int i = 0; i < n; i++)  {
		if(d[i] > the_max)  {
			the_max = d[i];
			the_index = i;
		}
	}
	*index = the_index;
	*maximum = the_max;
}

}
