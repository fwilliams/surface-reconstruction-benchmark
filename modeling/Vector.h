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

#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <math.h>

using namespace std;

namespace BLinAlg  {

class Vector  {
	public:
		Vector();
		Vector(int _n);
		Vector(double *_d, int _n);
		Vector(Vector* _vector);

		~Vector();

		int dim()  { return n; };
		double* data()  { return d; };
		void setEntry(int _n, double _entry)  { d[_n] = _entry; }
		void incrementEntry(int _n, double _entry)  { d[_n] += _entry; }
		void scaleEntry(int _n, double _scale)  { d[_n] *= _scale; }
		double getEntry(int _n)  { return d[_n]; }

		double* dataCopy()  {
			double* copy = new double[n];
			for(int i = 0; i < n; i++)
				copy[i] = d[i];
			return copy;
		}

		void add(Vector* _vec);
		void sub(Vector* _vec);
		void scale(double _scale);
		double dot(Vector* _vec);
		double l2Norm();
		void normalize();

		void minAndMax(int* min_index, double* minimum, int* max_index, double* maximum);
		void min(int* index, double* minimum);
		void max(int* index, double* maximum);

	private:
		double* d;
		int n;

};

}

#endif
