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

#ifndef RIGID_H
#define RIGID_H

#include "DenseEigenSystem.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

class Rigid  {
	public:
		Rigid();
		Rigid(string _filename);
		Rigid(vector<Vector3> _srcPoints, vector<Vector3> _dstPoints);
		~Rigid();

		void update_rigid(vector<Vector3> _srcPoints, vector<Vector3> _dstPoints);
		Vector3 apply_transformation(Vector3 _pt);
		Vector3 apply_rotation(Vector3 _vec);

		void write_to_file(string _filename)  {
			FILE* f_out = fopen(_filename.c_str(), "wb");
			double flattened_matrix[] = { rotation->getEntry(0,0), rotation->getEntry(1,0), rotation->getEntry(2,0), 
													rotation->getEntry(0,1), rotation->getEntry(1,1), rotation->getEntry(2,1), 
													rotation->getEntry(0,2), rotation->getEntry(1,2), rotation->getEntry(2,2) };
			double flattened_tran[] = { translation.x, translation.y, translation.z };
			fwrite(flattened_matrix, sizeof(double), 9, f_out);
			fwrite(flattened_tran, sizeof(double), 3, f_out);
			fclose(f_out);
		}

	private:
		DenseMatrix* rotation;
		Vector3 translation;

		void solve_rigid(vector<Vector3> _srcPoints, vector<Vector3> _dstPoints, DenseMatrix* rot, Vector3& tran);
};

#endif
