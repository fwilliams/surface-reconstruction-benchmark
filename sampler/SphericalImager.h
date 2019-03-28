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

#ifndef SPHERICALIMAGER_H
#define SPHERICALIMAGER_H

#include "../evaluator/ComputationTimer.h"

#include "RayCaster.h"
#include "OrientedRangeImage.h"
#include "RangeImage.h"

#include "NPTSReader.h"
#include "PointCloud.h"

#include "RangeScanner.h"
#include "../modeling/ImplicitSphere.h"
#include "../modeling/RotationMatrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

class SphericalImager  {
	public:
		SphericalImager(ImplicitFunction* _implicit, int _numScans, int _resX, int _resY);
		~SphericalImager();

		// all angle parameters in degrees

		void setMinRange(double _minRange)  { min_range = _minRange; }

		void setRandomSampleRotation(bool _toRotate)  {
			to_random_sample_rotation = _toRotate;
		}

		void sample(string _filebase);

	private:
		ImplicitFunction* implicit_function;

		double my_pi;

		// required
		int num_scans;
		int res_x, res_y;

		// optional

		// range scanner parameters
		double min_range;

		// misc
		bool to_random_sample_rotation;
};

#endif
