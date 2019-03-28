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

#ifndef SHAPELOADER_H
#define SHAPELOADER_H

#include "mpu_soup.h"

class ShapeLoader  {
	public:
		static ImplicitFunction* load_shape(string _filename)  {
			string file_type = "";
			bool hit_dot = false;

			for(int l = _filename.length()-1; l >= 0; l--)  {
				string next_char = _filename.substr(l,1);
				if(next_char.compare(".") == 0)  {
					hit_dot = true;
					break;
				}
				file_type = next_char + file_type;
			}

			if(!hit_dot)  {
				cerr << _filename << " has no file extension!" << endl;
				return 0;
			}

			// TODO: add superquadric here!
			if(file_type.compare("mpu") == 0)
				return new MPUSoup(_filename);
			else if(file_type.compare("mpui") == 0)
				return new InvMPUSoup(_filename);
			// else
			cerr << _filename << " is not a valid filename!" << endl;
			return 0;
		}

	private:
};

#endif
