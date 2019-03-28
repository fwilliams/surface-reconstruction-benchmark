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

#ifndef VECTOR3_H
#define VECTOR3_H

#include <iostream>
#include <math.h>

class Vector3  {

	public:
		Vector3() : x(0), y(0), z(0) {}
		Vector3(double* _vec) : x(_vec[0]), y(_vec[1]), z(_vec[2]) {}
		Vector3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
		Vector3(const Vector3& _vec) : x(_vec.x), y(_vec.y), z(_vec.z) {}
		Vector3(const Vector3& _to, const Vector3& _from) : x(_to.x-_from.x), y(_to.y-_from.y), z(_to.z-_from.z) {}

		Vector3 operator+(const Vector3& _vec)  {
			Vector3 result(x+_vec.x, y+_vec.y, z+_vec.z);
			return result;
		}

		Vector3 operator-(const Vector3& _vec)  {
			Vector3 result(x-_vec.x, y-_vec.y, z-_vec.z);
			return result;
		}

		double operator*(const Vector3& _vec)  {
			return x*_vec.x+y*_vec.y+z*_vec.z;
		}

		Vector3 operator*(double _scale)  {
			Vector3 result(x*_scale, y*_scale, z*_scale);
			return Vector3(x*_scale, y*_scale, z*_scale);
		}

		inline bool operator==(const Vector3& _otherVec) const  {
			return x == _otherVec.x && y == _otherVec.y && z == _otherVec.z;
		}

		inline bool operator<(const Vector3& _otherVec) const  {
			return x < _otherVec.x && y < _otherVec.y && z < _otherVec.z;
		}

		void scale(double _scale)  {
			x *= _scale;
			y *= _scale;
			z *= _scale;
		}

		double entry(int _i)  {
			if(_i == 0)
				return x;
			else if(_i == 1)
				return y;
			else if(_i == 2)
				return z;
			else
				return x;
		}

		double sqdDist() const {
			return x*x+y*y+z*z;
		}

		double length()  {
			return sqrt(x*x+y*y+z*z);
		}

		double sqdDist(Vector3 _to)  {
			double delX = _to.x-x, delY = _to.y-y, delZ = _to.z-z;
			return delX*delX+delY*delY+delZ*delZ;
		}

		double length(Vector3 _to)  {
			double delX = _to.x-x, delY = _to.y-y, delZ = _to.z-z;
			return sqrt(delX*delX+delY*delY+delZ*delZ);
		}

		void normalize()  {
			double inv_length = 1.0 / length();
			x *= inv_length;
			y *= inv_length;
			z *= inv_length;
		}

		double dotProduct(const Vector3& _vec)  {
			return x*_vec.x+y*_vec.y+z*_vec.z;
		}

		Vector3 cross(const Vector3& _toVec)  {
			double crossX = y*_toVec.z-_toVec.y*z;
			double crossY = _toVec.x*z-x*_toVec.z;
			double crossZ = x*_toVec.y-_toVec.x*y;
			return Vector3(crossX, crossY, crossZ);
		}

		double angle(const Vector3& _toVec)  {
			Vector3 from_vec(x,y,z), to_vec(_toVec);
			from_vec.normalize();
			to_vec.normalize();
			return acos(from_vec.dotProduct(to_vec));
		}

		void expand_min_bound(const Vector3& _vec)  {
			x = _vec.x < x ? _vec.x : x;
			y = _vec.y < y ? _vec.y : y;
			z = _vec.z < z ? _vec.z : z;
		}

		void expand_max_bound(const Vector3& _vec)  {
			x = _vec.x > x ? _vec.x : x;
			y = _vec.y > y ? _vec.y : y;
			z = _vec.z > z ? _vec.z : z;
		}

		friend std::ostream& operator <<(std::ostream &out, const Vector3& _vec)  {
			out << "<" << _vec.x << ", " << _vec.y << ", " << _vec.z << ">";
			return out;
		}

		double x, y, z;
};

#endif
