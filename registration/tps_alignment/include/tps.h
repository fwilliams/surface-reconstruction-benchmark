/*
Copyright 2007 Benedict Brown
Princeton University

tps.h
Prototypes, etc. for tps alignment

-----

This file is part of tps_alignment.

tps_alignment is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.

tps_alignment is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef TPS_H_INCLUDED
#define TPS_H_INCLUDED

#include "tnt_array2d.h"
typedef double matrix_t;
typedef TNT::Array2D<matrix_t> farr;


void update_transform_mat(const farr &x, const farr &y, matrix_t lambda,
                          farr &w, farr &A);
farr warp_points(const farr &pts, const farr &x, const farr &w, const farr &A);

#endif /* TPS_H_INCLUDED */
