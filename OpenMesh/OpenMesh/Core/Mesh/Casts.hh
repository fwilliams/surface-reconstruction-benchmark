/*===========================================================================*\
 *                                                                           *
 *                               OpenMesh                                    *
 *        Copyright (C) 2004 by Computer Graphics Group, RWTH Aachen         *
 *                           www.openmesh.org                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *                                License                                    *
 *                                                                           *
 *  This library is free software; you can redistribute it and/or modify it  *
 *  under the terms of the GNU Lesser General Public License as published    *
 *  by the Free Software Foundation, version 2.1.                            *
 *                                                                           *
 *  This library is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        *
 *  Lesser General Public License for more details.                          *
 *                                                                           *
 *  You should have received a copy of the GNU Lesser General Public         *
 *  License along with this library; if not, write to the Free Software      *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                *
 *                                                                           *
\*===========================================================================*/

#ifndef OPENMESH_CASTS_HH
#define OPENMESH_CASTS_HH
//== INCLUDES =================================================================

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

//== NAMESPACES ===============================================================
namespace OpenMesh 
{

template <class Traits>
inline TriMesh_ArrayKernelT<Traits>& TRIMESH_CAST(PolyMesh_ArrayKernelT<Traits>& _poly_mesh)
{ return reinterpret_cast< TriMesh_ArrayKernelT<Traits>& >(_poly_mesh); }

template <class Traits>
inline const TriMesh_ArrayKernelT<Traits>& TRIMESH_CAST(const PolyMesh_ArrayKernelT<Traits>& _poly_mesh)
{ return reinterpret_cast< const TriMesh_ArrayKernelT<Traits>& >(_poly_mesh); }

template <class Traits>
inline PolyMesh_ArrayKernelT<Traits>& POLYMESH_CAST(TriMesh_ArrayKernelT<Traits>& _tri_mesh)
{ return reinterpret_cast< PolyMesh_ArrayKernelT<Traits>& >(_tri_mesh); }

template <class Traits>
inline const PolyMesh_ArrayKernelT<Traits>& POLYMESH_CAST(const TriMesh_ArrayKernelT<Traits>& _tri_mesh)
{ return reinterpret_cast< const PolyMesh_ArrayKernelT<Traits>& >(_tri_mesh); }

};
#endif//OPENMESH_CASTS_HH
