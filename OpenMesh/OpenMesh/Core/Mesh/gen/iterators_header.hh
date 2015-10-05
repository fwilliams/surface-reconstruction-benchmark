/*===========================================================================*\
 *                                                                           *
 *                               OpenMesh                                    *
 *      Copyright (C) 2001-2003 by Computer Graphics Group, RWTH Aachen      *
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
#ifndef OPENMESH_ITERATORS_HH
#define OPENMESH_ITERATORS_HH
//=============================================================================
//
//  Iterators for PolyMesh/TriMesh
//
//=============================================================================



//== INCLUDES =================================================================

#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Attributes/Status.hh>
#include <assert.h>


//== NAMESPACES ===============================================================

namespace OpenMesh {
namespace Iterators {


//== FORWARD DECLARATIONS =====================================================


template <class Mesh> class VertexIterT;
template <class Mesh> class ConstVertexIterT;
template <class Mesh> class HalfedgeIterT;
template <class Mesh> class ConstHalfedgeIterT;
template <class Mesh> class EdgeIterT;
template <class Mesh> class ConstEdgeIterT;
template <class Mesh> class FaceIterT;
template <class Mesh> class ConstFaceIterT;




