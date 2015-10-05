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
#ifndef OPENMESH_CIRCULATORS_HH
#define OPENMESH_CIRCULATORS_HH
//=============================================================================
//
//  Vertex and Face circulators for PolyMesh/TriMesh
//
//=============================================================================



//== INCLUDES =================================================================

#include <OpenMesh/Core/System/config.h>
#include <assert.h>


//== NAMESPACES ===============================================================

namespace OpenMesh {
namespace Iterators {


//== FORWARD DECLARATIONS =====================================================


template <class Mesh> class VertexVertexIterT;
template <class Mesh> class VertexIHalfedgeIterT;
template <class Mesh> class VertexOHalfedgeIterT;
template <class Mesh> class VertexEdgeIterT;
template <class Mesh> class VertexFaceIterT;

template <class Mesh> class ConstVertexVertexIterT;
template <class Mesh> class ConstVertexIHalfedgeIterT;
template <class Mesh> class ConstVertexOHalfedgeIterT;
template <class Mesh> class ConstVertexEdgeIterT;
template <class Mesh> class ConstVertexFaceIterT;

template <class Mesh> class FaceVertexIterT;
template <class Mesh> class FaceHalfedgeIterT;
template <class Mesh> class FaceEdgeIterT;
template <class Mesh> class FaceFaceIterT;

template <class Mesh> class ConstFaceVertexIterT;
template <class Mesh> class ConstFaceHalfedgeIterT;
template <class Mesh> class ConstFaceEdgeIterT;
template <class Mesh> class ConstFaceFaceIterT;



