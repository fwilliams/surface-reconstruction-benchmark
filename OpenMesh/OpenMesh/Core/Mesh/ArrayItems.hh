//=============================================================================
//                                                                            
//                               OpenMesh                                     
//        Copyright (C) 2003 by Computer Graphics Group, RWTH Aachen          
//                           www.openmesh.org                                 
//                                                                            
//-----------------------------------------------------------------------------
//                                                                            
//                                License                                     
//                                                                            
//   This library is free software; you can redistribute it and/or modify it 
//   under the terms of the GNU Lesser General Public License as published   
//   by the Free Software Foundation, version 2.1.                           
//                                                                             
//   This library is distributed in the hope that it will be useful, but       
//   WITHOUT ANY WARRANTY; without even the implied warranty of                
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU         
//   Lesser General Public License for more details.                           
//                                                                            
//   You should have received a copy of the GNU Lesser General Public          
//   License along with this library; if not, write to the Free Software       
//   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 
//
//-----------------------------------------------------------------------------
//
//   $Revision: 1801 $
//   $Date: 2008-05-19 11:53:56 +0200 (Mon, 19 May 2008) $
//
//=============================================================================

#ifndef OPENMESH_ARRAY_ITEMS_HH
#define OPENMESH_ARRAY_ITEMS_HH


//== INCLUDES =================================================================


#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Utils/GenProg.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>


//== NAMESPACES ===============================================================

namespace OpenMesh {


//== CLASS DEFINITION =========================================================


/// Definition of mesh items for use in the ArrayKernel
struct ArrayItems
{

  //------------------------------------------------------ internal vertex type

  /// The vertex item
  class Vertex
  {
    friend class ArrayKernel;
    HalfedgeHandle  halfedge_handle_;
  };


  //---------------------------------------------------- internal halfedge type

#ifndef DOXY_IGNORE_THIS
  class Halfedge_without_prev
  {
    friend class ArrayKernel;
    FaceHandle      face_handle_;
    VertexHandle    vertex_handle_;
    HalfedgeHandle  next_halfedge_handle_;
  };
#endif

#ifndef DOXY_IGNORE_THIS
  class Halfedge_with_prev : public Halfedge_without_prev
  {
    friend class ArrayKernel;
    HalfedgeHandle  prev_halfedge_handle_;
  };
#endif

  //TODO: should be selected with config.h define
  typedef Halfedge_with_prev                Halfedge;
  typedef GenProg::Bool2Type<true>          HasPrevHalfedge;

  //-------------------------------------------------------- internal edge type
#ifndef DOXY_IGNORE_THIS
  class Edge
  {
    friend class ArrayKernel;
    Halfedge  halfedges_[2];
  };
#endif

  //-------------------------------------------------------- internal face type
#ifndef DOXY_IGNORE_THIS
  class Face
  {
    friend class ArrayKernel;
    HalfedgeHandle  halfedge_handle_;
  };
};
#endif

//=============================================================================
} // namespace OpenMesh
//=============================================================================
#endif // OPENMESH_ITEMS_HH defined
//=============================================================================
