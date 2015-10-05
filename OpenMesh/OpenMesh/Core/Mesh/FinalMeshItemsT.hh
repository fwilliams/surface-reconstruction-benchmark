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

#ifndef OPENMESH_MESH_ITEMS_HH
#define OPENMESH_MESH_ITEMS_HH


//== INCLUDES =================================================================


#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Utils/GenProg.hh>
#include <OpenMesh/Core/Utils/vector_traits.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>


//== NAMESPACES ===============================================================


namespace OpenMesh {


//== CLASS DEFINITION =========================================================

/// Definition of the mesh entities (items).
template <class Traits, bool IsTriMesh>
struct FinalMeshItemsT
{
  //--- build Refs structure ---
#ifndef DOXY_IGNORE_THIS
  struct Refs
  {
    typedef typename Traits::Point            Point;
    typedef typename vector_traits<Point>::value_type Scalar;

    typedef typename Traits::Normal           Normal;
    typedef typename Traits::Color            Color;
    typedef typename Traits::TexCoord1D       TexCoord1D;
    typedef typename Traits::TexCoord2D       TexCoord2D;
    typedef typename Traits::TexCoord3D       TexCoord3D;
    typedef OpenMesh::VertexHandle            VertexHandle;
    typedef OpenMesh::FaceHandle              FaceHandle;
    typedef OpenMesh::EdgeHandle              EdgeHandle;
    typedef OpenMesh::HalfedgeHandle          HalfedgeHandle;
  };
#endif
  //--- export Refs types ---
  typedef typename Refs::Point           Point;
  typedef typename Refs::Scalar          Scalar;
  typedef typename Refs::Normal          Normal;
  typedef typename Refs::Color           Color;
  typedef typename Refs::TexCoord1D      TexCoord1D;
  typedef typename Refs::TexCoord2D      TexCoord2D;
  typedef typename Refs::TexCoord3D      TexCoord3D;

  //--- get attribute bits from Traits ---
  enum Attribs
  {
    VAttribs = Traits::VertexAttributes,
    HAttribs = Traits::HalfedgeAttributes,
    EAttribs = Traits::EdgeAttributes,
    FAttribs = Traits::FaceAttributes
  };
  //--- merge internal items with traits items ---


/*
  typedef typename GenProg::IF<
    (bool)(HAttribs & Attributes::PrevHalfedge),
    typename InternalItems::Halfedge_with_prev,
    typename InternalItems::Halfedge_without_prev
  >::Result   InternalHalfedge;
*/
  //typedef typename InternalItems::Vertex                     InternalVertex;
  //typedef typename InternalItems::template Edge<Halfedge>      InternalEdge;
  //typedef typename InternalItems::template Face<IsTriMesh>     InternalFace;
  class ITraits
  {};

  typedef typename Traits::template VertexT<ITraits, Refs>      VertexData;
  typedef typename Traits::template HalfedgeT<ITraits, Refs>    HalfedgeData;
  typedef typename Traits::template EdgeT<ITraits, Refs>        EdgeData;
  typedef typename Traits::template FaceT<ITraits, Refs>        FaceData;
};


//=============================================================================
} // namespace OpenMesh
//=============================================================================
#endif // OPENMESH_MESH_ITEMS_HH defined
//=============================================================================

