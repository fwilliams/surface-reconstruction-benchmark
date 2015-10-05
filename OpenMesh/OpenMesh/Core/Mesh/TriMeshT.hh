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


//=============================================================================
//
//  CLASS TriMeshT
//
//=============================================================================


#ifndef OPENMESH_TRIMESH_HH
#define OPENMESH_TRIMESH_HH


//== INCLUDES =================================================================


#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/PolyMeshT.hh>
#include <vector>


//== NAMESPACES ===============================================================


namespace OpenMesh {


//== CLASS DEFINITION =========================================================


/** \class TriMeshT TriMeshT.hh <OpenMesh/Mesh/TriMeshT.hh>

    Base type for a triangle mesh.

    Base type for a triangle mesh, parameterized by a mesh kernel. The
    mesh inherits all methods from the kernel class and the
    more general polygonal mesh PolyMeshT. Therefore it provides
    the same types for items, handles, iterators and so on.

    \param Kernel: template argument for the mesh kernel
    \note You should use the predefined mesh-kernel combinations in
    \ref mesh_types_group
    \see \ref mesh_type
    \see OpenMesh::PolyMeshT
*/

template <class Kernel>
class TriMeshT : public PolyMeshT<Kernel>
{

public:


  // self
  typedef TriMeshT<Kernel>                      This;
  typedef PolyMeshT<Kernel>                     PolyMesh;

  //@{
  /// Determine whether this is a PolyMeshT or TriMeshT
  enum { IsPolyMesh = 0 };
  enum { IsTriMesh  = 1 };
  static bool is_polymesh() { return false; }
  static bool is_trimesh()  { return  true; }
  //@}

  //--- items ---

  typedef typename PolyMesh::Scalar             Scalar;
  typedef typename PolyMesh::Point              Point;
  typedef typename PolyMesh::Normal             Normal;
  typedef typename PolyMesh::Color              Color;
  typedef typename PolyMesh::TexCoord1D         TexCoord1D;
  typedef typename PolyMesh::TexCoord2D         TexCoord2D;
  typedef typename PolyMesh::TexCoord3D         TexCoord3D;
  typedef typename PolyMesh::Vertex             Vertex;
  typedef typename PolyMesh::Halfedge           Halfedge;
  typedef typename PolyMesh::Edge               Edge;
  typedef typename PolyMesh::Face               Face;


  //--- handles ---

  typedef typename PolyMesh::VertexHandle       VertexHandle;
  typedef typename PolyMesh::HalfedgeHandle     HalfedgeHandle;
  typedef typename PolyMesh::EdgeHandle         EdgeHandle;
  typedef typename PolyMesh::FaceHandle         FaceHandle;


  //--- iterators ---

  typedef typename PolyMesh::VertexIter         VertexIter;
  typedef typename PolyMesh::ConstVertexIter    ConstVertexIter;
  typedef typename PolyMesh::EdgeIter           EdgeIter;
  typedef typename PolyMesh::ConstEdgeIter      ConstEdgeIter;
  typedef typename PolyMesh::FaceIter           FaceIter;
  typedef typename PolyMesh::ConstFaceIter      ConstFaceIter;



  //--- circulators ---

  typedef typename PolyMesh::VertexVertexIter         VertexVertexIter;
  typedef typename PolyMesh::VertexOHalfedgeIter      VertexOHalfedgeIter;
  typedef typename PolyMesh::VertexIHalfedgeIter      VertexIHalfedgeIter;
  typedef typename PolyMesh::VertexEdgeIter           VertexEdgeIter;
  typedef typename PolyMesh::VertexFaceIter           VertexFaceIter;
  typedef typename PolyMesh::FaceVertexIter           FaceVertexIter;
  typedef typename PolyMesh::FaceHalfedgeIter         FaceHalfedgeIter;
  typedef typename PolyMesh::FaceEdgeIter             FaceEdgeIter;
  typedef typename PolyMesh::FaceFaceIter             FaceFaceIter;
  typedef typename PolyMesh::ConstVertexVertexIter    ConstVertexVertexIter;
  typedef typename PolyMesh::ConstVertexOHalfedgeIter ConstVertexOHalfedgeIter;
  typedef typename PolyMesh::ConstVertexIHalfedgeIter ConstVertexIHalfedgeIter;
  typedef typename PolyMesh::ConstVertexEdgeIter      ConstVertexEdgeIter;
  typedef typename PolyMesh::ConstVertexFaceIter      ConstVertexFaceIter;
  typedef typename PolyMesh::ConstFaceVertexIter      ConstFaceVertexIter;
  typedef typename PolyMesh::ConstFaceHalfedgeIter    ConstFaceHalfedgeIter;
  typedef typename PolyMesh::ConstFaceEdgeIter        ConstFaceEdgeIter;
  typedef typename PolyMesh::ConstFaceFaceIter        ConstFaceFaceIter;

  // --- constructor/destructor

  /// Default constructor
  TriMeshT() : PolyMesh() {}
  /// Destructor
  virtual ~TriMeshT() {}

  //--- halfedge collapse / vertex split ---

  /// Vertex Split: inverse operation to collapse().
  inline HalfedgeHandle vertex_split(Point _v0_point,  VertexHandle _v1,
                                     VertexHandle _vl, VertexHandle _vr)
  { return PolyMesh::vertex_split(add_vertex(_v0_point), _v1, _vl, _vr); }

  inline HalfedgeHandle vertex_split(VertexHandle _v0, VertexHandle _v1,
                                     VertexHandle _vl, VertexHandle _vr)
  { return PolyMesh::vertex_split(_v0, _v1, _vl, _vr); }

  /// Edge split (= 2-to-4 split)
  inline void split(EdgeHandle _eh, const Point& _p)
  { PolyMesh::split(_eh, add_vertex(_p)); }

  inline void split(EdgeHandle _eh, VertexHandle _vh)
  { PolyMesh::split(_eh, _vh); }

  /// Face split (= 1-to-3 split, calls corresponding PolyMeshT function).

  /// Face split (= 1-to-3 split, calls corresponding PolyMeshT function).
  inline void split(FaceHandle _fh, const Point& _p)
  { PolyMesh::split(_fh, _p); }

  inline void split(FaceHandle _fh, VertexHandle _vh)
  { PolyMesh::split(_fh, _vh); }
};


//=============================================================================
} // namespace OpenMesh
//=============================================================================
#if defined(OM_INCLUDE_TEMPLATES) && !defined(OPENMESH_TRIMESH_C)
#define OPENMESH_TRIMESH_TEMPLATES
#include "TriMeshT.cc"
#endif
//=============================================================================
#endif // OPENMESH_TRIMESH_HH defined
//=============================================================================
