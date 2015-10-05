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

#ifndef OPENMESH_POLYCONNECTIVITY_HH
#define OPENMESH_POLYCONNECTIVITY_HH

#include <OpenMesh/Core/Mesh/ArrayKernel.hh>
#include <OpenMesh/Core/Mesh/IteratorsT.hh>
#include <OpenMesh/Core/Mesh/CirculatorsT.hh>

namespace OpenMesh
{

class PolyConnectivity : public ArrayKernel
{
public:
  /// \name Mesh Handles
  //@{
  /// Invalid handle
  static const VertexHandle                           InvalidVertexHandle;
  /// Invalid handle
  static const HalfedgeHandle                         InvalidHalfedgeHandle;
  /// Invalid handle
  static const EdgeHandle                             InvalidEdgeHandle;
  /// Invalid handle
  static const FaceHandle                             InvalidFaceHandle;
  //@}

  typedef PolyConnectivity                            This;

  //--- iterators ---

  /** \name Mesh Iterators
      Refer to OpenMesh::Mesh::Iterators or \ref mesh_iterators for
      documentation.
  */
  //@{
  /// Linear iterator
  typedef Iterators::VertexIterT<This>                VertexIter;
  typedef Iterators::HalfedgeIterT<This>              HalfedgeIter;
  typedef Iterators::EdgeIterT<This>                  EdgeIter;
  typedef Iterators::FaceIterT<This>                  FaceIter;

  typedef Iterators::ConstVertexIterT<This>           ConstVertexIter;
  typedef Iterators::ConstHalfedgeIterT<This>         ConstHalfedgeIter;
  typedef Iterators::ConstEdgeIterT<This>             ConstEdgeIter;
  typedef Iterators::ConstFaceIterT<This>             ConstFaceIter;
  //@}

  //--- circulators ---

  /** \name Mesh Circulators
      Refer to OpenMesh::Mesh::Iterators or \ref mesh_iterators
      for documentation.
  */
  //@{
  /// Circulator
  typedef Iterators::VertexVertexIterT<This>          VertexVertexIter;
  typedef Iterators::VertexOHalfedgeIterT<This>       VertexOHalfedgeIter;
  typedef Iterators::VertexIHalfedgeIterT<This>       VertexIHalfedgeIter;
  typedef Iterators::VertexEdgeIterT<This>            VertexEdgeIter;
  typedef Iterators::VertexFaceIterT<This>            VertexFaceIter;
  typedef Iterators::FaceVertexIterT<This>            FaceVertexIter;
  typedef Iterators::FaceHalfedgeIterT<This>          FaceHalfedgeIter;
  typedef Iterators::FaceEdgeIterT<This>              FaceEdgeIter;
  typedef Iterators::FaceFaceIterT<This>              FaceFaceIter;

  typedef Iterators::ConstVertexVertexIterT<This>     ConstVertexVertexIter;
  typedef Iterators::ConstVertexOHalfedgeIterT<This>  ConstVertexOHalfedgeIter;
  typedef Iterators::ConstVertexIHalfedgeIterT<This>  ConstVertexIHalfedgeIter;
  typedef Iterators::ConstVertexEdgeIterT<This>       ConstVertexEdgeIter;
  typedef Iterators::ConstVertexFaceIterT<This>       ConstVertexFaceIter;
  typedef Iterators::ConstFaceVertexIterT<This>       ConstFaceVertexIter;
  typedef Iterators::ConstFaceHalfedgeIterT<This>     ConstFaceHalfedgeIter;
  typedef Iterators::ConstFaceEdgeIterT<This>         ConstFaceEdgeIter;
  typedef Iterators::ConstFaceFaceIterT<This>         ConstFaceFaceIter;
  //@}

  // --- shortcuts

  /** \name Typedef Shortcuts
      Provided for convenience only
  */
  //@{
  /// Alias typedef
  typedef VertexHandle    VHandle;
  typedef HalfedgeHandle  HHandle;
  typedef EdgeHandle      EHandle;
  typedef FaceHandle      FHandle;

  typedef VertexIter    VIter;
  typedef HalfedgeIter  HIter;
  typedef EdgeIter      EIter;
  typedef FaceIter      FIter;

  typedef ConstVertexIter    CVIter;
  typedef ConstHalfedgeIter  CHIter;
  typedef ConstEdgeIter      CEIter;
  typedef ConstFaceIter      CFIter;

  typedef VertexVertexIter      VVIter;
  typedef VertexOHalfedgeIter   VOHIter;
  typedef VertexIHalfedgeIter   VIHIter;
  typedef VertexEdgeIter        VEIter;
  typedef VertexFaceIter        VFIter;
  typedef FaceVertexIter        FVIter;
  typedef FaceHalfedgeIter      FHIter;
  typedef FaceEdgeIter          FEIter;
  typedef FaceFaceIter          FFIter;

  typedef ConstVertexVertexIter      CVVIter;
  typedef ConstVertexOHalfedgeIter   CVOHIter;
  typedef ConstVertexIHalfedgeIter   CVIHIter;
  typedef ConstVertexEdgeIter        CVEIter;
  typedef ConstVertexFaceIter        CVFIter;
  typedef ConstFaceVertexIter        CFVIter;
  typedef ConstFaceHalfedgeIter      CFHIter;
  typedef ConstFaceEdgeIter          CFEIter;
  typedef ConstFaceFaceIter          CFFIter;
  //@}

public:

  PolyConnectivity()  {}
  virtual ~PolyConnectivity() {}

  inline static bool is_triangles()
  { return false; }

  /** assign_connectivity() method. See ArrayKernel::assign_connectivity()
      for more details. */
  inline void assign_connectivity(const PolyConnectivity& _other)
  { ArrayKernel::assign_connectivity(_other); }
  
  /** \name Adding items to a mesh
  */
  //@{
  /// Add a new vertex 
  inline VertexHandle add_vertex()
  { return new_vertex(); }

  /// Add and connect a new face
  FaceHandle add_face(const std::vector<VertexHandle>& _vhandles)
  { return add_face(&_vhandles.front(), _vhandles.size()); }
  
  FaceHandle add_face(VertexHandle _vh0, VertexHandle _vh1, VertexHandle _vh2)
  { 
    VertexHandle vhs[3] = { _vh0, _vh1, _vh2 };
    return add_face(vhs, 3); 
  }
  
  FaceHandle add_face(VertexHandle _vh0, VertexHandle _vh1, VertexHandle _vh2, VertexHandle _vh3)
  { 
    VertexHandle vhs[4] = { _vh0, _vh1, _vh2, _vh3 };
    return add_face(vhs, 4); 
  }
  
  FaceHandle add_face(const VertexHandle* _vhandles, uint _vhs_size);
  //@}

  /// \name Deleting mesh items and other connectivity/topology modifications
  //@{

  /** Mark vertex and all incident edges and faces deleted.
      Items marked deleted will be removed by garbageCollection().
      \attention Needs the Attributes::Status attribute for vertices,
      edges and faces.
  */
  void delete_vertex(VertexHandle _vh, bool _delete_isolated_vertices = true);

  /** Mark edge (two opposite halfedges) and incident faces deleted.
      Resulting isolated vertices are marked deleted if
      _delete_isolated_vertices is true. Items marked deleted will be
      removed by garbageCollection().

      \attention Needs the Attributes::Status attribute for vertices,
      edges and faces.
  */
  void delete_edge(EdgeHandle _eh, bool _delete_isolated_vertices=true);

  /** Delete face _fh and resulting degenerated empty halfedges as
      well.  Resultling isolated vertices will be deleted if
      _delete_isolated_vertices is true.

      \attention All item will only be marked to be deleted. They will
      actually be removed by calling garbage_collection().

      \attention Needs the Attributes::Status attribute for vertices,
      edges and faces.
  */
  void delete_face(FaceHandle _fh, bool _delete_isolated_vertices=true);

  //@}
  
  /** \name Begin and end iterators
  */
  //@{

  /// Begin iterator for vertices
  VertexIter vertices_begin()
  { return VertexIter(*this, VertexHandle(0)); }
  /// Const begin iterator for vertices
  ConstVertexIter vertices_begin() const
  { return ConstVertexIter(*this, VertexHandle(0)); }
  /// End iterator for vertices
  VertexIter vertices_end()
  { return VertexIter(*this, VertexHandle(n_vertices())); }
  /// Const end iterator for vertices
  ConstVertexIter vertices_end() const
  { return ConstVertexIter(*this, VertexHandle(n_vertices())); }

  /// Begin iterator for halfedges
  HalfedgeIter halfedges_begin()
  { return HalfedgeIter(*this, HalfedgeHandle(0)); }
  /// Const begin iterator for halfedges
  ConstHalfedgeIter halfedges_begin() const
  { return ConstHalfedgeIter(*this, HalfedgeHandle(0)); }
  /// End iterator for halfedges
  HalfedgeIter halfedges_end()
  { return HalfedgeIter(*this, HalfedgeHandle(n_halfedges())); }
  /// Const end iterator for halfedges
  ConstHalfedgeIter halfedges_end() const
  { return ConstHalfedgeIter(*this, HalfedgeHandle(n_halfedges())); }

  /// Begin iterator for edges
  EdgeIter edges_begin()
  { return EdgeIter(*this, EdgeHandle(0)); }
  /// Const begin iterator for edges
  ConstEdgeIter edges_begin() const
  { return ConstEdgeIter(*this, EdgeHandle(0)); }
  /// End iterator for edges
  EdgeIter edges_end()
  { return EdgeIter(*this, EdgeHandle(n_edges())); }
  /// Const end iterator for edges
  ConstEdgeIter edges_end() const
  { return ConstEdgeIter(*this, EdgeHandle(n_edges())); }

  /// Begin iterator for faces
  FaceIter faces_begin()
  { return FaceIter(*this, FaceHandle(0)); }
  /// Const begin iterator for faces
  ConstFaceIter faces_begin() const
  { return ConstFaceIter(*this, FaceHandle(0)); }
  /// End iterator for faces
  FaceIter faces_end()
  { return FaceIter(*this, FaceHandle(n_faces())); }
  /// Const end iterator for faces
  ConstFaceIter faces_end() const
  { return ConstFaceIter(*this, FaceHandle(n_faces())); }

  //@}



  /** \name Begin for skipping iterators
  */
  //@{

  /// Begin iterator for vertices
  VertexIter vertices_sbegin()
  { return VertexIter(*this, VertexHandle(0), true); }
  /// Const begin iterator for vertices
  ConstVertexIter vertices_sbegin() const
  { return ConstVertexIter(*this, VertexHandle(0), true); }

  /// Begin iterator for halfedges
  HalfedgeIter halfedges_sbegin()
  { return HalfedgeIter(*this, HalfedgeHandle(0), true); }
  /// Const begin iterator for halfedges
  ConstHalfedgeIter halfedges_sbegin() const
  { return ConstHalfedgeIter(*this, HalfedgeHandle(0), true); }

  /// Begin iterator for edges
  EdgeIter edges_sbegin()
  { return EdgeIter(*this, EdgeHandle(0), true); }
  /// Const begin iterator for edges
  ConstEdgeIter edges_sbegin() const
  { return ConstEdgeIter(*this, EdgeHandle(0), true); }

  /// Begin iterator for faces
  FaceIter faces_sbegin()
  { return FaceIter(*this, FaceHandle(0), true); }
  /// Const begin iterator for faces
  ConstFaceIter faces_sbegin() const
  { return ConstFaceIter(*this, FaceHandle(0), true); }

  //@}

  //--- circulators ---

  /** \name Vertex and Face circulators
  */
  //@{

  /// vertex - vertex circulator
  VertexVertexIter vv_iter(VertexHandle _vh) {
    return VertexVertexIter(*this, _vh); }
  /// vertex - incoming halfedge circulator
  VertexIHalfedgeIter vih_iter(VertexHandle _vh)
  { return VertexIHalfedgeIter(*this, _vh); }
  /// vertex - outgoing halfedge circulator
  VertexOHalfedgeIter voh_iter(VertexHandle _vh)
  { return VertexOHalfedgeIter(*this, _vh); }
  /// vertex - edge circulator
  VertexEdgeIter ve_iter(VertexHandle _vh)
  { return VertexEdgeIter(*this, _vh); }
  /// vertex - face circulator
  VertexFaceIter vf_iter(VertexHandle _vh)
  { return VertexFaceIter(*this, _vh); }

  /// const vertex circulator
  ConstVertexVertexIter cvv_iter(VertexHandle _vh) const
  { return ConstVertexVertexIter(*this, _vh); }
  /// const vertex - incoming halfedge circulator
  ConstVertexIHalfedgeIter cvih_iter(VertexHandle _vh) const
  { return ConstVertexIHalfedgeIter(*this, _vh); }
  /// const vertex - outgoing halfedge circulator
  ConstVertexOHalfedgeIter cvoh_iter(VertexHandle _vh) const
  { return ConstVertexOHalfedgeIter(*this, _vh); }
  /// const vertex - edge circulator
  ConstVertexEdgeIter cve_iter(VertexHandle _vh) const
  { return ConstVertexEdgeIter(*this, _vh); }
  /// const vertex - face circulator
  ConstVertexFaceIter cvf_iter(VertexHandle _vh) const
  { return ConstVertexFaceIter(*this, _vh); }

  /// face - vertex circulator
  FaceVertexIter fv_iter(FaceHandle _fh)
  { return FaceVertexIter(*this, _fh); }
  /// face - halfedge circulator
  FaceHalfedgeIter fh_iter(FaceHandle _fh)
  { return FaceHalfedgeIter(*this, _fh); }
  /// face - edge circulator
  FaceEdgeIter fe_iter(FaceHandle _fh)
  { return FaceEdgeIter(*this, _fh); }
  /// face - face circulator
  FaceFaceIter ff_iter(FaceHandle _fh)
  { return FaceFaceIter(*this, _fh); }

  /// const face - vertex circulator
  ConstFaceVertexIter cfv_iter(FaceHandle _fh) const
  { return ConstFaceVertexIter(*this, _fh); }
  /// const face - halfedge circulator
  ConstFaceHalfedgeIter cfh_iter(FaceHandle _fh) const
  { return ConstFaceHalfedgeIter(*this, _fh); }
  /// const face - edge circulator
  ConstFaceEdgeIter cfe_iter(FaceHandle _fh) const
  { return ConstFaceEdgeIter(*this, _fh); }
  /// const face - face circulator
  ConstFaceFaceIter cff_iter(FaceHandle _fh) const
  { return ConstFaceFaceIter(*this, _fh); }
  //@}

  /** \name Boundary and manifold tests
  */
  //@{
  bool is_boundary(HalfedgeHandle _heh) const
  { return ArrayKernel::is_boundary(_heh); }

  /** Is the edge _eh a boundary edge, i.e. is one of its halfedges
      a boundary halfege ? */
  bool is_boundary(EdgeHandle _eh) const
  {
    return (is_boundary(halfedge_handle(_eh, 0)) ||
            is_boundary(halfedge_handle(_eh, 1)));
  }
  /// Is vertex _vh a boundary vertex ?
  bool is_boundary(VertexHandle _vh) const
  {
    HalfedgeHandle heh(halfedge_handle(_vh));
    return (!(heh.is_valid() && face_handle(heh).is_valid()));
  }

  /** Is face _fh at boundary, i.e. is one of its edges (or vertices)
   *   a boundary edge?
   *  \param _fh Check this face
   *  \param _check_vertex If \c true, check the corner vertices of
   *         the face, too.
   */
  bool is_boundary(FaceHandle _fh, bool _check_vertex=false) const;
  /// Is (the mesh at) vertex _vh  two-manifold ?
  bool is_manifold(VertexHandle _vh) const;
  //@}

  // --- shortcuts ---
  
  /// returns the face handle of the opposite halfedge 
  inline FaceHandle opposite_face_handle(HalfedgeHandle _heh) const
  { return face_handle(opposite_halfedge_handle(_heh)); }
    
  // --- misc ---

  /** Adjust outgoing halfedge handle for vertices, so that it is a
      boundary halfedge whenever possible. 
  */
  void adjust_outgoing_halfedge(VertexHandle _vh);
  /// Find halfedge from _vh0 to _vh1. Returns invalid handle if not found.
  HalfedgeHandle find_halfedge(VertexHandle _start_vh, VertexHandle _end_vh) const;
  /// Vertex valence
  uint valence(VertexHandle _vh) const;
  /// Face valence
  uint valence(FaceHandle _fh) const;
  
  // --- connectivity operattions 
    
  /** Halfedge collapse: collapse the from-vertex of halfedge _heh
      into its to-vertex.

      \attention Needs vertex/edge/face status attribute in order to
      delete the items that degenerate.

      \note This function does not perform a garbage collection. It
      only marks degenerate items as deleted.

      \attention A halfedge collapse may lead to topological inconsistencies.
      Therefore you should check this using is_collapse_ok().  
      \TODO: implement is_collapse_ok() const for polygonal/triangle meshes
  */
  void collapse(HalfedgeHandle _heh);
  /** return true if the this the only link between the faces adjacent to _eh.
      _eh is allowed to be boundary, in which case true is returned iff _eh is 
      the only boundary edge of its ajdacent face.
  */
  bool is_simple_link(EdgeHandle _eh) const;
  /** return true if _fh shares only one edge with all of its adjacent faces.
      Boundary is treated as one face, i.e., the function false if _fh has more
      than one boundary edge.
  */
  bool is_simply_connected(FaceHandle _fh) const;
  /** Removes the edge _eh. Its adjacent faces are merged. _eh and one of the 
      adjacent faces are set deleted. The handle of the remaining face is 
      returned (InvalidFaceHandle is returned if _eh is a boundary edge).
      
      \precondition is_simple_link(_eh). This ensures that there are no hole faces
      or isolated vertices appearing in result of the operation.
      
      \attention Needs the Attributes::Status attribute for edges and faces.
      
      \note This function does not perform a garbage collection. It
      only marks items as deleted.
  */
  FaceHandle remove_edge(EdgeHandle _eh);
  /** Inverse of remove_edge. _eh should be the handle of the edge and the
      vertex and halfedge handles pointed by edge(_eh) should be valid. 
  */
  void reinsert_edge(EdgeHandle _eh);
  /** Inserts an edge between to_vh(_prev_heh) and from_vh(_next_heh).
      A new face is created started at heh0 of the inserted edge and
      its halfedges loop includes both _prev_heh and _next_heh. If an 
      old face existed which includes the argument halfedges, it is 
      split at the new edge. heh0 is returned. 
      
      \note assumes _prev_heh and _next_heh are either boundary or pointed
      to the same face
  */
  HalfedgeHandle insert_edge(HalfedgeHandle _prev_heh, HalfedgeHandle _next_heh);
    
  /// Face split (= 1-to-n split)
  void split(FaceHandle _fh, VertexHandle _vh);
  /// triangulate the face _fh
  void triangulate(FaceHandle _fh);
  /// triangulate the entire mesh
  void triangulate();
  

  /** \name Generic handle derefertiation.
      Calls the respective vertex(), halfedge(), edge(), face()
      method of the mesh kernel.
  */
  //@{
  /// Get item from handle
  const Vertex&    deref(VertexHandle _h)   const { return vertex(_h); }
  Vertex&          deref(VertexHandle _h)         { return vertex(_h); }
  const Halfedge&  deref(HalfedgeHandle _h) const { return halfedge(_h); }
  Halfedge&        deref(HalfedgeHandle _h)       { return halfedge(_h); }
  const Edge&      deref(EdgeHandle _h)     const { return edge(_h); }
  Edge&            deref(EdgeHandle _h)           { return edge(_h); }
  const Face&      deref(FaceHandle _h)     const { return face(_h); }
  Face&            deref(FaceHandle _h)           { return face(_h); }
  //@}

protected:  
  /// Helper for halfedge collapse
  void collapse_edge(HalfedgeHandle _hh);
  /// Helper for halfedge collapse
  void collapse_loop(HalfedgeHandle _hh);
};

}//namespace OpenMesh

#endif//OPENMESH_POLYCONNECTIVITY_HH
