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

#ifndef OPENMESH_ATTRIBKERNEL_HH
#define OPENMESH_ATTRIBKERNEL_HH


//== INCLUDES =================================================================

#include <OpenMesh/Core/Mesh/Attributes.hh>
#include <OpenMesh/Core/Utils/GenProg.hh>
#include <OpenMesh/Core/Utils/vector_traits.hh>
#include <vector>
#include <algorithm>

//== NAMESPACES ===============================================================

namespace OpenMesh {


//== CLASS DEFINITION =========================================================

/// This class adds the standard properties to the mesh type.
///
/// The attribute kernel adds all standard properties to the kernel. Therefore
/// the functions/types defined here provide a subset of the kernel
/// interface as described in Concepts::KernelT.
///
/// \see Concepts::KernelT
template <class MeshItems, class Connectivity>
class AttribKernelT : public Connectivity
{
public:

  //---------------------------------------------------------------- item types

  typedef typename Connectivity::Vertex     Vertex;
  typedef typename Connectivity::Halfedge   Halfedge;
  typedef typename Connectivity::Edge       Edge;
  typedef typename Connectivity::Face       Face;

  typedef typename MeshItems::Point         Point;
  typedef typename MeshItems::Normal        Normal;
  typedef typename MeshItems::Color         Color;
  typedef typename MeshItems::TexCoord1D    TexCoord1D;
  typedef typename MeshItems::TexCoord2D    TexCoord2D;
  typedef typename MeshItems::TexCoord3D    TexCoord3D;
  typedef typename MeshItems::Scalar        Scalar;

  typedef typename MeshItems::VertexData    VertexData;
  typedef typename MeshItems::HalfedgeData  HalfedgeData;
  typedef typename MeshItems::EdgeData      EdgeData;
  typedef typename MeshItems::FaceData      FaceData;

  typedef AttribKernelT<MeshItems,Connectivity>  AttribKernel;

  enum Attribs  {
    VAttribs = MeshItems::VAttribs,
    HAttribs = MeshItems::HAttribs,
    EAttribs = MeshItems::EAttribs,
    FAttribs = MeshItems::FAttribs,
  };

  typedef VPropHandleT<VertexData>              DataVPropHandle;
  typedef HPropHandleT<HalfedgeData>            DataHPropHandle;
  typedef EPropHandleT<EdgeData>                DataEPropHandle;
  typedef FPropHandleT<FaceData>                DataFPropHandle;

public:

  //-------------------------------------------------- constructor / destructor

  AttribKernelT()
  : refcount_vnormals_(0),
    refcount_vcolors_(0),
    refcount_vtexcoords1D_(0),
    refcount_vtexcoords2D_(0),
    refcount_vtexcoords3D_(0),
    refcount_htexcoords1D_(0),
    refcount_htexcoords2D_(0),
    refcount_htexcoords3D_(0),
    refcount_fnormals_(0),
    refcount_fcolors_(0)
  {
    add_property( points_, "v:points" );

    if (VAttribs & Attributes::Normal)
      request_vertex_normals();

    if (VAttribs & Attributes::Color)
      request_vertex_colors();

    if (VAttribs & Attributes::TexCoord1D)
      request_vertex_texcoords1D();

    if (VAttribs & Attributes::TexCoord2D)
      request_vertex_texcoords2D();

    if (VAttribs & Attributes::TexCoord3D)
      request_vertex_texcoords3D();

    if (HAttribs & Attributes::TexCoord1D)
      request_halfedge_texcoords1D();

    if (HAttribs & Attributes::TexCoord2D)
      request_halfedge_texcoords2D();

    if (HAttribs & Attributes::TexCoord3D)
      request_halfedge_texcoords3D();

    if (VAttribs & Attributes::Status)
      Connectivity::request_vertex_status();

    if (HAttribs & Attributes::Status)
      Connectivity::request_halfedge_status();

    if (EAttribs & Attributes::Status)
      Connectivity::request_edge_status();

    if (FAttribs & Attributes::Normal)
      request_face_normals();

    if (FAttribs & Attributes::Color)
      request_face_colors();

    if (FAttribs & Attributes::Status)
      Connectivity::request_face_status();
    //FIXME: data properties might actually cost storage even 
    //if there are no data traits??
    add_property(data_vpph_);
    add_property(data_fpph_);
    add_property(data_hpph_);
    add_property(data_epph_);
  }

  virtual ~AttribKernelT()
  {
    // should remove properties, but this will be done in
    // BaseKernel's destructor anyway...
  }

// Martin: This below does not make any sense, right?
  // -------------------------------------------------------- copy & assignment
//   AttribKernelT(const AttribKernelT& _rhs)
//     : BaseKernel(_rhs)
//   { operator=(_rhs); }
// 
//   AttribKernelT& operator=(const AttribKernelT& _rhs)
//   {
//     // remove old properties
//     remove_property(points_);
//     remove_property(vertex_normals_);
//     remove_property(vertex_colors_);
//     remove_property(vertex_texcoords_);
//     remove_property(vertex_status_);
//     remove_property(halfedge_status_);
//     remove_property(edge_status_);
//     remove_property(face_normals_);
//     remove_property(face_colors_);
//     remove_property(face_status_);
// 
//     // parent deep-copies properties
//     BaseKernel::operator=(_rhs);
// 
//     // copy property handles
//     points_            = _rhs.points_;
//     vertex_normals_    = _rhs.vertex_normals_;
//     vertex_colors_     = _rhs.vertex_colors_;
//     vertex_texcoords_  = _rhs.vertex_texcoords_;
//     vertex_status_     = _rhs.vertex_status_;
//     halfedge_status_   = _rhs.halfedge_status_;
//     edge_status_       = _rhs.edge_status_;
//     face_normals_      = _rhs.face_normals_;
//     face_colors_       = _rhs.face_colors_;
//     face_status_       = _rhs.face_status_;
// 
//     // copy ref-counts
//     refcount_vnormals_   = _rhs.refcount_vnormals_;
//     refcount_vcolors_    = _rhs.refcount_vcolors_;
//     refcount_vtexcoords_ = _rhs.refcount_vtexcoords_;
//     refcount_vstatus_    = _rhs.refcount_vstatus_;
//     refcount_hstatus_    = _rhs.refcount_hstatus_;
//     refcount_estatus_    = _rhs.refcount_estatus_;
//     refcount_fnormals_   = _rhs.refcount_fnormals_;
//     refcount_fcolors_    = _rhs.refcount_fcolors_;
//     refcount_fstatus_    = _rhs.refcount_fstatus_;
//
//     return *this;
//   }

  /** Assignment from another mesh of \em another type.
      \note All that's copied is connectivity and vertex positions.
      All other information (like e.g. attributes or additional
      elements from traits classes) is not copied.
      \note If you want to copy all information, including *custom* properties,
      use PolyMeshT::operator=() instead.  
      TODO: version which copies standard properties specified by the user
  */
  template <class _AttribKernel>
  void assign(const _AttribKernel& _other)
  {
    assign_connectivity(_other);
    for (typename Connectivity::VertexIter v_it = Connectivity::vertices_begin(); 
         v_it != Connectivity::vertices_end(); ++v_it)
    {//assumes Point constructor supports cast from _AttribKernel::Point
      set_point(v_it, (Point)_other.point(v_it));
    }    
  }

  //-------------------------------------------------------------------- points

  const Point* points() const 
  { return property(points_).data(); }

  const Point& point(VertexHandle _vh) const 
  { return property(points_, _vh); }

  Point& point(VertexHandle _vh) 
  { return property(points_, _vh); }

  void set_point(VertexHandle _vh, const Point& _p) 
  { property(points_, _vh) = _p; }


  //------------------------------------------------------------ vertex normals

  const Normal* vertex_normals() const 
  { return property(vertex_normals_).data(); }

  const Normal& normal(VertexHandle _vh) const 
  { return property(vertex_normals_, _vh); }

  void set_normal(VertexHandle _vh, const Normal& _n) 
  { property(vertex_normals_, _vh) = _n; }


  //------------------------------------------------------------- vertex colors

  const Color* vertex_colors() const 
  { return property(vertex_colors_).data(); }

  const Color& color(VertexHandle _vh) const
  { return property(vertex_colors_, _vh); }

  void set_color(VertexHandle _vh, const Color& _c) 
  { property(vertex_colors_, _vh) = _c; }


  //------------------------------------------------------- vertex 1D texcoords

  const TexCoord1D* texcoords1D() const {
    return property(vertex_texcoords1D_).data();
  }

  const TexCoord1D& texcoord1D(VertexHandle _vh) const {
    return property(vertex_texcoords1D_, _vh);
  }

  void set_texcoord1D(VertexHandle _vh, const TexCoord1D& _t) {
    property(vertex_texcoords1D_, _vh) = _t;
  }


  //------------------------------------------------------- vertex 2D texcoords

  const TexCoord2D* texcoords2D() const {
    return property(vertex_texcoords2D_).data();
  }

  const TexCoord2D& texcoord2D(VertexHandle _vh) const {
    return property(vertex_texcoords2D_, _vh);
  }

  void set_texcoord2D(VertexHandle _vh, const TexCoord2D& _t) {
    property(vertex_texcoords2D_, _vh) = _t;
  }


  //------------------------------------------------------- vertex 3D texcoords

  const TexCoord3D* texcoords3D() const {
    return property(vertex_texcoords3D_).data();
  }

  const TexCoord3D& texcoord3D(VertexHandle _vh) const {
    return property(vertex_texcoords3D_, _vh);
  }

  void set_texcoord3D(VertexHandle _vh, const TexCoord3D& _t) {
    property(vertex_texcoords3D_, _vh) = _t;
  }
 
  //.------------------------------------------------------ halfedge 1D texcoords
  
  const TexCoord1D* htexcoords1D() const {
    return property(halfedge_texcoords1D_).data();
  }

  const TexCoord1D& texcoord1D(HalfedgeHandle _heh) const {
    return property(halfedge_texcoords1D_, _heh);
  }

  void set_texcoord1D(HalfedgeHandle _heh, const TexCoord1D& _t) {
    property(halfedge_texcoords1D_, _heh) = _t;
  }


  //------------------------------------------------------- halfedge 2D texcoords

  const TexCoord2D* htexcoords2D() const {
    return property(halfedge_texcoords2D_).data();
  }

  const TexCoord2D& texcoord2D(HalfedgeHandle _heh) const {
    return property(halfedge_texcoords2D_, _heh);
  }

  void set_texcoord2D(HalfedgeHandle _heh, const TexCoord2D& _t) {
    property(halfedge_texcoords2D_, _heh) = _t;
  }


  //------------------------------------------------------- halfedge 3D texcoords

  const TexCoord3D* htexcoords3D() const {
    return property(halfedge_texcoords3D_).data();
  }

  const TexCoord3D& texcoord3D(HalfedgeHandle _heh) const {
    return property(halfedge_texcoords3D_, _heh);
  }

  void set_texcoord3D(HalfedgeHandle _heh, const TexCoord3D& _t) {
    property(halfedge_texcoords3D_, _heh) = _t;
  }
  
  //-------------------------------------------------------------- face normals

  const Normal& normal(FaceHandle _fh) const 
  { return property(face_normals_, _fh); }

  void set_normal(FaceHandle _fh, const Normal& _n) 
  { property(face_normals_, _fh) = _n; }

  //--------------------------------------------------------------- face colors

  const Color& color(FaceHandle _fh) const
  { return property(face_colors_, _fh); }

  void set_color(FaceHandle _fh, const Color& _c)
  { property(face_colors_, _fh) = _c; }

  //------------------------------------------------ request / alloc properties

  void request_vertex_normals()
  {
    if (!refcount_vnormals_++)
      add_property( vertex_normals_, "v:normals" );
  }

  void request_vertex_colors()
  {
    if (!refcount_vcolors_++)
      add_property( vertex_colors_, "v:colors" );
  }

  void request_vertex_texcoords1D() 
  {
    if (!refcount_vtexcoords1D_++)
      add_property( vertex_texcoords1D_, "v:texcoords1D" );
  }

  void request_vertex_texcoords2D() 
  {
    if (!refcount_vtexcoords2D_++)
      add_property( vertex_texcoords2D_, "v:texcoords2D" );
  }

  void request_vertex_texcoords3D() 
  {
    if (!refcount_vtexcoords3D_++)
      add_property( vertex_texcoords3D_, "v:texcoords3D" );
  }

  void request_halfedge_texcoords1D() 
  {
    if (!refcount_htexcoords1D_++)
      add_property( halfedge_texcoords1D_, "h:texcoords1D" );
  }

  void request_halfedge_texcoords2D() 
  {
    if (!refcount_htexcoords2D_++)
      add_property( halfedge_texcoords2D_, "h:texcoords2D" );
  }

  void request_halfedge_texcoords3D() 
  {
    if (!refcount_htexcoords3D_++)
      add_property( halfedge_texcoords3D_, "h:texcoords3D" );
  }

  void request_face_normals()
  {
    if (!refcount_fnormals_++)
      add_property( face_normals_, "f:normals" );
  }

  void request_face_colors()
  {
    if (!refcount_fcolors_++)
      add_property( face_colors_, "f:colors" );
  }

  //------------------------------------------------- release / free properties

  void release_vertex_normals()
  {
    if ((refcount_vnormals_ > 0) && (! --refcount_vnormals_))
      remove_property(vertex_normals_);
  }

  void release_vertex_colors()
  {
    if ((refcount_vcolors_ > 0) && (! --refcount_vcolors_))
      remove_property(vertex_colors_);
  }

  void release_vertex_texcoords1D() {
    if ((refcount_vtexcoords1D_ > 0) && (! --refcount_vtexcoords1D_))
      remove_property(vertex_texcoords1D_);
  }

  void release_vertex_texcoords2D() {
    if ((refcount_vtexcoords2D_ > 0) && (! --refcount_vtexcoords2D_))
      remove_property(vertex_texcoords2D_);
  }

  void release_vertex_texcoords3D() {
    if ((refcount_vtexcoords3D_ > 0) && (! --refcount_vtexcoords3D_))
      remove_property(vertex_texcoords3D_);
  }

  void release_halfedge_texcoords1D() {
    if ((refcount_htexcoords1D_ > 0) && (! --refcount_htexcoords1D_))
      remove_property(halfedge_texcoords1D_);
  }

  void release_halfedge_texcoords2D() {
    if ((refcount_htexcoords2D_ > 0) && (! --refcount_htexcoords2D_))
      remove_property(halfedge_texcoords2D_);
  }

  void release_halfedge_texcoords3D() {
    if ((refcount_htexcoords3D_ > 0) && (! --refcount_htexcoords3D_))
      remove_property(halfedge_texcoords3D_);
  }

  void release_face_normals()
  {
    if ((refcount_fnormals_ > 0) && (! --refcount_fnormals_))
      remove_property(face_normals_);
  }

  void release_face_colors()
  {
    if ((refcount_fcolors_ > 0) && (! --refcount_fcolors_))
      remove_property(face_colors_);
  }

  //---------------------------------------------- dynamic check for properties

  bool has_vertex_normals()   const { return vertex_normals_.is_valid();   }
  bool has_vertex_colors()    const { return vertex_colors_.is_valid();    }
  bool has_vertex_texcoords1D() const { return vertex_texcoords1D_.is_valid();}
  bool has_vertex_texcoords2D() const { return vertex_texcoords2D_.is_valid();}
  bool has_vertex_texcoords3D() const { return vertex_texcoords3D_.is_valid();}
  bool has_halfedge_texcoords1D() const { return halfedge_texcoords1D_.is_valid();}
  bool has_halfedge_texcoords2D() const { return halfedge_texcoords2D_.is_valid();}
  bool has_halfedge_texcoords3D() const { return halfedge_texcoords3D_.is_valid();}
  bool has_face_normals()     const { return face_normals_.is_valid();     }
  bool has_face_colors()      const { return face_colors_.is_valid();      }

public:

  typedef VPropHandleT<Point>               PointsPropertyHandle;
  typedef VPropHandleT<Normal>              VertexNormalsPropertyHandle;
  typedef VPropHandleT<Color>               VertexColorsPropertyHandle;
  typedef VPropHandleT<TexCoord1D>          VertexTexCoords1DPropertyHandle;
  typedef VPropHandleT<TexCoord2D>          VertexTexCoords2DPropertyHandle;
  typedef VPropHandleT<TexCoord3D>          VertexTexCoords3DPropertyHandle;
  typedef HPropHandleT<TexCoord1D>          HalfedgeTexCoords1DPropertyHandle;
  typedef HPropHandleT<TexCoord2D>          HalfedgeTexCoords2DPropertyHandle;
  typedef HPropHandleT<TexCoord3D>          HalfedgeTexCoords3DPropertyHandle;
  typedef FPropHandleT<Normal>              FaceNormalsPropertyHandle;
  typedef FPropHandleT<Color>               FaceColorsPropertyHandle;

public:
  //standard vertex properties
  PointsPropertyHandle                      points_pph() const
  { return points_; }

  VertexNormalsPropertyHandle               vertex_normals_pph() const
  { return vertex_normals_; }

  VertexColorsPropertyHandle                vertex_colors_pph() const
  { return vertex_colors_; }

  VertexTexCoords1DPropertyHandle           vertex_texcoords1D_pph() const
  { return vertex_texcoords1D_; }

  VertexTexCoords2DPropertyHandle           vertex_texcoords2D_pph() const
  { return vertex_texcoords2D_; }

  VertexTexCoords3DPropertyHandle           vertex_texcoords3D_pph() const
  { return vertex_texcoords3D_; }

  //standard halfedge properties
  HalfedgeTexCoords1DPropertyHandle           halfedge_texcoords1D_pph() const
  { return halfedge_texcoords1D_; }

  HalfedgeTexCoords2DPropertyHandle           halfedge_texcoords2D_pph() const
  { return halfedge_texcoords2D_; }

  HalfedgeTexCoords3DPropertyHandle           halfedge_texcoords3D_pph() const
  { return halfedge_texcoords3D_; }
 
	//standard face properties
  FaceNormalsPropertyHandle                 face_normals_pph() const
  { return face_normals_; }

  FaceColorsPropertyHandle                  face_colors_pph() const
  { return face_colors_; }

  VertexData&                               data(VertexHandle _vh)
  { return property(data_vpph_, _vh); }

  const VertexData&                         data(VertexHandle _vh) const
  { return property(data_vpph_, _vh); }

  FaceData&                                 data(FaceHandle _fh)
  { return property(data_fpph_, _fh); }

  const FaceData&                           data(FaceHandle _fh) const
  { return property(data_fpph_, _fh); }

  EdgeData&                                 data(EdgeHandle _eh)
  { return property(data_epph_, _eh); }

  const EdgeData&                           data(EdgeHandle _eh) const
  { return property(data_epph_, _eh); }

  HalfedgeData&                             data(HalfedgeHandle _heh)
  { return property(data_hpph_, _heh); }

  const HalfedgeData&                       data(HalfedgeHandle _heh) const
  { return property(data_hpph_, _heh); }

private:
  //standard vertex properties
  PointsPropertyHandle                      points_;
  VertexNormalsPropertyHandle               vertex_normals_;
  VertexColorsPropertyHandle                vertex_colors_;
  VertexTexCoords1DPropertyHandle           vertex_texcoords1D_;
  VertexTexCoords2DPropertyHandle           vertex_texcoords2D_;
  VertexTexCoords3DPropertyHandle           vertex_texcoords3D_;
  //standard halfedge properties
  HalfedgeTexCoords1DPropertyHandle         halfedge_texcoords1D_;
  HalfedgeTexCoords2DPropertyHandle         halfedge_texcoords2D_;
  HalfedgeTexCoords3DPropertyHandle         halfedge_texcoords3D_;
  //standard face properties
  FaceNormalsPropertyHandle                 face_normals_;
  FaceColorsPropertyHandle                  face_colors_;
  //data properties handles
  DataVPropHandle                           data_vpph_;
  DataHPropHandle                           data_hpph_;
  DataEPropHandle                           data_epph_;
  DataFPropHandle                           data_fpph_;

  unsigned int                              refcount_vnormals_;
  unsigned int                              refcount_vcolors_;
  unsigned int                              refcount_vtexcoords1D_;
  unsigned int                              refcount_vtexcoords2D_;
  unsigned int                              refcount_vtexcoords3D_;
  unsigned int                              refcount_htexcoords1D_;
  unsigned int                              refcount_htexcoords2D_;
  unsigned int                              refcount_htexcoords3D_;
  unsigned int                              refcount_fnormals_;
  unsigned int                              refcount_fcolors_;
};

//=============================================================================
} // namespace OpenMesh
//=============================================================================
#endif // OPENMESH_ATTRIBKERNEL_HH defined
//=============================================================================
