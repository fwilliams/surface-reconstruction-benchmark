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
//   by the Free Software Foundation, version 2.
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
//  CLASS PolyMeshT - IMPLEMENTATION
//
//=============================================================================


#define OPENMESH_POLYMESH_C


//== INCLUDES =================================================================

#include <OpenMesh/Core/Mesh/PolyMeshT.hh>
#include <OpenMesh/Core/Geometry/LoopSchemeMaskT.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/System/omstream.hh>
#include <vector>


//== NAMESPACES ===============================================================


namespace OpenMesh {

//== IMPLEMENTATION ==========================================================

template <class Kernel>
uint PolyMeshT<Kernel>::find_feature_edges(Scalar _angle_tresh)
{
  assert(Kernel::has_edge_status());//this function needs edge status property
  uint n_feature_edges = 0;
  for (EdgeIter e_it = Kernel::edges_begin(); e_it != Kernel::edges_end(); ++e_it)
  {
    if (fabs(calc_dihedral_angle(e_it)) > _angle_tresh)
    {//note: could be optimized by comparing cos(dih_angle) vs. cos(_angle_tresh)
      status(e_it).set_feature(true);
      n_feature_edges++;
    }
    else
    {
      status(e_it).set_feature(false);
    }
  }
  return n_feature_edges;
}

//-----------------------------------------------------------------------------

template <class Kernel>
typename PolyMeshT<Kernel>::Normal
PolyMeshT<Kernel>::
calc_face_normal(FaceHandle _fh) const
{
  assert(halfedge_handle(_fh).is_valid());
  ConstFaceVertexIter fv_it(cfv_iter(_fh));

  const Point& p0(point(fv_it));  ++fv_it;
  const Point& p1(point(fv_it));  ++fv_it;
  const Point& p2(point(fv_it));

  return calc_face_normal(p0, p1, p2);
}


//-----------------------------------------------------------------------------


template <class Kernel>
typename PolyMeshT<Kernel>::Normal
PolyMeshT<Kernel>::
calc_face_normal(const Point& _p0,
     const Point& _p1,
     const Point& _p2) const
{
#if 1
  // The OpenSG <Vector>::operator -= () does not support the type Point
  // as rhs. Therefore use vector_cast at this point!!!
  // Note! OpenSG distinguishes between Normal and Point!!!
  Normal p1p0(_p0);  p1p0 -= vector_cast<Normal>(_p1);
  Normal p1p2(_p2);  p1p2 -= vector_cast<Normal>(_p1);

  Normal n    = cross(p1p2, p1p0);
  Scalar norm = n.length();

  // The expression ((n *= (1.0/norm)),n) is used because the OpenSG
  // vector class does not return self after component-wise
  // self-multiplication with a scalar!!!
  return (norm != Scalar(0)) ? ((n *= (Scalar(1)/norm)),n) : Normal(0,0,0);
#else
  Point p1p0 = _p0;  p1p0 -= _p1;
  Point p1p2 = _p2;  p1p2 -= _p1;

  Normal n = vector_cast<Normal>(cross(p1p2, p1p0));
  Scalar norm = n.length();

  return (norm != 0.0) ? n *= (1.0/norm) : Normal(0,0,0);
#endif
}

//-----------------------------------------------------------------------------

template <class Kernel>
void
PolyMeshT<Kernel>::
calc_face_centroid(FaceHandle _fh, Point& _pt) const
{
  _pt.vectorize(0);
  uint valence = 0;
  for (ConstFaceVertexIter cfv_it = cfv_iter(_fh); cfv_it; ++cfv_it, ++valence)
  {
    _pt += point(cfv_it);
  }
  _pt /= valence;
}
//-----------------------------------------------------------------------------


template <class Kernel>
void
PolyMeshT<Kernel>::
update_normals()
{
  if (Kernel::has_face_normals())    update_face_normals();
  if (Kernel::has_vertex_normals())  update_vertex_normals();
}


//-----------------------------------------------------------------------------


template <class Kernel>
void
PolyMeshT<Kernel>::
update_face_normals()
{
  FaceIter f_it(Kernel::faces_begin()), f_end(Kernel::faces_end());

  for (; f_it != f_end; ++f_it)
    set_normal(f_it.handle(), calc_face_normal(f_it.handle()));
}


//-----------------------------------------------------------------------------


template <class Kernel>
typename PolyMeshT<Kernel>::Normal
PolyMeshT<Kernel>::
calc_vertex_normal(VertexHandle _vh) const
{
  Normal n;
  calc_vertex_normal_fast(_vh,n);

  Scalar norm = n.length();
  if (norm != 0.0) n *= (1.0/norm);

  return n;
}

//-----------------------------------------------------------------------------
template <class Kernel>
void PolyMeshT<Kernel>::
calc_vertex_normal_fast(VertexHandle _vh, Normal& _n) const
{
  _n.vectorize(0.0);
  for (ConstVertexFaceIter vf_it=cvf_iter(_vh); vf_it; ++vf_it)
    _n += normal(vf_it.handle());
}

//-----------------------------------------------------------------------------
template <class Kernel>
void PolyMeshT<Kernel>::
calc_vertex_normal_correct(VertexHandle _vh, Normal& _n) const
{
  _n.vectorize(0.0);
  ConstVertexIHalfedgeIter cvih_it = cvih_iter(_vh);
  if (!cvih_it)
  {//don't crash on isolated vertices
    return;
  }
  Normal in_he_vec;
  calc_edge_vector(cvih_it, in_he_vec);
  for ( ; cvih_it; ++cvih_it)
  {//calculates the sector normal defined by cvih_it and adds it to _n
    if (is_boundary(cvih_it))
    {
      continue;
    }
    HalfedgeHandle out_heh(next_halfedge_handle(cvih_it));
    Normal out_he_vec;
    calc_edge_vector(out_heh, out_he_vec);
    _n += cross(in_he_vec, out_he_vec);//sector area is taken into account
    in_he_vec = out_he_vec;
    in_he_vec *= -1;//change the orientation
  }
}

//-----------------------------------------------------------------------------
template <class Kernel>
void PolyMeshT<Kernel>::
calc_vertex_normal_loop(VertexHandle _vh, Normal& _n) const
{
  static const LoopSchemeMaskDouble& loop_scheme_mask__ =
                  LoopSchemeMaskDoubleSingleton::Instance();

  Normal t_v(0.0,0.0,0.0), t_w(0.0,0.0,0.0);
  unsigned int vh_val = valence(_vh);
  unsigned int i = 0;
  for (ConstVertexOHalfedgeIter cvoh_it = cvoh_iter(_vh); cvoh_it; ++cvoh_it, ++i)
  {
    VertexHandle r1_v(to_vertex_handle(cvoh_it));
    t_v += (typename Point::value_type)(loop_scheme_mask__.tang0_weight(vh_val, i))*point(r1_v);
    t_w += (typename Point::value_type)(loop_scheme_mask__.tang1_weight(vh_val, i))*point(r1_v);
  }
  _n = cross(t_w, t_v);//hack: should be cross(t_v, t_w), but then the normals are reversed?
}

//-----------------------------------------------------------------------------


template <class Kernel>
void
PolyMeshT<Kernel>::
update_vertex_normals()
{
  VertexIter  v_it(Kernel::vertices_begin()), v_end(Kernel::vertices_end());

  for (; v_it!=v_end; ++v_it)
    set_normal(v_it.handle(), calc_vertex_normal(v_it.handle()));
}

//=============================================================================
} // namespace OpenMesh
//=============================================================================
