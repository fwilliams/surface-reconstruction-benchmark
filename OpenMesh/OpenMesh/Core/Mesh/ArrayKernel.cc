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

#include <OpenMesh/Core/Mesh/ArrayKernel.hh>

namespace OpenMesh
{

ArrayKernel::ArrayKernel()
: refcount_vstatus_(0), refcount_hstatus_(0),
  refcount_estatus_(0), refcount_fstatus_(0)
{
  init_bit_masks(); //Status bit masks initialization
}

ArrayKernel::~ArrayKernel()
{
  clear();
}

// ArrayKernel::ArrayKernel(const ArrayKernel& _rhs)
// : BaseKernel(_rhs),
//   vertices_(_rhs.vertices_), edges_(_rhs.edges_), faces_(_rhs.faces_),
//   vertex_status_(_rhs.vertex_status_), halfedge_status_(_rhs.halfedge_status_),
//   edge_status_(_rhs.edge_status_), face_status_(_rhs.face_status_),
//   refcount_vstatus_(_rhs.refcount_vstatus_), refcount_hstatus_(_rhs.refcount_hstatus_),
//   refcount_estatus_(_rhs.refcount_estatus_), refcount_fstatus_(_rhs.refcount_fstatus_)
// {}


void ArrayKernel::assign_connectivity(const ArrayKernel& _other)
{
  vertices_ = _other.vertices_;
  edges_ = _other.edges_;
  faces_ = _other.faces_;
  
  vprops_resize(n_vertices());
  hprops_resize(n_halfedges());
  eprops_resize(n_edges());
  fprops_resize(n_faces());
  
#define COPY_STATUS_PROPERTY(ENTITY) \
  if (_other.ENTITY##_status_.is_valid()) \
  {   \
    if (!ENTITY##_status_.is_valid()) \
    { \
      request_##ENTITY##_status(); \
    } \
    property(ENTITY##_status_) = _other.property(_other.ENTITY##_status_); \
  }
  COPY_STATUS_PROPERTY(vertex)
  COPY_STATUS_PROPERTY(halfedge)
  COPY_STATUS_PROPERTY(edge)
  COPY_STATUS_PROPERTY(face)
  
#undef COPY_STATUS_PROPERTY
}

uint ArrayKernel::delete_isolated_vertices()
{
  assert(has_vertex_status());//this function requires vertex status property
  uint n_isolated = 0;
  for (KernelVertexIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it)
  {
    if (is_isolated(handle(*v_it)))
    {
      status(handle(*v_it)).set_deleted(true);
      n_isolated++;
    }
  }
  return n_isolated;
}

void ArrayKernel::garbage_collection(bool _v, bool _e, bool _f)
{
  int i, i0, i1, nV(n_vertices()), nE(n_edges()), nH(2*n_edges()), nF(n_faces());

  std::vector<VertexHandle>    vh_map;
  std::vector<HalfedgeHandle>  hh_map;
  std::vector<FaceHandle>      fh_map;

  // setup handle mapping:
  vh_map.reserve(nV);
  for (i=0; i<nV; ++i) vh_map.push_back(VertexHandle(i));

  hh_map.reserve(nH);
  for (i=0; i<nH; ++i) hh_map.push_back(HalfedgeHandle(i));

  fh_map.reserve(nF);
  for (i=0; i<nF; ++i) fh_map.push_back(FaceHandle(i));

  // remove deleted vertices
  if (_v && n_vertices() > 0)
  {
    i0=0;  i1=nV-1;

    while (1)
    {
      // find 1st deleted and last un-deleted
      while (!status(VertexHandle(i0)).deleted() && i0 < i1)  ++i0;
      while ( status(VertexHandle(i1)).deleted() && i0 < i1)  --i1;
      if (i0 >= i1) break;

      // swap
      std::swap(vertices_[i0], vertices_[i1]);
      std::swap(vh_map[i0],  vh_map[i1]);
      vprops_swap(i0, i1);
    };

    vertices_.resize(status(VertexHandle(i0)).deleted() ? i0 : i0+1);
    vprops_resize(n_vertices());
  }


  // remove deleted edges
  if (_e && n_edges() > 0)
  {
    i0=0;  i1=nE-1;

    while (1)
    {
      // find 1st deleted and last un-deleted
      while (!status(EdgeHandle(i0)).deleted() && i0 < i1)  ++i0;
      while ( status(EdgeHandle(i1)).deleted() && i0 < i1)  --i1;
      if (i0 >= i1) break;

      // swap
      std::swap(edges_[i0], edges_[i1]);
      std::swap(hh_map[2*i0], hh_map[2*i1]);
      std::swap(hh_map[2*i0+1], hh_map[2*i1+1]);
      eprops_swap(i0, i1);
      hprops_swap(2*i0,   2*i1);
      hprops_swap(2*i0+1, 2*i1+1);
    };

    edges_.resize(status(EdgeHandle(i0)).deleted() ? i0 : i0+1);
    eprops_resize(n_edges());
    hprops_resize(n_halfedges());
  }


  // remove deleted faces
  if (_f && n_faces() > 0)
  {
    i0=0;  i1=nF-1;

    while (1)
    {
      // find 1st deleted and last un-deleted
      while (!status(FaceHandle(i0)).deleted() && i0 < i1)  ++i0;
      while ( status(FaceHandle(i1)).deleted() && i0 < i1)  --i1;
      if (i0 >= i1) break;

      // swap
      std::swap(faces_[i0], faces_[i1]);
      std::swap(fh_map[i0], fh_map[i1]);
      fprops_swap(i0, i1);
    };

    faces_.resize(status(FaceHandle(i0)).deleted() ? i0 : i0+1);
    fprops_resize(n_faces());
  }


  // update handles of vertices
  if (_e)
  {
    KernelVertexIter v_it(vertices_begin()), v_end(vertices_end());
    VertexHandle     vh;

    for (; v_it!=v_end; ++v_it)
    {
      vh = handle(*v_it);
      if (!is_isolated(vh))
      {
        set_halfedge_handle(vh, hh_map[halfedge_handle(vh).idx()]);
      }
    }
  }

  HalfedgeHandle hh;
  // update handles of halfedges
  for (KernelEdgeIter e_it(edges_begin()); e_it != edges_end(); ++e_it)
  {//in the first pass update the (half)edges vertices
    hh = halfedge_handle(handle(*e_it), 0);
    set_vertex_handle(hh, vh_map[to_vertex_handle(hh).idx()]);
    hh = halfedge_handle(handle(*e_it), 1);
    set_vertex_handle(hh, vh_map[to_vertex_handle(hh).idx()]);
  }
  for (KernelEdgeIter e_it(edges_begin()); e_it != edges_end(); ++e_it)
  {//in the second pass update the connectivity of the (half)edges
    hh = halfedge_handle(handle(*e_it), 0);
    set_next_halfedge_handle(hh, hh_map[next_halfedge_handle(hh).idx()]);
    if (!is_boundary(hh))
    {
      set_face_handle(hh, fh_map[face_handle(hh).idx()]);
    }
    hh = halfedge_handle(handle(*e_it), 1);
    set_next_halfedge_handle(hh, hh_map[next_halfedge_handle(hh).idx()]);
    if (!is_boundary(hh))
    {
      set_face_handle(hh, fh_map[face_handle(hh).idx()]);
    }
  }

  // update handles of faces
  if (_e)
  {
    KernelFaceIter  f_it(faces_begin()), f_end(faces_end());
    FaceHandle      fh;

    for (; f_it!=f_end; ++f_it)
    {
      fh = handle(*f_it);
      set_halfedge_handle(fh, hh_map[halfedge_handle(fh).idx()]);
    }
  }
}

void ArrayKernel::clear()
{
  vertices_.clear();
  edges_.clear();
  faces_.clear();

  vprops_resize(0);
  eprops_resize(0);
  hprops_resize(0);
  fprops_resize(0);
}

void ArrayKernel::resize( uint _n_vertices, uint _n_edges, uint _n_faces )
{
  vertices_.resize(_n_vertices);
  edges_.resize(_n_edges);
  faces_.resize(_n_faces);

  vprops_resize(n_vertices());
  hprops_resize(n_halfedges());
  eprops_resize(n_edges());
  fprops_resize(n_faces());
}

void ArrayKernel::reserve(uint _n_vertices, uint _n_edges, uint _n_faces )
{
  vertices_.reserve(_n_vertices);
  edges_.reserve(_n_edges);
  faces_.reserve(_n_faces);

  vprops_reserve(_n_vertices);
  hprops_reserve(_n_edges*2);
  eprops_reserve(_n_edges);
  fprops_reserve(_n_faces);
}

// Status Sets API
void ArrayKernel::init_bit_masks(BitMaskContainer& _bmc)
{
  for (uint i = Attributes::UNUSED; i != 0; i <<= 1)
  {
    _bmc.push_back(i);
  }
}

void ArrayKernel::init_bit_masks()
{
  init_bit_masks(vertex_bit_masks_);
  edge_bit_masks_ = vertex_bit_masks_;//init_bit_masks(edge_bit_masks_);
  face_bit_masks_ = vertex_bit_masks_;//init_bit_masks(face_bit_masks_);
  halfedge_bit_masks_= vertex_bit_masks_;//init_bit_masks(halfedge_bit_masks_);
}


};

