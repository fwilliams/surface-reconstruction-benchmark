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
//  CLASS PolyMesh_ArrayKernelT
//
//=============================================================================


#ifndef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
#define OPENMESH_POLY_MESH_ARRAY_KERNEL_HH


//== INCLUDES =================================================================


#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/Mesh/FinalMeshItemsT.hh>
#include <OpenMesh/Core/Mesh/AttribKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMeshT.hh>


//== NAMESPACES ===============================================================


namespace OpenMesh {


//== CLASS DEFINITION =========================================================

/// Helper class to build a PolyMesh-type
template <class Traits>
struct PolyMesh_ArrayKernel_GeneratorT
{
  typedef FinalMeshItemsT<Traits, false>              MeshItems;
  typedef AttribKernelT<MeshItems, PolyConnectivity>  AttribKernel;
  typedef PolyMeshT<AttribKernel>                     Mesh;
};


/** \class PolyMesh_ArrayKernelT PolyMesh_ArrayKernelT.hh <OpenMesh/Mesh/Types/PolyMesh_ArrayKernelT.hh>

    \ingroup mesh_types_group
    Polygonal mesh based on the ArrayKernel.
    \see OpenMesh::PolyMeshT
    \see OpenMesh::ArrayKernel
*/
template <class Traits = DefaultTraits>
class PolyMesh_ArrayKernelT
  : public PolyMesh_ArrayKernel_GeneratorT<Traits>::Mesh
{};


//=============================================================================
} // namespace OpenMesh
//=============================================================================
#endif // OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
//=============================================================================
