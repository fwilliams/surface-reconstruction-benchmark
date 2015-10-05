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
//  CLASS TriMesh_ArrayKernelT
//
//=============================================================================


#ifndef OPENMESH_TRIMESH_ARRAY_KERNEL_HH
#define OPENMESH_TRIMESH_ARRAY_KERNEL_HH


//== INCLUDES =================================================================


#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/TriConnectivity.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/Mesh/FinalMeshItemsT.hh>
#include <OpenMesh/Core/Mesh/AttribKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMeshT.hh>


//== NAMESPACES ===============================================================


namespace OpenMesh {


//== CLASS DEFINITION =========================================================


/// Helper class to create a TriMesh-type based on ArrayKernelT
template <class Traits>
struct TriMesh_ArrayKernel_GeneratorT
{
  typedef FinalMeshItemsT<Traits, true>               MeshItems;
  typedef AttribKernelT<MeshItems, TriConnectivity>   AttribKernel;
  typedef TriMeshT<AttribKernel>                      Mesh;
};



/** \ingroup mesh_types_group
    Triangle mesh based on the ArrayKernel.
    \see OpenMesh::TriMeshT
    \see OpenMesh::ArrayKernelT
*/
template <class Traits = DefaultTraits>
class TriMesh_ArrayKernelT
  : public TriMesh_ArrayKernel_GeneratorT<Traits>::Mesh
{};


//=============================================================================
} // namespace OpenMesh
//=============================================================================
#endif // OPENMESH_TRIMESH_ARRAY_KERNEL_HH
//=============================================================================
