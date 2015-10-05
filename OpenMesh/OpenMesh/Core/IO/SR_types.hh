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
//  Helper Functions for binary reading / writing
//
//=============================================================================

#ifndef OPENMESH_SR_TYPES_HH
#define OPENMESH_SR_TYPES_HH


//== INCLUDES =================================================================

#include <OpenMesh/Core/System/config.h>


//== NAMESPACES ===============================================================

namespace OpenMesh {
namespace IO {


//=============================================================================


/** \name Handling binary input/output.
    These functions take care of swapping bytes to get the right Endian.
*/
//@{

//-----------------------------------------------------------------------------

typedef unsigned char    uchar;
typedef unsigned short   ushort;
typedef unsigned long    ulong;

typedef char           int8_t;  typedef unsigned char      uint8_t;
typedef short          int16_t; typedef unsigned short     uint16_t;
typedef long           int32_t; typedef unsigned long      uint32_t;
#if defined(OM_CC_MSVC)
typedef __int64        int64_t; typedef unsigned __int64   uint64_t;
#else
typedef long long      int64_t; typedef unsigned long long uint64_t;
#endif

typedef float          float32_t;
typedef double         float64_t;

typedef uint8_t        rgb_t[3];
typedef uint8_t        rgba_t[4];
   
//@}


//=============================================================================
} // namespace IO
} // namespace OpenMesh
//=============================================================================
#endif // OPENMESH_MESHREADER_HH defined
//=============================================================================

