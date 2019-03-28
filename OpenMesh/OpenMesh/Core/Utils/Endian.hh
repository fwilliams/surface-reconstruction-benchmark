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


#ifndef OPENMESH_UTILS_ENDIAN_HH
#define OPENMESH_UTILS_ENDIAN_HH


//== INCLUDES =================================================================


#include <OpenMesh/Core/System/config.h>
#include <iostream>


//== NAMESPACES ===============================================================


namespace OpenMesh {


//=============================================================================


/**  Determine byte order of host system.
 */
class Endian
{
public:
   
  enum Type {
    LSB = 1, ///< Little endian (Intel family and clones)
    MSB      ///< big endian (Motorola's 68x family, DEC Alpha, MIPS)
  };
  
  /// Return endian type of host system.
  static Type local() { return local_; }

  /// Return type _t as string.
  static const char * as_string(Type _t);
   
private:
   static int one_;
   static const Type local_;
};

//=============================================================================
} // namespace OpenMesh
//=============================================================================
#endif // OPENMESH_MESHREADER_HH defined
//=============================================================================

