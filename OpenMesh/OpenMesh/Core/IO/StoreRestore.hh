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

#ifndef OPENMESH_STORERESTORE_HH
#define OPENMESH_STORERESTORE_HH


//== INCLUDES =================================================================

#include <stdexcept>
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/IO/SR_binary.hh>

//== NAMESPACES ===============================================================

namespace OpenMesh {
namespace IO {


//=============================================================================


/** \name Handling binary input/output.
    These functions take care of swapping bytes to get the right Endian.
*/
//@{


//-----------------------------------------------------------------------------
// StoreRestore definitions

template <typename T> inline
bool is_streamable(void)
{ return binary< T >::is_streamable; }

template <typename T> inline
bool is_streamable( const T& ) 
{ return binary< T >::is_streamable; }

template <typename T> inline
size_t size_of( const T& _v ) 
{ return binary< T >::size_of(_v); }

template <typename T> inline
size_t size_of(void) 
{ return binary< T >::size_of(); }

template <typename T> inline
size_t store( std::ostream& _os, const T& _v, bool _swap=false)
{ return binary< T >::store( _os, _v, _swap ); }

template <typename T> inline
size_t restore( std::istream& _is, T& _v, bool _swap=false)
{ return binary< T >::restore( _is, _v, _swap ); }

//@}


//=============================================================================
} // namespace IO
} // namespace OpenMesh
//=============================================================================
#endif // OPENMESH_MESHREADER_HH defined
//=============================================================================
