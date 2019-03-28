/*===========================================================================*\
 *                                                                           *
 *                               OpenMesh                                    *
 *        Copyright (C) 2003 by Computer Graphics Group, RWTH Aachen         *
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


//=============================================================================
//
//  Implements a reader module for OFF files
//
//=============================================================================


#ifndef __OMREADER_HH__
#define __OMREADER_HH__


//=== INCLUDES ================================================================

// OpenMesh
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Utils/SingletonT.hh>
#include <OpenMesh/Core/IO/OMFormat.hh>
#include <OpenMesh/Core/IO/IOManager.hh>
#include <OpenMesh/Core/IO/importer/BaseImporter.hh>
#include <OpenMesh/Core/IO/reader/BaseReader.hh>

// STD C++
#include <iostream>
#include <string>


//== NAMESPACES ===============================================================


namespace OpenMesh {
namespace IO {


//== IMPLEMENTATION ===========================================================


/** 
    Implementation of the OM format reader. This class is singleton'ed by 
    SingletonT to OMReader.
*/
class _OMReader_ : public BaseReader
{
public:

  _OMReader_();
  virtual ~_OMReader_() { }

  std::string get_description() const { return "OpenMesh File Format"; }
  std::string get_extensions()  const { return "om"; }
  std::string get_magic()       const { return "OM"; }
   
  bool read(const std::string& _filename, 
	    BaseImporter& _bi, 
	    Options& _opt );

  virtual bool can_u_read(const std::string& _filename) const;
  virtual bool can_u_read(std::istream& _is) const;

  
private:

  bool supports( const OMFormat::uint8 version ) const;

  bool read(std::istream& _is, BaseImporter& _bi, Options& _opt ) const;
  bool read_ascii(std::istream& _is, BaseImporter& _bi, Options& _opt) const;
  bool read_binary(std::istream& _is, BaseImporter& _bi, Options& _opt) const;

  typedef OMFormat::Header              Header;
  typedef OMFormat::Chunk::Header       ChunkHeader;
  typedef OMFormat::Chunk::PropertyName PropertyName;

  // initialized/updated by read_binary*/read_ascii*
  mutable size_t       bytes_;
  mutable Header       header_;
  mutable ChunkHeader  chunk_header_;
  mutable PropertyName property_name_;

  bool read_binary_vertex_chunk(   std::istream      &_is, 
				   BaseImporter      &_bi, 
				   Options           &_opt,
				   bool              _swap) const;

  bool read_binary_face_chunk(     std::istream      &_is, 
			           BaseImporter      &_bi, 
			           Options           &_opt,
				   bool              _swap) const;

  bool read_binary_edge_chunk(     std::istream      &_is, 
			           BaseImporter      &_bi, 
			           Options           &_opt,
				   bool              _swap) const;

  bool read_binary_halfedge_chunk( std::istream      &_is, 
				   BaseImporter      &_bi, 
				   Options           &_opt,
				   bool              _swap) const;

  bool read_binary_mesh_chunk(     std::istream      &_is, 
				   BaseImporter      &_bi, 
				   Options           &_opt,
				   bool              _swap) const;

  size_t restore_binary_custom_data( std::istream& _is, 
				     BaseProperty* _bp,
				     size_t _n_elem, 
				     bool _swap) const;

};


//== TYPE DEFINITION ==========================================================


/// Declare the single entity of the OM reader.
extern _OMReader_  __OMReaderInstance;
_OMReader_&  OMReader();


//=============================================================================
} // namespace IO
} // namespace OpenMesh
//=============================================================================
#endif
//=============================================================================
