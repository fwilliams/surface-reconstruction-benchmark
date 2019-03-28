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


#ifndef __OFFREADER_HH__
#define __OFFREADER_HH__


//=== INCLUDES ================================================================


#include <iostream>
#include <string>
#include <stdio.h>
#include <fstream>

#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Utils/SingletonT.hh>
#include <OpenMesh/Core/IO/reader/BaseReader.hh>

#ifndef WIN32
#include <string.h>
#endif

//== NAMESPACES ===============================================================


namespace OpenMesh {
namespace IO {


//== FORWARDS =================================================================


class BaseImporter;


//== IMPLEMENTATION ===========================================================


/** 
    Implementation of the OFF format reader. This class is singleton'ed by 
    SingletonT to OFFReader.
*/
class _OFFReader_ : public BaseReader
{
public:

  _OFFReader_();

  std::string get_description() const { return "Object File Format"; }
  std::string get_extensions()  const { return "off"; }
  std::string get_magic()       const { return "OFF"; }
   
  bool read(const std::string& _filename, 
	    BaseImporter& _bi, 
	    Options& _opt);

  bool can_u_read(const std::string& _filename) const;



private:

  bool can_u_read(std::istream& _is) const;

  bool read(std::fstream& _in, BaseImporter& _bi, Options& _opt ) const;
  bool read_ascii(std::fstream& _in, BaseImporter& _bi) const;
  bool read_binary(std::fstream& _in, BaseImporter& _bi, bool swap) const;
   
  unsigned int getInt(std::fstream& _in) const;
  float getFloat(std::fstream& _in) const;
  
  /// will be initialized, when calling can_u_parse()
  mutable struct {
    unsigned vertex_dim:8; // allow up to 256 dimensional vectors
    unsigned vertex_has_normal:1;
    unsigned vertex_has_color:1;
    unsigned vertex_has_texcoord:1;
    unsigned vertex_has_hcoord:1;
    unsigned is_binary:1;
  } options_;
};


//== TYPE DEFINITION ==========================================================


/// Declare the single entity of the OFF reader
extern _OFFReader_  __OFFReaderInstance;
_OFFReader_&  OFFReader();


//=============================================================================
} // namespace IO
} // namespace OpenMesh
//=============================================================================
#endif
//=============================================================================
