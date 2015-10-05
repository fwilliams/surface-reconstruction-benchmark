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
//  Implements a writer module for OFF files
//
//=============================================================================


#ifndef __OFFWRITER_HH__
#define __OFFWRITER_HH__


//=== INCLUDES ================================================================


#include <stdio.h>
#include <string>
#include <fstream>

#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Utils/SingletonT.hh>
#include <OpenMesh/Core/IO/exporter/BaseExporter.hh>
#include <OpenMesh/Core/IO/writer/BaseWriter.hh>


//== NAMESPACES ===============================================================


namespace OpenMesh {
namespace IO {


//=== IMPLEMENTATION ==========================================================


/** 
    Implementation of the OFF format writer. This class is singleton'ed by 
    SingletonT to OFFWriter.
*/
class _OFFWriter_ : public BaseWriter
{
public:

  _OFFWriter_();

  std::string get_description() const { return "no description"; }
  std::string get_extensions() const  { return "off"; }

  bool write(const std::string&, BaseExporter&, Options) const;

  size_t binary_size(BaseExporter& _be, Options _opt) const;


protected:
  void writeInt(std::fstream& _out, unsigned int value) const;
  void writeFloat(std::fstream& _out, float value) const;

  bool write_ascii(std::fstream& _in, BaseExporter&, Options) const;
  bool write_binary(std::fstream& _in, BaseExporter&, Options) const;
};


//== TYPE DEFINITION ==========================================================


/// Declare the single entity of the OFF writer.
extern _OFFWriter_  __OFFWriterInstance;
_OFFWriter_& OFFWriter();


//=============================================================================
} // namespace IO
} // namespace OpenMesh
//=============================================================================
#endif
//=============================================================================
