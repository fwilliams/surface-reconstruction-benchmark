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


//== INCLUDES =================================================================


#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/System/omstream.hh>
#include <OpenMesh/Core/Utils/Endian.hh>
#include <OpenMesh/Core/IO/IOManager.hh>
#include <OpenMesh/Core/IO/BinaryHelper.hh>
#include <OpenMesh/Core/IO/writer/OFFWriter.hh>


//=== NAMESPACES ==============================================================


namespace OpenMesh {
namespace IO {


//=== INSTANCIATE =============================================================


// register the OFFLoader singleton with MeshLoader
_OFFWriter_  __OFFWriterInstance;
_OFFWriter_& OFFWriter() { return __OFFWriterInstance; }


//=== IMPLEMENTATION ==========================================================


_OFFWriter_::_OFFWriter_() { IOManager().register_module(this); }


//-----------------------------------------------------------------------------


bool
_OFFWriter_::
write(const std::string& _filename, BaseExporter& _be, Options _opt) const
{
  // check exporter features
  if ( !check( _be, _opt ) )
    return false;


  // check writer features
  if ( _opt.check(Options::FaceNormal)  || // not supported by format
       _opt.check(Options::FaceColor)  )   // not supported by module
    return false;
     
  // open file
  std::fstream out(_filename.c_str(), (_opt.check(Options::Binary) ? std::ios_base::binary | std::ios_base::out
                                                         : std::ios_base::out) );
  if (!out)
  {
    omerr() << "[OFFWriter] : cannot open file "
	  << _filename
	  << std::endl;
    return false;
  }


  // write header line
  if (_opt.check(Options::VertexTexCoord)) out << "ST";
  if (_opt.check(Options::VertexColor))    out << "C";
  if (_opt.check(Options::VertexNormal))   out << "N";
  out << "OFF";
  if (_opt.check(Options::Binary)) out << " BINARY";
  out << "\n";


  // write to file
  bool result = (_opt.check(Options::Binary) ?
		 write_binary(out, _be, _opt) :
		 write_ascii(out, _be, _opt));


  // return result
  out.close();
  return result;
}


//-----------------------------------------------------------------------------


bool
_OFFWriter_::
write_ascii(std::fstream& _out, BaseExporter& _be, Options _opt) const
{
  omlog() << "[OFFWriter] : write ascii file\n";
  _out.precision(10);


  unsigned int i, j, nV, nF;
  Vec3f v, n;
  Vec2f t;
  OpenMesh::Vec3i c;
  VertexHandle vh;
  std::vector<VertexHandle> vhandles;


  // #vertices, #faces
  _out << _be.n_vertices() << " ";
  _out << _be.n_faces() << " ";
  _out << 0 << "\n";


  // vertex data (point, normals, texcoords)
  for (i=0, nV=_be.n_vertices(); i<nV; ++i)
  {
    vh = VertexHandle(i);
    v  = _be.point(vh);
    
    _out << v[0] << " " << v[1] << " " << v[2];
    
    if ( _opt.check(Options::VertexNormal) ) {
      n  = _be.normal(vh);
      _out << " " << n[0] << " " << n[1] << " " << n[2];
    }
    
    if ( _opt.check(Options::VertexColor) ) {
      c  = _be.color(vh);
      _out << " " << c[0] << " " << c[1] << " " << c[2];
    }
    
    if (_opt.check(Options::VertexTexCoord) ) {
      t  = _be.texcoord(vh);
      _out << " " << t[0] << " " << t[1];
    }
    
    _out << "\n";
    
  }
  
  // faces (indices starting at 0)
  if (_be.is_triangle_mesh())
  {
    for (i=0, nF=_be.n_faces(); i<nF; ++i)
    {
      _be.get_vhandles(FaceHandle(i), vhandles);
      _out << 3 << " "; 
      _out << vhandles[0].idx()  << " ";
      _out << vhandles[1].idx()  << " "; 
      _out << vhandles[2].idx()  << "\n";
    }
  }
  else
  {
    for (i=0, nF=_be.n_faces(); i<nF; ++i)
    {
      nV = _be.get_vhandles(FaceHandle(i), vhandles);
      _out << nV << " ";
      for (j=0; j<vhandles.size(); ++j)
	       _out << vhandles[j].idx() << " ";
      _out << "\n";
    }
  }


  return true;
}


//-----------------------------------------------------------------------------

void _OFFWriter_::writeInt(std::fstream& _out, unsigned int value) const {
  _out.write(reinterpret_cast<char *>(&value),sizeof(value));
}

void _OFFWriter_::writeFloat(std::fstream& _out, float value) const {
  _out.write(reinterpret_cast<char *>(&value),sizeof(value));
}

bool 
_OFFWriter_::
write_binary(std::fstream& _out, BaseExporter& _be, Options _opt) const
{
  omlog() << "[OFFWriter] : write ascii file\n";


  unsigned int i, j, nV, nF;
  Vec3f v, n;
  Vec2f t;
  OpenMesh::Vec3i c;
  VertexHandle vh;
  std::vector<VertexHandle> vhandles;


  // #vertices, #faces
  writeInt(_out, _be.n_vertices() );
  writeInt(_out, _be.n_faces() );
  writeInt(_out, 0 );

  
  // vertex data (point, normals, texcoords)
  for (i=0, nV=_be.n_vertices(); i<nV; ++i)
  {
    vh = VertexHandle(i);
    v  = _be.point(vh);
    
    writeFloat(_out, v[0]);
    writeFloat(_out, v[1]);
    writeFloat(_out, v[2]);
    
    if ( _opt.check(Options::VertexNormal) ) {
      n  = _be.normal(vh);
      writeFloat(_out, n[0]);
      writeFloat(_out, n[1]);
      writeFloat(_out, n[2]);
    }
    
    if ( _opt.check(Options::VertexColor) ) {
      c  = _be.color(vh);
      writeInt(_out, c[0]);
      writeInt(_out, c[1]);
      writeInt(_out, c[2]);
    }
    
    if (_opt.check(Options::VertexTexCoord) ) {
      t  = _be.texcoord(vh);
      writeFloat(_out, t[0]);
      writeFloat(_out, t[1]);
    }
    
  }

  // faces (indices starting at 0)
  if (_be.is_triangle_mesh())
  {
    for (i=0, nF=_be.n_faces(); i<nF; ++i)
    {
      _be.get_vhandles(FaceHandle(i), vhandles);
      writeInt(_out, 3);
      writeInt(_out, vhandles[0].idx());
      writeInt(_out, vhandles[1].idx());
      writeInt(_out, vhandles[2].idx());
    }
  }
  else
  {
    for (i=0, nF=_be.n_faces(); i<nF; ++i)
    {
      nV = _be.get_vhandles(FaceHandle(i), vhandles);
      writeInt(_out, nV);
      for (j=0; j<vhandles.size(); ++j)
        writeInt(_out, vhandles[j].idx() );
    }
  }

  return true;
}

// ----------------------------------------------------------------------------


size_t
_OFFWriter_::
binary_size(BaseExporter& _be, Options _opt) const
{
  size_t header(0);
  size_t data(0);
  size_t _3longs(3*sizeof(long));
  size_t _3floats(3*sizeof(float));
  
  if ( !_opt.check(Options::Binary) )
    return 0;
  else
  {
    header += 11;                             // 'OFF BINARY\n'
    header += _3longs;                        // #V #F #E
    data   += _be.n_vertices() * _3floats;    // vertex data
  }

  if ( _opt.check(Options::VertexNormal) && _be.has_vertex_normals() )
  {
    header += 1; // N
    data   += _be.n_vertices() * _3floats;
  }
  
  if ( _opt.check(Options::VertexColor) && _be.has_vertex_colors() )
  {
    header += 1; // C
    data   += _be.n_vertices() * _3floats;
  }

  if ( _opt.check(Options::VertexTexCoord) && _be.has_vertex_texcoords() )
  {
    size_t _2floats(2*sizeof(float));
    header += 2; // ST
    data   += _be.n_vertices() * _2floats;
  }
   

  // topology
  if (_be.is_triangle_mesh())
  {
    size_t _4ui(4*sizeof(unsigned int));                  
    data += _be.n_faces() * _4ui;
  }
  else
  {
    unsigned int i, nV, nF;
    std::vector<VertexHandle> vhandles;

    for (i=0, nF=_be.n_faces(); i<nF; ++i)
    {
      nV = _be.get_vhandles(FaceHandle(i), vhandles);
      data += nV * sizeof(unsigned int);
    }
  }
   
  return header+data;
}


//=============================================================================
} // namespace IO
} // namespace OpenMesh
//=============================================================================
