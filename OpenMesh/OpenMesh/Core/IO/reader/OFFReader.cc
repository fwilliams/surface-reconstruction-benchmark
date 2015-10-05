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


#define LINE_LEN 4096


//== INCLUDES =================================================================

// OpenMesh
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/System/omstream.hh>
#include <OpenMesh/Core/IO/reader/OFFReader.hh>
#include <OpenMesh/Core/IO/importer/BaseImporter.hh>
#include <OpenMesh/Core/IO/IOManager.hh>
// #include <OpenMesh/Core/IO/BinaryHelper.hh>

//STL
#include <fstream>
#include <memory>

#ifndef WIN32
#include <string.h>
#endif

//=== NAMESPACES ==============================================================


namespace OpenMesh {
namespace IO {


//=============================================================================

//=== INSTANCIATE =============================================================


_OFFReader_  __OFFReaderInstance;
_OFFReader_&  OFFReader() { return __OFFReaderInstance; }


//=== IMPLEMENTATION ==========================================================



_OFFReader_::_OFFReader_() 
{ 
  IOManager().register_module(this); 
}


//-----------------------------------------------------------------------------


bool 
_OFFReader_::read(const std::string& _filename, BaseImporter& _bi, 
                  Options& _opt)
{
  std::fstream in(_filename.c_str(), (options_.is_binary ? std::ios_base::binary | std::ios_base::in
                                                         : std::ios_base::in) );

  if (!in)
  {
    omerr() << "[OFFReader] : cannot not open file " 
	  << _filename
	  << std::endl;
    return false;
  }


  bool result = read(in, _bi, _opt);


  in.close();
  return result;
}


//-----------------------------------------------------------------------------


bool 
_OFFReader_::read(std::fstream& _in, BaseImporter& _bi, Options& _opt ) const
{
   // filter relevant options for reading
   bool swap = _opt.check( Options::Swap );
  

   // build options to be returned
   _opt.clear();
   if (options_.vertex_has_normal)   _opt += Options::VertexNormal;
   if (options_.vertex_has_texcoord) _opt += Options::VertexTexCoord;
   if (options_.vertex_has_color)    _opt += Options::VertexColor;
   if (options_.is_binary)           _opt += Options::Binary;     


    return (options_.is_binary ?
 	   read_binary(_in, _bi, swap) :
	   read_ascii(_in, _bi));

}



//-----------------------------------------------------------------------------

bool 
_OFFReader_::read_ascii(std::fstream& _in, BaseImporter& _bi) const
{


omlog() << "[OFFReader] : read ascii file\n";
   
  unsigned int            i, j, k, l, idx;
  unsigned int            nV, nF, dummy;
  OpenMesh::Vec3f         v, n;
  OpenMesh::Vec2f         t;
  OpenMesh::Vec3i         c;
  BaseImporter::VHandles  vhandles;
  VertexHandle            vh;

  // read header line
  std::string header;
  std::getline(_in,header);

  // + optional space dimension
  if (options_.vertex_dim)
  {
    _in >> i;
    if (i != 3) // cannot process more 3 dimensions!
      return false;
  }


  // + #Vertice, #Faces, #Edges
  _in >> nV;
  _in >> nF;
  _in >> dummy;

  _bi.reserve(nV, 3*nV, nF);

  // read vertices: coord [hcoord] [normal] [color] [texcoord]
  for (i=0; i<nV && !_in.eof(); ++i)
  {
    // Always read Vertex
    _in >> v[0]; _in >> v[1]; _in >> v[2];
    
    vh = _bi.add_vertex(v);
    
    
    if ( options_.vertex_has_normal ) {
      _in >> n[0]; _in >> n[1]; _in >> n[2];
      _bi.set_normal(vh, n);
    }
    
    if ( options_.vertex_has_color ) {
      _in >> c[0]; _in >> c[1]; _in >> c[2];
      _bi.set_color( vh, Vec3uc( c ) );
    }
    
    if ( options_.vertex_has_texcoord) {
      _in >> t[0]; _in >> t[1];
      _bi.set_texcoord(vh, t);
    }
  }
  
  
  // faces
  // #N <v1> <v2> .. <v(n-1)> [color spec]
  // So far color spec is unsupported!
  for (i=0; i<nF; ++i)
  {
    _in >> nV;


    if (nV == 3)
    {
      vhandles.resize(3);
      _in >> j;
      _in >> k;
      _in >> l;

      vhandles[0] = VertexHandle(j);
      vhandles[1] = VertexHandle(k);
      vhandles[2] = VertexHandle(l);
    }
    else 
    {
      vhandles.clear();
      for (j=0; j<nV; ++j)
      {
         _in >> idx;
         vhandles.push_back(VertexHandle(idx));
      }
    }
    
    _bi.add_face(vhandles);
  }

  // File was successfully parsed.
  return true;
}


//-----------------------------------------------------------------------------

unsigned int _OFFReader_::getInt(std::fstream& _in) const {
  unsigned int value;
  _in.read(reinterpret_cast < char * > (&value), sizeof(value));
  return value;
}

float _OFFReader_::getFloat(std::fstream& _in) const{
  float value;
  _in.read(reinterpret_cast < char * > (&value), sizeof(value));
  return value;
}

bool 
_OFFReader_::read_binary(std::fstream& _in, BaseImporter& _bi, bool _swap) const
{
  omlog() << "[OFFReader] : read ascii file\n";
   
  unsigned int            i, j, k, l, idx;
  unsigned int            nV, nF, dummy;
  OpenMesh::Vec3f         v, n;
  OpenMesh::Vec3i         c;
  OpenMesh::Vec2f         t;
  BaseImporter::VHandles  vhandles;
  VertexHandle            vh;

  // read header line
  std::string header;
  std::getline(_in,header);

  // + optional space dimension
  if (options_.vertex_dim)
  {
    _in >> i;
    if (i != 3) // cannot process more 3 dimensions!
      return false;
  }


  // + #Vertice, #Faces, #Edges
  nV = getInt(_in);
  nF = getInt(_in);
  dummy = getInt(_in);

  _bi.reserve(nV, 3*nV, nF);

  // read vertices: coord [hcoord] [normal] [color] [texcoord] 
  for (i=0; i<nV && !_in.eof(); ++i)
  {
    // Always read Vertex
    v[0] = getFloat(_in);
    v[1] = getFloat(_in);
    v[2] = getFloat(_in);
    
    vh = _bi.add_vertex(v);
    
    if ( options_.vertex_has_normal ) {
      n[0] = getFloat(_in);
      n[1] = getFloat(_in);
      n[2] = getFloat(_in);
      
      _bi.set_normal(vh, n);
    }
    
    if ( options_.vertex_has_color ) {
      c [0] = getInt(_in); 
      c [1] = getInt(_in);
      c [2] = getInt(_in);
      _bi.set_color( vh, Vec3uc( c ) );
    }
    
    if ( options_.vertex_has_texcoord) {
      t[0] = getFloat(_in);
      t[1] = getFloat(_in);
      
      _bi.set_texcoord(vh, t);
    }
  }
  
  // faces
  // #N <v1> <v2> .. <v(n-1)> [color spec]
  // So far color spec is unsupported!
  for (i=0; i<nF; ++i)
  {
    nV = getInt(_in);


    if (nV == 3)
    {
      vhandles.resize(3);
      j = getInt(_in);
      k = getInt(_in);
      l = getInt(_in);

      vhandles[0] = VertexHandle(j);
      vhandles[1] = VertexHandle(k);
      vhandles[2] = VertexHandle(l);
    }
    else 
    {
      vhandles.clear();
      for (j=0; j<nV; ++j)
      {
         idx = getInt(_in);
         vhandles.push_back(VertexHandle(idx));
      }
    }
    
    _bi.add_face(vhandles);
  }

  // File was successfully parsed.
  return true;

}


//-----------------------------------------------------------------------------


bool _OFFReader_::can_u_read(const std::string& _filename) const
{
  // !!! Assuming BaseReader::can_u_parse( std::string& )
  // does not call BaseReader::read_magic()!!!
  if (BaseReader::can_u_read(_filename))
  {
    std::ifstream ifs(_filename.c_str());
    if (ifs.is_open() && can_u_read(ifs))
    {
      ifs.close();
      return true;
    }      
  }
  return false;
}


//-----------------------------------------------------------------------------

  
bool _OFFReader_::can_u_read(std::istream& _is) const
{
  // read 1st line
  char line[LINE_LEN], *p;
  _is.getline(line, LINE_LEN);
  p = line;

  int remainingChars = _is.gcount();

  

  // check header: [ST][C][N][4][n]OFF BINARY
  memset( &options_, 0, sizeof(options_) );

  if ( ( remainingChars > 1 ) && ( p[0] == 'S' && p[1] == 'T') )
  { options_.vertex_has_texcoord = true; p += 2; remainingChars -= 2; }

  if ( ( remainingChars > 0 ) && ( p[0] == 'C') )
  { options_.vertex_has_color = true; ++p; --remainingChars; }

  if ( ( remainingChars > 0 ) && ( p[0] == 'N') )
  { options_.vertex_has_normal = true; ++p; --remainingChars; }

  if ( ( remainingChars > 0 ) && (p[0] == '4' ) )
  { options_.vertex_has_hcoord = true; ++p; --remainingChars; }

  if ( ( remainingChars > 0 ) && ( p[0] == 'n') )
  { options_.vertex_dim = true; ++p; --remainingChars; }
  
  if ( ( remainingChars < 3 ) || (!(p[0] == 'O' && p[1] == 'F' && p[2] == 'F')  ) )
    return false;

  p += 4;
  remainingChars -= 4;

  if ( ( remainingChars >= 6 ) && ( strncmp(p, "BINARY", 6) == 0 ) )
    options_.is_binary = true;

  // currently these options are not supported!!!
  if (options_.vertex_dim  || options_.vertex_has_hcoord)
    return false;
   

  return true;
}


//=============================================================================
} // namespace IO
} // namespace OpenMesh
//=============================================================================
