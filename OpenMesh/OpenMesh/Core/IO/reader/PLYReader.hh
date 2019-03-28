#ifndef PLYREADER_HH
#define PLYREADER_HH

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <stdio.h>
#include <stdlib.h>

#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Utils/SingletonT.hh>
#include <OpenMesh/Core/IO/IOManager.hh>
#include <OpenMesh/Core/IO/importer/BaseImporter.hh>
#include <OpenMesh/Core/IO/reader/BaseReader.hh>

#include <OpenMesh/Core/IO/reader/ply.hh>

namespace OpenMesh  {
namespace IO  {

class _PLYReader_ : public BaseReader  {
	public:

		_PLYReader_();

		virtual ~_PLYReader_() { }

		std::string get_description() const { return "PLY"; }
		std::string get_extensions()  const { return "ply"; }

		bool read(const std::string& _filename, 
			 BaseImporter& _bi, 
			 Options& _opt);

	private:

		bool read(FILE* _in, BaseImporter& _bi, Options& _opt);

		std::string path_;

		typedef struct Vertex {
		  float x,y,z;
		  float nx,ny,nz;
		  void *other_props;       /* other properties */
		} PLYVertex;

		typedef struct Face {
		  unsigned char nverts;    /* number of vertex indices in list */
		  int *verts;              /* vertex index list */
		  void *other_props;       /* other properties */
		} PLYFace;

};

extern _PLYReader_  __PLYReaderInstance;
_PLYReader_& PLYReader();

} // end namespace IO
} // end namespace OpenMesh

#endif
