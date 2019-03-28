#include <OpenMesh/Core/IO/reader/PLYReader.hh>

namespace OpenMesh  {
namespace IO  {

_PLYReader_  __PLYReaderInstance;
_PLYReader_& PLYReader() { return __PLYReaderInstance; }


_PLYReader_::_PLYReader_()  {
	IOManager().register_module(this); 
}

bool _PLYReader_::read(const std::string& _filename, 
	 BaseImporter& _bi, 
	 Options& _opt)  {

	FILE* in = fopen(_filename.c_str(), "r");
	if (!in)
	{
		omerr() << "[PLYReader] : cannot not open file " 
			<< _filename 
			<< std::endl;
		return false;
	}

	{
#if defined(WIN32)
		std::string::size_type dot = _filename.find_last_of("\\/");
#else
		std::string::size_type dot = _filename.rfind("/");
#endif
		path_ = (dot == std::string::npos) 
			? "./" 
			: std::string(_filename.substr(0,dot+1));
	}

	bool ok = read(in, _bi, _opt);

	fclose(in);
	return ok;
}

bool _PLYReader_::read(FILE* _in, BaseImporter& _bi, Options& _opt)  {
	PlyFile* plyFile = read_ply(_in);

	char x_char[2], y_char[2], z_char[2], nx_char[3], ny_char[3], nz_char[3];
	x_char[0] = 'x'; x_char[1] = NULL;
	y_char[0] = 'y'; y_char[1] = NULL;
	z_char[0] = 'z'; z_char[1] = NULL;
	nx_char[0] = 'n'; nx_char[1] = 'x'; nx_char[2] = NULL;
	ny_char[0] = 'n'; ny_char[1] = 'y'; ny_char[2] = NULL;
	nz_char[0] = 'n'; nz_char[1] = 'z'; nz_char[2] = NULL;

	PlyProperty vert_props[] = {
	  {x_char, Float32, Float32, offsetof(Vertex,x), 0, 0, 0, 0},
	  {y_char, Float32, Float32, offsetof(Vertex,y), 0, 0, 0, 0},
	  {z_char, Float32, Float32, offsetof(Vertex,z), 0, 0, 0, 0},
	  {nx_char, Float32, Float32, offsetof(Vertex,nx), 0, 0, 0, 0},
	  {ny_char, Float32, Float32, offsetof(Vertex,ny), 0, 0, 0, 0},
	  {nz_char, Float32, Float32, offsetof(Vertex,nz), 0, 0, 0, 0},
	};

	char vertex_indices[15];
	vertex_indices[0] = 'v';
	vertex_indices[1] = 'e';
	vertex_indices[2] = 'r';
	vertex_indices[3] = 't';
	vertex_indices[4] = 'e';
	vertex_indices[5] = 'x';
	vertex_indices[6] = '_';
	vertex_indices[7] = 'i';
	vertex_indices[8] = 'n';
	vertex_indices[9] = 'd';
	vertex_indices[10] = 'i';
	vertex_indices[11] = 'c';
	vertex_indices[12] = 'e';
	vertex_indices[13] = 's';
	vertex_indices[14] = NULL;

	PlyProperty face_props[] = {
	  {vertex_indices, Int32, Int32, offsetof(Face,verts),
		1, Uint8, Uint8, offsetof(Face,nverts)},
	};

	int nverts, nfaces;
	PLYVertex **vlist;
	PLYFace **flist;
	PlyOtherProp *vert_other,*face_other;
	int has_nx = 0, has_ny = 0, has_nz = 0;

	for(int i = 0; i < plyFile->num_elem_types; i++)  {
		PlyElement* element = plyFile->elems[i];
		int elem_count;
		const char* element_name = setup_element_read_ply (plyFile, i, &elem_count);
		int prop_count = element->nprops;

		if(equal_strings(element_name, "vertex"))  {

			/* create a vertex list to hold all the vertices */
			vlist = (Vertex **) malloc (sizeof (Vertex *) * elem_count);
			nverts = elem_count;

			/* set up for getting vertex elements */
			/* (we want x,y,z) */

			setup_property_ply (plyFile, &vert_props[0]);
			setup_property_ply (plyFile, &vert_props[1]);
			setup_property_ply (plyFile, &vert_props[2]);

			/* we also want normal information if it is there (nx,ny,nz) */

			for (int j = 0; j < prop_count; j++) {
				PlyProperty *prop;
				prop = element->props[j];
				if (equal_strings ("nx", prop->name)) {
				  setup_property_ply (plyFile, &vert_props[3]);
				  has_nx = 1;
				}
				if (equal_strings ("ny", prop->name)) {
				  setup_property_ply (plyFile, &vert_props[4]);
				  has_ny = 1;
				}
				if (equal_strings ("nz", prop->name)) {
				  setup_property_ply (plyFile, &vert_props[5]);
				  has_nz = 1;
				}
			}

			/* also grab anything else that we don't need to know about */

			vert_other = get_other_properties_ply (plyFile, 
							  offsetof(Vertex,other_props));

			/* grab the vertex elements and store them in our list */

			for (int j = 0; j < elem_count; j++) {
			  vlist[j] = (Vertex *) malloc (sizeof (Vertex));
			  get_element_ply (plyFile, (void *) vlist[j]);
			}
		}

		else if(equal_strings(element_name, "face"))  {

			/* create a list to hold all the face elements */
			flist = (Face **) malloc (sizeof (Face *) * elem_count);
			nfaces = elem_count;

			/* set up for getting face elements */
			/* (all we need are vertex indices) */

			setup_property_ply (plyFile, &face_props[0]);
			face_other = get_other_properties_ply (plyFile, 
							  offsetof(Face,other_props));

			/* grab all the face elements and place them in our list */

			for (int j = 0; j < elem_count; j++) {
			  flist[j] = (Face *) malloc (sizeof (Face));
			  get_element_ply (plyFile, (void *) flist[j]);
			}
		}

		else  /* all non-vertex and non-face elements are grabbed here */
			get_other_element_ply (plyFile);
	}

	bool hasNormals = (has_nx == 1 && has_ny == 1 && has_nz == 1);

	// load vertices
	for(int v = 0; v < nverts; v++)  {
		Vertex* vertex = vlist[v];
	   _bi.add_vertex(OpenMesh::Vec3f(vertex->x,vertex->y,vertex->z));
	}

	BaseImporter::VHandles vhandles;
	for(int f = 0; f < nfaces; f++)  {
		vhandles.clear();

		Face* face = flist[f];
		unsigned char num_verts = face->nverts;
		int *vert_list = face->verts;

		for(int v = 0; v < num_verts; v++)  {
			Vertex* next_vert = vlist[vert_list[v]];
		   vhandles.push_back(VertexHandle(vert_list[v]));
			if(hasNormals)
				_bi.set_normal(vhandles.back(), OpenMesh::Vec3f(next_vert->nx,next_vert->ny,next_vert->nz));
		}

      FaceHandle fh = _bi.add_face(vhandles);
	}

	free(vert_other);
	free(face_other);
	for(int v = 0; v < nverts; v++)  {
		free(vlist[v]);
	}
	free(vlist);
	for(int f = 0; f < nfaces; f++)  {
		free(flist[f]->verts);
		free(flist[f]);
	}
	free(flist);

	//printf("closing ply...\n");
	free_ply(plyFile);

	return true;
}

} // end namespace IO
} // end namespace OpenMesh
