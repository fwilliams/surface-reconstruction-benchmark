/*
Copyright (c) 2010, Matt Berger
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list
  of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list
  of conditions and the following disclaimer in the documentation and/or other
  materials provided with the distribution.
- Neither the name of the University of Utah nor the names of its contributors may be used
  to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "PNGImage.h"

PNGImage::PNGImage(string _filename, int _resx, int _resy)  {
	filename = _filename;
	res_x = _resx;
	res_y = _resy;

	if(!this->init_file())
		cerr << "bad initialization!" << endl;

	img_pointer = (png_bytep*) malloc(sizeof(png_bytep) * res_y);
	for (int y=0; y < res_y; y++)
		img_pointer[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr,info_ptr));
}

PNGImage::PNGImage(string _filein)  {
	filename = _filein;
	if(!this->read_png_file())
		cerr << "bad png file!" << endl;
}

PNGImage::~PNGImage()  {
	/* cleanup heap allocation */
	for (int y=0; y<res_y; y++)
		free(img_pointer[y]);
	free(img_pointer);
	fclose(fp);
}

bool PNGImage::init_file()  {
	/* create file */
	fp = fopen(filename.c_str(), "wb");
	if (!fp)  {
		cerr << "bad 1" << endl;
		return false;
	}

	/* initialize stuff */
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

	if (!png_ptr)  {
		cerr << "bad png ptr" << endl;
		return false;
	}

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)  {
		cerr << "bad info ptr" << endl;
		return false;
	}

	if (setjmp(png_jmpbuf(png_ptr)))  {
		cerr << "bad setjmp" << endl;
		return false;
	}

	png_init_io(png_ptr, fp);

	/* write header */
	if (setjmp(png_jmpbuf(png_ptr)))  {
		cerr << "bad write header" << endl;
		return false;
	}

	png_set_IHDR(png_ptr, info_ptr, res_x, res_y,
					8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
					PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
	png_write_info(png_ptr, info_ptr);

	return true;
}

bool PNGImage::read_png_file()  {
	png_byte header[8];    // 8 is the maximum size that can be checked

	/* open file and test for it being a png */
	fp = fopen(filename.c_str(), "rb");
	if (!fp)
		return false;
	fread(header, 1, 8, fp);
	if (png_sig_cmp(header, 0, 8))
		return false;

	/* initialize stuff */
	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

	if (!png_ptr)
		return false;

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
		return false;

	if (setjmp(png_jmpbuf(png_ptr)))
		return false;

	png_init_io(png_ptr, fp);
	png_set_sig_bytes(png_ptr, 8);

	png_read_info(png_ptr, info_ptr);

	res_x = png_get_image_width(png_ptr, info_ptr);
	res_y = png_get_image_height(png_ptr, info_ptr);
	png_byte color_type = png_get_color_type(png_ptr, info_ptr);
	png_byte bit_depth = png_get_bit_depth(png_ptr, info_ptr);

	int number_of_passes = png_set_interlace_handling(png_ptr);
	png_read_update_info(png_ptr, info_ptr);

	/* read file */
	if (setjmp(png_jmpbuf(png_ptr)))
		return false;

	img_pointer = (png_bytep*) malloc(sizeof(png_bytep) * res_y);
	for (int y=0; y < res_y; y++)
		img_pointer[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr,info_ptr));

	png_read_image(png_ptr, img_pointer);

	switch(png_get_color_type(png_ptr, info_ptr))  {
		case PNG_COLOR_TYPE_RGB:
			cout << "RGB" << endl;
			break;
		case PNG_COLOR_TYPE_RGBA:
			cout << "RGBA" << endl;
			break;
	}

	return true;
}

bool PNGImage::write_to_file()  {
	/* write bytes */
	if (setjmp(png_jmpbuf(png_ptr)))  {
		cerr << "bad write bytes" << endl;
		return false;
	}

	png_write_image(png_ptr, img_pointer);

	/* end write */
	if (setjmp(png_jmpbuf(png_ptr)))  {
		cerr << "bad write end" << endl;
		return false;
	}

	png_write_end(png_ptr, NULL);

	return true;
}

void PNGImage::get_png_pixel(int _x, int _y, png_byte& pix_r, png_byte& pix_g, png_byte& pix_b)  {
	png_byte* row = img_pointer[_y];
	png_byte* pixel_ptr = &(row[_x*3]);
	pix_r = pixel_ptr[0];
	pix_g = pixel_ptr[1];
	pix_b = pixel_ptr[2];
}

void PNGImage::set_png_pixel(int _x, int _y, double _r, double _g, double _b)  {
	png_byte r_byte = (png_byte)(255.0*_r);
	png_byte g_byte = (png_byte)(255.0*_g);
	png_byte b_byte = (png_byte)(255.0*_b);

	png_byte* row = img_pointer[_y];
	png_byte* pixel_ptr = &(row[_x*3]);
	pixel_ptr[0] = r_byte;
	pixel_ptr[1] = g_byte;
	pixel_ptr[2] = b_byte;
}
