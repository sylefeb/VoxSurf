/* --------------------------------------------------------------------
Author: Sylvain Lefebvre    sylvain.lefebvre@sophia.inria.fr

                  Simple Library for Graphics (LibSL)

This software is a computer program whose purpose is to offer a set of
tools to simplify programming real-time computer graphics applications
under OpenGL and DirectX.

This software is governed by the CeCILL-C license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL-C
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-C license and that you accept its terms.
-------------------------------------------------------------------- */

#include <stdio.h>

extern "C" {
#include "jpeglib.h"
}
//#endif
#include <setjmp.h>

//---------------------------------------------------------------------------

#include "ImageFormat_JPG.h"
using namespace LibSL::Image;

#include <LibSL/Errors/Errors.h>
using namespace LibSL::Errors;
#include <LibSL/Memory/Array.h>
using namespace LibSL::Memory::Array;
#include <LibSL/Memory/Pointer.h>
using namespace LibSL::Memory::Pointer;

//---------------------------------------------------------------------------

#define NAMESPACE LibSL::Image

//---------------------------------------------------------------------------

/// Declaring a global will automatically register the plugin
namespace {
  NAMESPACE::ImageFormat_JPG s_JPG;  /// FIXME: this mechanism does not work with VC++
}                                    ///        see also ImageFormatManager constructor

//---------------------------------------------------------------------------

int NAMESPACE::ImageFormat_JPG::s_JPEGQuality = 90;

//---------------------------------------------------------------------------

NAMESPACE::ImageFormat_JPG::ImageFormat_JPG()
{
  try {
    // register plugin
    IMAGE_FORMAT_MANAGER.registerPlugin(this);
  } catch (LibSL::Errors::Fatal& e) {
    std::cerr << e.message() << std::endl;
  }
}

//---------------------------------------------------------------------------

NAMESPACE::Image *NAMESPACE::ImageFormat_JPG::load(const char *name) const
{
  FILE *infile;
	fopen_s(&infile, name, "rb");
  if (infile == NULL)
    throw Fatal("ImageFormat_JPG::load - cannot open %s",name);

  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr         jerr;
  JSAMPARRAY     buffer;     // output row buffer
  int            row_stride; // physical row width in output buffer
  int            i,line = 0;
  unsigned int   width  = 0;
  unsigned int   height = 0;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, infile);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  width      = cinfo.output_width;
  height     = cinfo.output_height;
  row_stride = cinfo.output_width * cinfo.output_components;
  //row_stride;
  buffer = (*cinfo.mem->alloc_sarray)
    ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

  if (cinfo.output_components == 3) {
    ImageRGB *img = new ImageRGB(width,height);
    line=0;
    while (cinfo.output_scanline < cinfo.output_height)  {
      (void) jpeg_read_scanlines(&cinfo, buffer, 1);
      for (i=0;i<row_stride/3;i++) {
        img->pixel(i,line)[0] = buffer[0][i*3  ];
        img->pixel(i,line)[1] = buffer[0][i*3+1];
        img->pixel(i,line)[2] = buffer[0][i*3+2];
      }
      line++;
    }
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    return (img);
  } else if (cinfo.output_components == 1) {
    ImageL8 *img=new ImageL8(width,height);
    line=0;
    while (cinfo.output_scanline < cinfo.output_height)  {
      (void) jpeg_read_scanlines(&cinfo, buffer, 1);
      for (i=0;i<row_stride;i++) {
        img->pixel(i,line)[0] = buffer[0][i];
      }
      line++;
    }
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);    
    return (img);
  } else
    throw Fatal("ImageFormat_JPG::load - unsupported number of components (%d, %s)",cinfo.output_components,name);
}

//---------------------------------------------------------------------------

void NAMESPACE::ImageFormat_JPG::save(const char *fname,const NAMESPACE::Image *img) const
{
  // makes sure Tuple have the correct size (load scheme relies on pointers)
  sl_assert(sizeof(ImageRGB::t_Pixel) == sizeof(ImageRGB::t_Pixel::t_Element)*ImageRGB::t_Pixel::e_Size);
  
  const ImageRGBA *rgba = dynamic_cast<const ImageRGBA *>(img);
  const ImageRGB  *rgb  = dynamic_cast<const ImageRGB  *>(img);
  const ImageL8   *l8   = dynamic_cast<const ImageL8   *>(img);
  if (rgb == NULL && rgba == NULL && l8 == NULL) {
    throw Fatal("ImageFormat_JPG::save - PNG plugin only supports RGB, RGBA, L8 images (while saving '%s')",fname);
  }
  
  FILE *outfile;
	fopen_s(&outfile, fname, "wb");
  
  if (outfile == NULL) {
    throw Fatal("ImageFormat_JPG::save - cannot open %s",fname);
  }
  
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr       jerr;
  
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
  
  if (rgba) {
    cinfo.in_color_space   = JCS_RGB;
    cinfo.input_components = 4;
  } else if (rgb) {
    cinfo.in_color_space   = JCS_RGB;
    cinfo.input_components = 3;
  } else if (l8) {
    cinfo.in_color_space   = JCS_GRAYSCALE;
    cinfo.input_components = 1;
  } else {
    throw Fatal("ImageFormat_JPG::save - internal error");
  }
  cinfo.image_width  = img->w();
  cinfo.image_height = img->h();
  
  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo,s_JPEGQuality,TRUE);
  jpeg_default_colorspace(&cinfo);
	jpeg_stdio_dest(&cinfo, outfile);

  cinfo.dct_method         = JDCT_ISLOW;
  cinfo.optimize_coding    = TRUE;
  cinfo.write_Adobe_marker = FALSE;
  cinfo.write_JFIF_header  = FALSE;
  
	jpeg_start_compress(&cinfo, TRUE);
  
  JSAMPARRAY     buffer;
  buffer = (*cinfo.mem->alloc_sarray)((j_common_ptr) &cinfo, JPOOL_IMAGE, img->numComp()*img->w() , 1);
  while (cinfo.next_scanline < cinfo.image_height) {
    
    ForIndex(i,img->w()) {
      ForIndex(c,img->numComp()) {
        buffer[0][ i * img->numComp() + c ] = 
          img->raw()[ c + i * img->numComp() + cinfo.next_scanline * img->w() * img->numComp() ];
      }
    }
    jpeg_write_scanlines( &cinfo, buffer, 1 );
  }
  
  jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);
  
	fclose(outfile);
}

//---------------------------------------------------------------------------
