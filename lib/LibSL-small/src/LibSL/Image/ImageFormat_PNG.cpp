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
//---------------------------------------------------------------------------
#include "LibSL.precompiled.h"
//---------------------------------------------------------------------------

#include "ImageFormat_PNG.h"
using namespace LibSL::Image;

#include <LibSL/Errors/Errors.h>
using namespace LibSL::Errors;
#include <LibSL/Memory/Array.h>
using namespace LibSL::Memory::Array;
#include <LibSL/Memory/Pointer.h>
using namespace LibSL::Memory::Pointer;

//---------------------------------------------------------------------------

extern "C" {
#include "png.h"
}

//---------------------------------------------------------------------------

#define NAMESPACE LibSL::Image

//---------------------------------------------------------------------------

/// Declaring a global will automatically register the plugin
namespace {
  NAMESPACE::ImageFormat_PNG s_PNG;  /// FIXME: this mechanism does not work with VC++
}                                    ///        see also ImageFormatManager constructor

//---------------------------------------------------------------------------

NAMESPACE::ImageFormat_PNG::ImageFormat_PNG()
{
  try {
    // register plugin
    IMAGE_FORMAT_MANAGER.registerPlugin(this);
  } catch (LibSL::Errors::Fatal& e) {
    std::cerr << e.message() << std::endl;
  }
}

//---------------------------------------------------------------------------

NAMESPACE::Image *NAMESPACE::ImageFormat_PNG::load(const char *name) const
{
  FILE *file;
	fopen_s(&file, name, "rb");

  if (file == NULL) {
    throw Fatal("ImageFormat_PNG::load - cannot open %s",name);
  }

  png_structp png_ptr =
    png_create_read_struct(PNG_LIBPNG_VER_STRING,
    0, // (png_voidp)user_error_ptr
    0, // user_error_fn
    0  // user_warning_fn
    );

  png_infop info_ptr = png_create_info_struct(png_ptr);
  png_infop end_info = png_create_info_struct(png_ptr);
  png_init_io(png_ptr, file);

  png_read_info(png_ptr, info_ptr);
  uint w         = png_get_image_width(png_ptr,info_ptr);
  uint h         = png_get_image_height(png_ptr,info_ptr);
  int ncomp      = png_get_channels(png_ptr,info_ptr);
  int bit_depth  = png_get_bit_depth(png_ptr,info_ptr);
  int color_type = png_get_color_type(png_ptr,info_ptr);
  if (w <= 0 || h <=0)
    throw Fatal("ImageFormat_PNG::load - corrupted PNG file (%s)",name);
  if (bit_depth !=8 )
    throw Fatal("ImageFormat_PNG::load - unsupported bit depth (%d bit in %s)",ncomp,name);
  if (color_type == PNG_COLOR_TYPE_PALETTE)
    throw Fatal("ImageFormat_PNG::load - palette not supported (%s)",name);
  //  if (color_type == PNG_COLOR_TYPE_GRAY)
  //    png_set_gray_to_rgb(png_ptr);

  Array<unsigned char> tmp;
  tmp.allocate(h*w*ncomp);
  Array<png_bytep> rowptrs;
  rowptrs.allocate(h);
  for (uint j=0;j<h;j++) {
    rowptrs[j]=(png_bytep)&(tmp[(j*w)*ncomp]);
  }
  png_read_image(png_ptr, rowptrs.raw());
  png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
  fclose(file);

  if (ncomp == 3) {
    ImageRGB *rgb=new ImageRGB(w,h);
    ForPixels(rgb,i,j) {
      for (uint c=0;c<3;c++)
        rgb->pixel(i,j)[c]=tmp[(i+j*w)*ncomp+c];
    }
    return (rgb);
  } else if (ncomp == 4) {
    ImageRGBA *rgba=new ImageRGBA(w,h);
    ForPixels(rgba,i,j) {
      for (uint c=0;c<4;c++)
        rgba->pixel(i,j)[c]=tmp[(i+j*w)*ncomp+c];
    }
    return (rgba);
  } else if (ncomp == 1) {
    ImageL8 *l8=new ImageL8(w,h);
    ForPixels(l8,i,j) {
      l8->pixel(i,j)[0]=tmp[(i+j*w)*ncomp];
    }
    return (l8);
  } else
    throw Fatal("ImageFormat_PNG::load - unsupported number of components (%d components in %s)",ncomp,name);
}

//---------------------------------------------------------------------------


void NAMESPACE::ImageFormat_PNG::save(const char *fname,const NAMESPACE::Image *img) const
{
  // makes sure Tuple have the correct size (load scheme relies on pointers)
  sl_assert(sizeof(ImageRGB::t_Pixel) == sizeof(ImageRGB::t_Pixel::t_Element)*ImageRGB::t_Pixel::e_Size);

  const ImageRGBA *rgba = dynamic_cast<const ImageRGBA *>(img);
  const ImageRGB  *rgb  = dynamic_cast<const ImageRGB  *>(img);
  const ImageL8   *l8   = dynamic_cast<const ImageL8   *>(img);
  if (rgb == NULL && rgba == NULL && l8 == NULL) {
    throw Fatal("ImageFormat_PNG::save - PNG plugin only supports RGB, RGBA, L8 images (while saving '%s')",fname);
  }

  FILE *file;
	fopen_s(&file, fname, "wb");

  if (file == NULL)
    throw Fatal("Unable to save %s",fname);

  png_structp png_ptr  = png_create_write_struct(PNG_LIBPNG_VER_STRING,0,0,0);
  png_infop   info_ptr = png_create_info_struct(png_ptr);
  png_init_io(png_ptr, file);
  png_set_IHDR(png_ptr, info_ptr, img->w(), img->h(),
	       8,
	       (rgb  != NULL)  ? PNG_COLOR_TYPE_RGB
	       :(rgba != NULL) ? PNG_COLOR_TYPE_RGBA
	       :                 PNG_COLOR_TYPE_GRAY,
	       PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_BASE,
	       PNG_FILTER_TYPE_DEFAULT);
  // set max compression
  png_set_compression_level(png_ptr, 9);

  png_write_info(png_ptr, info_ptr);
  for (uint j=0;j<img->h();j++)  {
    png_bytep rowptr;
    if (rgb != NULL) {
      rowptr=(png_bytep)&(rgb->pixels().get(0,j)[0]); // NOTE: this could be dangerous (Tuple might not be byte aligned)
    } else if (rgba != NULL) {                        //       (see sl_assert above)
      rowptr=(png_bytep)&(rgba->pixels().get(0,j)[0]);// NOTE: same as above
    } else {
      rowptr=(png_bytep)&(l8->pixels().get(0,j)[0]);  // NOTE: same as above
    }
    png_write_row(png_ptr, rowptr);
  }
  png_write_end(png_ptr, NULL);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  fclose(file);
}

//---------------------------------------------------------------------------
