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

#include "ImageFormat_TGA.h"
using namespace LibSL::Image;

#include <LibSL/Errors/Errors.h>
using namespace LibSL::Errors;
#include <LibSL/Memory/Array.h>
using namespace LibSL::Memory::Array;
#include <LibSL/Memory/Pointer.h>
using namespace LibSL::Memory::Pointer;
#include <LibSL/Math/Tuple.h>
using namespace LibSL::Math;

#include <iostream>
#include <string>
#include <cstring>
using namespace std;

#include "tga.h"

//---------------------------------------------------------------------------

#define NAMESPACE LibSL::Image

//---------------------------------------------------------------------------

/// Declaring a global will automatically register the plugin
namespace {
  NAMESPACE::ImageFormat_TGA s_TGA;  /// FIXME: this mechanism does not work with VC++
}                                    ///        see also ImageFormatManager constructor

//---------------------------------------------------------------------------

NAMESPACE::ImageFormat_TGA::ImageFormat_TGA()
{
  try {
    // register plugin
    IMAGE_FORMAT_MANAGER.registerPlugin(this);
  } catch (LibSL::Errors::Fatal& e) {
    std::cerr << e.message() << std::endl;
  }
}

//---------------------------------------------------------------------------

NAMESPACE::Image *NAMESPACE::ImageFormat_TGA::load(const char *name) const
{
  t_image_nfo *tga = ReadTGAFile(name);
  if (tga == NULL) {
    throw Fatal("ImageFormat_TGA::load - Sorry, cannot open file '%s'",name);
  }
  // read image data
  uint w = tga->width;
  uint h = tga->height;
  Image *img = NULL;
  if (tga->depth == 24) {
    // create image
    ImageRGB *rgb = new ImageRGB(w,h);
    img = rgb;
    // read pixels
    memcpy(rgb->pixels().raw(),tga->pixels,w*h*3);
  } else if (tga->depth == 32) {
    // create image
    ImageRGBA *rgba = new ImageRGBA(w,h);
    img = rgba;
    // read pixels
    memcpy(rgba->pixels().raw(),tga->pixels,w*h*4);
  } else {
    throw Fatal("ImageFormat_TGA::load - Sorry, unknown tga format (file '%s')",name);
  }
  // done
  delete[](tga->pixels);
  delete  (tga);
  return  (img);
}

//---------------------------------------------------------------------------

void NAMESPACE::ImageFormat_TGA::save(const char *name,const NAMESPACE::Image *img) const
{
#pragma pack(push, 1)
  /* TGA header */
  struct tga_header_t
  {
    uchar id_lenght;          /* size of image id */
    uchar colormap_type;      /* 1 if has a colormap */
    uchar image_type;         /* compression type */

    short	cm_first_entry;       /* colormap origin */
    short	cm_length;            /* colormap length */
    uchar cm_size;            /* colormap size */

    short	x_origin;             /* bottom left x coord origin */
    short	y_origin;             /* bottom left y coord origin */

    short	width;                /* picture width (in pixels) */
    short	height;               /* picture height (in pixels) */

    uchar pixel_depth;        /* bits per pixel: 8, 16, 24 or 32 */
    uchar image_descriptor;   /* 24 bits = 0x00; 32 bits = 0x80 */
  };
#pragma pack(pop)

  const ImageRGBA *rgba = dynamic_cast<const ImageRGBA *>(img);
  const ImageRGB  *rgb  = dynamic_cast<const ImageRGB  *>(img);
  const ImageL8   *lum  = dynamic_cast<const ImageL8   *>(img);
  if (rgba == NULL && rgb == NULL && lum == NULL) {
    throw Fatal("ImageFormat_TGA::save - Cannot save this format in file '%s'",name);
  }
  FILE *f = NULL;
	fopen_s(&f, name, "wb");
  if (f == NULL) {
    throw Fatal("ImageFormat_TGA::save - Sorry, cannot open file '%s'",name);
  }
  struct tga_header_t h;
  h.id_lenght        = 0;
  h.colormap_type    = 0;
  if (lum) {
    h.image_type     = 3;
  } else {
    h.image_type     = 2;
  }
  h.cm_first_entry   = 0;
  h.cm_length        = 0;
  h.cm_size          = 0;
  h.x_origin         = 0;
  h.y_origin         = 0;
  h.width            = img->w();
  h.height           = img->h();
  h.pixel_depth      = img->numComp() * 8;
  h.image_descriptor = (1<<5);
  fwrite(&h,sizeof(struct tga_header_t),1,f);
  const uchar *data = img->raw();
  if (lum) {
    fwrite(data,img->w()*img->h(),1,f);
  } else {
    ForImage(img,i,j) {
      ForIndex(c,3) { 
        fwrite(data + ((i + j* img->w()) * img->numComp() + (2-c)) ,1,1,f); // BGR ...
      }
      if (img->numComp() == 4) {
        fwrite(data + ((i + j* img->w()) * img->numComp() + 3) ,1,1,f); // alpha
      }
    }
  }
  fclose(f);

}

//---------------------------------------------------------------------------
