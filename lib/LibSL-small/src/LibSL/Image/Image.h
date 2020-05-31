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
// ------------------------------------------------------
// LibSL::Image
// ------------------------------------------------------
//
// Image class
//
// Coordinate system origin is at top left corner
// pixel(x,y) retrieves pixel at row y, column x
// ( or Array2D::get(x,y) )
//
// ------------------------------------------------------
// Sylvain Lefebvre - 2006-03-07
// ------------------------------------------------------

#pragma once

//#define FORCE_LINK_THIS(x) int force_link_##x = 0;

//#define FORCE_LINK_THAT(x) { extern int force_link_##x; force_link_##x = 1; }

//FORCE_LINK_THAT(s_PNG)

#include <LibSL/LibSL.common.h>

#include <LibSL/Errors/Errors.h>
#include <LibSL/Memory/Array2D.h>
#include <LibSL/Memory/Pointer.h>
#include <LibSL/Math/Tuple.h>
#include <LibSL/Math/Math.h>
#include <LibSL/Math/Vertex.h>
#include <LibSL/Memory/Pointer.h>
#include <LibSL/System/Types.h>

#include <map>

#define ForPixels(IMG,I,J)   for (uint J=0;J<IMG->h();J++) for (uint I=0;I<IMG->w();I++)
#define ForImage(IMG,I,J)    ForPixels(IMG,I,J)

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

namespace LibSL {
  namespace Image {

    /// Base Image class
    class Image
      //  : public LibSL::Memory::TraceLeaks::LeakProbe<Image>
    {
    public:

      virtual ~Image() {}

      virtual uint   w()           const =0;
      virtual uint   h()           const =0;
      virtual uint   numComp()     const =0;
      virtual uint   sizeOfComp()  const =0;
      virtual uchar*       raw()         =0;
      virtual const uchar* raw()   const =0;
    };

    typedef LibSL::Memory::Pointer::AutoPtr<Image> Image_Ptr;

    /// Image format plugin (abstract)
    class LIBSL_DLL ImageFormat_plugin
      //  : public LibSL::Memory::TraceLeaks::LeakProbe<ImageFormat_plugin>
    {
    public:
      virtual void        save(const char *,const Image *) const =0;
      virtual Image      *load(const char *)               const =0;
      virtual const char *signature()                      const =0;
      virtual ~ImageFormat_plugin() {}
    };

    /// Image format manager (singleton)
    class LIBSL_DLL ImageFormatManager
    {
    private:

      std::map<std::string,ImageFormat_plugin*> m_Plugins;

      static ImageFormatManager *s_Manager;
      ImageFormatManager();

    public:

      ~ImageFormatManager();

      void registerPlugin(ImageFormat_plugin *plugin);
      ImageFormat_plugin *getPlugin(const char *signature) const;

      static ImageFormatManager *getUniqueInstance();

    };

    /// Generic Image class
    template <typename T_Type,unsigned int T_NumComp>
    class Image_generic : public Image
    {
    public:

      typedef T_Type                                         t_Component;
      typedef LibSL::Math::Tuple<T_Type,T_NumComp>           t_Pixel;
      typedef LibSL::Memory::Array::Array2D<t_Pixel>         t_PixelArray;
      typedef LibSL::Memory::Pointer::AutoPtr<Image_generic> t_AutoPtr;

      enum {e_NumComp = T_NumComp};

    protected:

      t_PixelArray m_PixelArray;

    public:

      Image_generic()
      {
      }

      Image_generic(uint w, uint h, const T_Type& init = 0)
      {
        m_PixelArray.allocate(w,h);
        t_Pixel z; z = init;
        m_PixelArray.fill(z);
      }

      Image_generic(uint w, uint h, const t_Pixel& init)
      {
        m_PixelArray.allocate(w,h);
        m_PixelArray.fill(init);
      }

      Image_generic(const Image_generic& img)
      {
        m_PixelArray=img.pixels();
      }

      ~Image_generic()
      {

      }

      /// Accessors to single pixel
      const t_Pixel& pixel(uint i, uint j) const { return m_PixelArray.at(i,j); }
      t_Pixel&       pixel(uint i, uint j)       { return m_PixelArray.at(i,j); }


      /// Accessors to single pixel with access Policy
      template <class T_AccessPolicy>
      const t_Pixel& pixel(uint i, uint j) const {
        return m_PixelArray. template at<T_AccessPolicy>(i,j);
      }
      template <class T_AccessPolicy>
      t_Pixel& pixel(uint i, uint j)       {
        return m_PixelArray. template at<T_AccessPolicy>(i,j);
      }

      /// Accessors to pixel array
      const t_PixelArray&  pixels() const {return (m_PixelArray);}
      t_PixelArray&        pixels()       {return (m_PixelArray);}

      /// Image sizes
      uint w() const            {return (m_PixelArray.xsize());}
      uint h() const            {return (m_PixelArray.ysize());}
      uint numComp()    const   {return (T_NumComp);           }
      uint sizeOfComp() const   {return (sizeof(t_Component)); }

      /// low level memory access
      uchar*        raw()       {return reinterpret_cast<uchar*      >(pixels().raw());}
      const uchar*  raw() const {return reinterpret_cast<const uchar*>(pixels().raw());}

      /// Copy an image at the given location
      void copy(int x, int y, const Image_generic* img)
      {
        ForImage(img, i, j) {
          int xd = i + x;
          int yd = j + y;
          if ((xd >= 0 && xd < int(w())) && (yd >= 0 && yd < int(h()))) {
            pixel(xd, yd) = img->pixel(i, j);
          }
        }
      }

      /// Copy an image at the given location
      void copy(int x,int y,const t_AutoPtr& img)
      {
        copy(x, y, img.raw());
      }

      /// Extract a sub-image
      template <class T_ImageAccessPolicy>
      Image_generic *extract(int xe,int ye,uint we,uint he) const
      {
        sl_assert(we > 0 && he > 0);
        Image_generic *sub = new Image_generic(we,he);
        ForImage(sub,i,j) {
          int xs = xe+i;
          int ys = ye+j;
          sub->pixel(i,j) = pixel<T_ImageAccessPolicy>(xs,ys);
        }
        return (sub);
      }

      Image_generic *extract(int xe,int ye,uint we,uint he) const
      {
        return (extract<Memory::Array::Wrap>(xe,ye,we,he));
      }


      /// Flip horizontaly
      void flipH()
      {
        t_Pixel tmp;
        ForIndex(j, h() / 2) {
          ForIndex(i, w()) {
            tmp = pixels().at(i, j);
            pixels().at(i, j) = pixels().at(i, h() - 1 - j);
            pixels().at(i, h() - 1 - j) = tmp;
          }
        }
      }

      /// Flip vertically
      void flipV()
      {
        t_Pixel tmp;
        ForIndex(j, h()) {
          ForIndex(i, w() / 2) {
            tmp = pixels().at(i, j);
            pixels().at(i, j) = pixels().at(w() - 1 - i, j);
            pixels().at(w() - 1 - i, j) = tmp;
          }
        }
      }

      /// Clamp pixels
      void clamp(const T_Type& minClamp,const T_Type& maxClamp)
      {
        ForImage(this,i,j) {
          pixel(i,j) = LibSL::Math::clamp(pixel(i,j),Math::Tuple<T_Type,e_NumComp>(minClamp),Math::Tuple<T_Type,e_NumComp>(maxClamp));
        }
      }

      /// Compute min/max over image
      void findMinMax(
        t_Pixel&       minImage,
        t_Pixel&       maxImage)
      {
        minImage = t_Component( 1e30f);
        maxImage = t_Component(-1e30f);
        ForImage(this,x,y) {
          minImage = LibSL::Math::tupleMin(pixel(x,y),minImage);
          maxImage = LibSL::Math::tupleMax(pixel(x,y),maxImage);
        }
      }

      /// Remap pixel values between minVal and maxVal
      void remap(
        const t_Pixel& minVal,
        const t_Pixel& maxVal)
      {
        t_Pixel minImage;
        t_Pixel maxImage;
        findMinMax(minImage,maxImage);
		sl_assert( sqLength( maxImage - minImage ) > 0 );
        ForImage(this,x,y) {
          pixel(x,y) = minVal + ((maxVal - minVal) * (pixel(x,y) - minImage))/(maxImage - minImage);
        }
      }

      // FIXME: vulnerability to 0 divide
      void remapFromTo(const t_Pixel& minImage,const t_Pixel& maxImage,const t_Pixel& minVal,const t_Pixel& maxVal)
      {
        ForImage(this,x,y) {
          pixel(x,y) = minVal + (maxVal - minVal) * (pixel(x,y) - minImage)/(maxImage - minImage);
        }
      }

      /// Rotate 90 degrees clockwise
      void rotateCW()
      {
        t_PixelArray rotatedArray;
        rotatedArray.allocate(h(),w());
        int hm = h()-1;

        ForIndex(j,h()) {
          ForIndex(i,w()) {
            rotatedArray.at(hm-j,i) = pixels().at(i,j);
          }
        }
        m_PixelArray.erase();
        m_PixelArray = rotatedArray;
      }

      /// Rotate 90 degrees counter clockwise
      void rotateCCW()
      {
        t_PixelArray rotatedArray;
        rotatedArray.allocate(h(),w());
        int wm = w()-1;

        ForIndex(j,h()) {
          ForIndex(i,w()) {
            rotatedArray.at(j,wm-i) = pixels().at(i,j);
          }
        }
        m_PixelArray.erase();
        m_PixelArray = rotatedArray;
      }

      /// Bilinearly interpolate the image at a given location
      //  - coordinate range is [0..1] similarly to GPU textures
      //  - texels are at location (i+0.5)/w() and (j+0.5)/h() where i,j are integers
      //  - access policies are configurable with policies Clamp, Warp, ...
      //
      // Usage: img->bilinear<Clamp>(u,v);
      // Note : floor(-i), with i>0, is expected to be negative (as should be!)
      template <class T_ImageAccessPolicy>
      t_Pixel bilinear(const float u, const float v) const
      {
        // translate into "Array coordinate space"
        float i  = u * w() - 0.5f;
        float j  = v * h() - 0.5f;

        // get corners of interpolation cell
        int i0   = int(floor(i));
        int i1   = i0 + 1;
        int j0   = int(floor(j));
        int j1   = j0 + 1;

        // apply access policy
        i0 = T_ImageAccessPolicy::access(i0,w());
        i1 = T_ImageAccessPolicy::access(i1,w());
        j0 = T_ImageAccessPolicy::access(j0,h());
        j1 = T_ImageAccessPolicy::access(j1,h());

        // compute interpolation fractions
        float fi = (i - int(floor(i)));
        float fj = (j - int(floor(j)));

        // interpolate each component and return
        t_Pixel pix;
        ForIndex(c, t_Pixel::e_Size){
          pix[c] = t_Component(
            (1.0f-fi)*((1.0f-fj)*pixel(i0,j0)[c]
         +  (     fj)*           pixel(i0,j1)[c])
         +  (     fi)*((1.0f-fj)*pixel(i1,j0)[c]
         +  (     fj)*           pixel(i1,j1)[c]));
        }
        return pix;
      }

      t_Pixel bilinear(const float u, const float v) const
      {
        return (bilinear<Memory::Array::Clamp>(u,v));
      }



      /// Nearest Neighbor interpolation of the image at a given location
      //  - coordinate range is [0..1] similarly to GPU textures
      //  - texels are at location (i+0.5)/w() and (j+0.5)/h() where i,j are integers
      //  - access policies are configurable with policies Clamp, Warp, ...
      // Usage: img->nearest<Clamp>(u,v);
      // Note:  floor(-i), with i>0, is expected to be negative (as should be!)
      template <class T_ImageAccessPolicy>
      t_Pixel nearest(const float u, const float v) const
      {
        // translate into "Array coordinate space"
        int i  = int(u * w());// equivalent to round(u * w() - 0.5f)
        int j  = int(v * h());

        // apply access policy
        i = T_ImageAccessPolicy::access(i,w());
        j = T_ImageAccessPolicy::access(j,h());

        return pixel(i,j);
      }

      t_Pixel nearest(const float u, const float v) const
      {
        return (nearest<Memory::Array::Clamp>(u,v));
      }



      /// Bicubic interpolate the image at a given location
      //  - coordinate range is [0..1] similarly to GPU textures
      //  - texels are at location (i+0.5)/w() and (j+0.5)/h() where i,j are integers
      //  - access policies are configurable with policies Clamp, Warp, ...
      // Usage: img->bicubic<Clamp>(u,v);
      // Note:  floor(-i), with i>0, is expected to be negative (as should be!)
      template <class T_ImageAccessPolicy>
      t_Pixel bicubic(const float u, const float v) const
      {
        // translate into "Array coordinate space"
        float i  = u * w() - 0.5f;
        float j  = v * h() - 0.5f;

        // get corners of interpolation cell
        int i0   = int(floor(i))-1;
        int j0   = int(floor(j))-1;
        int i1   = i0+1;
        int j1   = j0+1;
        int i2   = i0+2;
        int j2   = j0+2;
        int i3   = i0+3;
        int j3   = j0+3;

        // apply access policy
        i0 = T_ImageAccessPolicy::access(i0,w());
        j0 = T_ImageAccessPolicy::access(j0,h());
        i1 = T_ImageAccessPolicy::access(i1,w());
        j1 = T_ImageAccessPolicy::access(j1,h());
        i2 = T_ImageAccessPolicy::access(i2,w());
        j2 = T_ImageAccessPolicy::access(j2,h());
        i3 = T_ImageAccessPolicy::access(i3,w());
        j3 = T_ImageAccessPolicy::access(j3,h());

        // compute interpolation fractions
        float fi = (i - int(floor(i)));
        float fj = (j - int(floor(j)));

        float s0 = ((2-fi)*fi-1)*fi;    // -1
        float s1 = (3*fi-5)*fi*fi+2;    //  0
        float s2 = ((4-3*fi)*fi+1)*fi;  // +1
        float s3 = (fi-1)*fi*fi;        // +2

        float t0 = ((2-fj)*fj-1)*fj;
        float t1 = (3*fj-5)*fj*fj+2;
        float t2 = ((4-3*fj)*fj+1)*fj;
        float t3 = (fj-1)*fj*fj;

        LibSL::Math::Tuple<float,e_NumComp> pi0 =
          s0*pixel(i0,j0) + s1*pixel(i1,j0) + s2*pixel(i2,j0) + s3*pixel(i3,j0);
        LibSL::Math::Tuple<float,e_NumComp> pi1 =
          s0*pixel(i0,j1) + s1*pixel(i1,j1) + s2*pixel(i2,j1) + s3*pixel(i3,j1);
        LibSL::Math::Tuple<float,e_NumComp> pi2 =
          s0*pixel(i0,j2) + s1*pixel(i1,j2) + s2*pixel(i2,j2) + s3*pixel(i3,j2);
        LibSL::Math::Tuple<float,e_NumComp> pi3 =
          s0*pixel(i0,j3) + s1*pixel(i1,j3) + s2*pixel(i2,j3) + s3*pixel(i3,j3);

        return t_Pixel(0.25f * ( pi0*t0 + pi1*t1 + pi2*t2 + pi3*t3 ));
      }

      t_Pixel bicubic(const float u, const float v) const
      {
        return (bicubic<Memory::Array::Clamp>(u,v));
      }


      /// Cast into another image type
      template <class T_Image> T_Image *cast() const
      {
        T_Image *nimg = new T_Image(w(),h());
        ForImage(nimg,i,j) {
          ForIndex(c,Math::min(int(e_NumComp),int(T_Image::e_NumComp))) {
            nimg->pixel(i,j)[c]=(typename T_Image::t_Component)pixel(i,j)[c];
          }
        }
        return (nimg);
      }

      /// Clone image
      Image_generic *clone() const
      {
        Image_generic *nimg = new Image_generic(w(),h());
        ForImage(nimg,i,j) {
          nimg->pixel(i,j) = pixel(i,j);
        }
        return (nimg);
      }

    };

    /// Load and save global methods

    LIBSL_DLL Image    *loadImage(const char *);
    LIBSL_DLL void      saveImage(const char *,const Image_Ptr&);
    LIBSL_DLL void      saveImage(const char *,const Image *);

    template<class T_Image>
    LIBSL_DLL void  saveImage(const char *name, const LibSL::Memory::Pointer::AutoPtr<T_Image>& ptr)
    {
      saveImage(name, ptr.raw());
    }

    /// Load an image and cast it

    template<class T_Image> T_Image*
      loadImage(const char *fname)
    {
      Image   *img = loadImage(fname);
      T_Image *ptr = dynamic_cast<T_Image *>(img);
      if (ptr == NULL) {
        delete (img); // delete image!
        throw LibSL::Errors::Fatal("loadImage<> - Image format mismatch! (%s)",fname);
        return NULL;
      }
      return (ptr);
    }

    //! random colors
    LIBSL_DLL LibSL::Math::v3f randomColorFromIndex(uint idx);

    /// Standard image types

    typedef Image_generic<unsigned char,3>            ImageRGB;
    typedef ImageRGB::t_AutoPtr                       ImageRGB_Ptr;

    typedef Image_generic<unsigned char,4>            ImageRGBA;
    typedef ImageRGBA::t_AutoPtr                      ImageRGBA_Ptr;

    typedef Image_generic<unsigned char,1>            ImageL8;
    typedef ImageL8::t_AutoPtr                        ImageL8_Ptr;

    typedef Image_generic<unsigned char,2>            ImageUV8;
    typedef ImageUV8::t_AutoPtr                       ImageUV8_Ptr;

    typedef Image_generic<float,3>                    ImageRGB32F;
    typedef ImageRGB32F::t_AutoPtr                    ImageRGB32F_Ptr;
    typedef Image_generic<float,3>                    ImageFloat3;
    typedef ImageFloat3::t_AutoPtr                    ImageFloat3_Ptr;

    typedef Image_generic<float,4>                    ImageRGBA32F;
    typedef ImageRGBA32F::t_AutoPtr                   ImageRGBA32F_Ptr;
    typedef Image_generic<float,4>                    ImageFloat4;
    typedef ImageFloat4::t_AutoPtr                    ImageFloat4_Ptr;

    typedef Image_generic<float,1>                    ImageL32F;
    typedef ImageL32F::t_AutoPtr                      ImageL32F_Ptr;
    typedef Image_generic<float,1>                    ImageFloat1;
    typedef ImageFloat1::t_AutoPtr                    ImageFloat1_Ptr;

    typedef Image_generic<float,2>                    ImageFloat2;
    typedef ImageFloat2::t_AutoPtr                    ImageFloat2_Ptr;

    typedef Image_generic<half,3>                     ImageRGB16F;
    typedef ImageRGB16F::t_AutoPtr                    ImageRGB16F_Ptr;

    typedef Image_generic<half,4>                     ImageRGBA16F;
    typedef ImageRGBA16F::t_AutoPtr                   ImageRGBA16F_Ptr;

    typedef Image_generic<half,1>                     ImageL16F;
    typedef ImageL16F::t_AutoPtr                      ImageL16F_Ptr;

    typedef Image_generic<bool,1>                     ImageBool1;
    typedef ImageBool1::t_AutoPtr                     ImageBool1_Ptr;

    typedef Image_generic<double,1>                    ImageDouble1;
    typedef ImageDouble1::t_AutoPtr                    ImageDouble1_Ptr;

    typedef Image_generic<double,2>                    ImageDouble2;
    typedef ImageDouble2::t_AutoPtr                    ImageDouble2_Ptr;

    typedef Image_generic<double,3>                    ImageDouble3;
    typedef ImageDouble1::t_AutoPtr                    ImageDouble3_Ptr;

    typedef Image_generic<double,4>                    ImageDouble4;
    typedef ImageDouble2::t_AutoPtr                    ImageDouble4_Ptr;


  } //namespace LibSL::Image
} //namespace LibSL

// ------------------------------------------------------

#define IMAGE_FORMAT_MANAGER (*LibSL::Image::ImageFormatManager::getUniqueInstance())

// ------------------------------------------------------
