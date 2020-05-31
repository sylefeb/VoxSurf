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
// LibSL::Image::ImagePyramid
// ------------------------------------------------------
//
// ------------------------------------------------------
// Sylvain Lefebvre     - 2008-02-01
// ------------------------------------------------------

#pragma once

#include <LibSL/Image/Image.h>
#include <LibSL/System/Types.h>
#include <LibSL/Memory/Array.h>
#include <LibSL/Image/Filter.h>

namespace LibSL {
  namespace Image {

    //! ImagePyramid
    template <class T_Img>
    class ImagePyramid
    {
    protected:

      LibSL::Memory::Array::Array<typename T_Img::t_AutoPtr> m_Levels;

    public:

      typedef T_Img                                           t_Image;
      typedef typename T_Img::t_Pixel                         t_Pixel;
      typedef LibSL::Memory::Pointer::AutoPtr<ImagePyramid>   t_AutoPtr;

      ImagePyramid()                 
      {
        
      }

      ImagePyramid(const uint levels)
      {
        m_Levels.allocate(levels);
      }
      
      ImagePyramid(const uint w, const uint h)
      {
        uint levels = uint(LibSL::Math::log2((double)LibSL::Math::min<uint>(w, h))+1);
        m_Levels.allocate(levels);
      }

      ~ImagePyramid() {};

      uint numLevels() const {return m_Levels.size();}

      //! Accessors to single Level 
      const typename T_Img::t_AutoPtr& level(const uint l) const { return m_Levels[l];}
      typename T_Img::t_AutoPtr&       level(const uint l)       { return m_Levels[l];}
    };


    //! PyramidFactory
    template <class T_Img,template <class,uint,uint> class T_Filter,uint T_FilterSize=2>
    class ImagePyramidFactory
    {
    public:

      typedef LibSL::Filter::SeparableFilter2D< typename T_Img::t_PixelArray,T_Filter,T_FilterSize,2 > t_Filter;
    
    protected:

    public:

      ImagePyramidFactory() {};

      ImagePyramid<T_Img> *newPyramid(const typename T_Img::t_PixelArray& pixs)const
      {
        sl_assert((T_FilterSize % 2) == 0);

        t_Filter filter;

        ImagePyramid<T_Img> *pyr = new ImagePyramid<T_Img>(pixs.xsize(), pixs.ysize());
        pyr->level(0) = LibSL::Memory::Pointer::AutoPtr<T_Img>(new T_Img());
        pyr->level(0)->pixels() = pixs;
        ForRange(l, 1, pyr->numLevels() - 1) {
          LIBSL_BEGIN
            pyr->level(l) = LibSL::Memory::Pointer::AutoPtr<T_Img>(new T_Img());
          filter.filter(pyr->level(l - 1)->pixels(), pyr->level(l)->pixels());
          LIBSL_END
        }
        return (pyr);
      }
    };

  } //end namespace Image
} //end namespace LibSL
