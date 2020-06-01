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
// LibSL::Memory::Array
// ------------------------------------------------------
//
// Array tools:
//  - save / load from file
// 
// ------------------------------------------------------
// Sylvain Lefebvre - 2009-10-05
// ------------------------------------------------------

#pragma once

// ------------------------------------------------------

#include <LibSL/Errors/Errors.h>
#include <LibSL/CppHelpers/CppHelpers.h>
#include <LibSL/System/Types.h>
#include <LibSL/Math/Tuple.h>
using namespace LibSL::System::Types;

#include <vector>

// ------------------------------------------------------

namespace LibSL  {
  namespace Memory {
    namespace Array {

      template <
        typename T_Type,
        template <typename> class P_Init,
        class P_Check
      >
      void saveArray(const Array<T_Type,P_Init,P_Check>& array,const char *fname)
      {
         FILE *f = NULL;
         fopen_s(&f,fname,"wb");
         if (f == NULL) {
           throw LibSL::Errors::Fatal("saveArray - Cannot open file '%s'",fname);
         }
         int    n  = array.size();
         size_t sz = fwrite(&n,sizeof(int),1,f);
         sl_assert( sz == 1 );
         sz        = fwrite(array.raw(),sizeof(T_Type),array.size(),f);
         sl_assert( sz == array.size() );
         fclose(f);
      }

      template <
        typename T_Type,
        template <typename> class P_Init,
        class P_Check
      >
      void loadArray(Array<T_Type,P_Init,P_Check>& _array,const char *fname)
      {
         FILE *f = NULL;
         fopen_s(&f,fname,"rb");
         if (f == NULL) {
           throw LibSL::Errors::Fatal("loadArray - Cannot open file '%s'",fname);
         }
         int    n  = 0;
         size_t sz = fread(&n,sizeof(int),1,f);
         sl_assert( sz == 1 );
         _array.allocate(n);
         sz        = fread(_array.raw(),sizeof(T_Type),_array.size(),f);
         sl_assert( sz == _array.size() );
         fclose(f);
      }

      template <
        typename T_Type,
        template <typename> class P_Init,
        class P_Check
      >
      void saveArray2D(const Array2D<T_Type,P_Init,P_Check>& array,const char *fname)
      {
         FILE *f = NULL;
         fopen_s(&f,fname,"wb");
         if (f == NULL) {
           throw LibSL::Errors::Fatal("saveArray - Cannot open file '%s'",fname);
         }
         int    x  = array.xsize();
         int    y  = array.ysize();
         size_t sx = fwrite(&x,sizeof(int),1,f);
         size_t sy = fwrite(&y,sizeof(int),1,f);
         sl_assert( sx == 1 );
         sl_assert( sy == 1 );
         size_t sz = fwrite(array.raw(),sizeof(T_Type),array.xsize()*array.ysize(),f);
         sl_assert( sz == array.xsize()*array.ysize() );
         fclose(f);
      }

      template <
        typename T_Type,
        template <typename> class P_Init,
        class P_Check
      >
      void loadArray2D(Array2D<T_Type,P_Init,P_Check>& _array,const char *fname)
      {
         FILE *f = NULL;
         fopen_s(&f,fname,"rb");
         if (f == NULL) {
           throw LibSL::Errors::Fatal("loadArray - Cannot open file '%s'",fname);
         }
         int    x  = 0;
         int    y  = 0;
         size_t sx = fread(&x,sizeof(int),1,f);
         size_t sy = fread(&y,sizeof(int),1,f);
         sl_assert( sx == 1 );
         sl_assert( sy == 1 );
         _array.allocate(x,y);
         size_t sz = fread(_array.raw(),sizeof(T_Type),_array.xsize()*_array.ysize(),f);
         sl_assert( sz == _array.xsize()*_array.ysize() );
         fclose(f);
      }


        template <
          typename T_Type,
          template <typename> class P_Init,
        class P_Check
        >
        void saveArray3D(const Array3D<T_Type, P_Init, P_Check>& array, const char *fname)
        {
          FILE *f = NULL;
          fopen_s(&f, fname, "wb");
          if (f == NULL) {
            throw LibSL::Errors::Fatal("saveArray - Cannot open file '%s'", fname);
          }
          int    x = array.xsize();
          int    y = array.ysize();
          int    z = array.zsize();
          size_t sx = fwrite(&x, sizeof(int), 1, f);
          size_t sy = fwrite(&y, sizeof(int), 1, f);
          size_t sz = fwrite(&z, sizeof(int), 1, f);
          sl_assert(sx == 1);
          sl_assert(sy == 1);
          sl_assert(sz == 1);
          size_t size = fwrite(array.raw(), sizeof(T_Type), array.xsize()*array.ysize()*array.zsize(), f);
          sl_assert(size == array.xsize()*array.ysize()*array.zsize());
          fclose(f);
        }

        template <
          typename T_Type,
          template <typename> class P_Init,
        class P_Check
        >
        void loadArray3D(Array3D<T_Type, P_Init, P_Check>& _array, const char *fname)
        {
          FILE *f = NULL;
          fopen_s(&f, fname, "rb");
          if (f == NULL) {
            throw LibSL::Errors::Fatal("loadArray - Cannot open file '%s'", fname);
          }
          int    x = 0;
          int    y = 0;
          int    z = 0;
          size_t sx = fread(&x, sizeof(int), 1, f);
          size_t sy = fread(&y, sizeof(int), 1, f);
          size_t sz = fread(&z, sizeof(int), 1, f);
          sl_assert(sx == 1);
          sl_assert(sy == 1);
          sl_assert(sz == 1);
          _array.allocate(x, y, z);
          size_t size = fread(_array.raw(), sizeof(T_Type), _array.xsize()*_array.ysize()*_array.zsize(), f);
          sl_assert(size == _array.xsize()*_array.ysize()*_array.ysize());
          fclose(f);
        }
        
        // TODO: ArrayND

    }
  }
}

// ----------------------------------------------------
