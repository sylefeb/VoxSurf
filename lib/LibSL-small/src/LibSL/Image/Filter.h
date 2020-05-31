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
// LibSL::Filter
// ------------------------------------------------------ 
// 
// ------------------------------------------------------
// Sylvain Lefebvre     - 2008-02-01
// ------------------------------------------------------

#pragma once
#include <LibSL/Image/Image.h>
#include <LibSL/System/Types.h>
#include <LibSL/Memory/Array.h>

namespace LibSL {
  namespace Filter {

    // Filters along 1 dimension
    template <class T_Array,uint T_DownSampleFactor=1>
    class Filter1D
    {
    public:
      typedef T_Array t_Array;
      enum {e_DownSampleFactor = T_DownSampleFactor};
    public:
      virtual void filter(const T_Array& source,T_Array& _dest) const = 0;
    };

    //! Weighted 1D filter (dual filter)
    template <class T_Array,uint T_DownSampleFactor=1>
    class WeightedFilter1D
    {
    protected:
      LibSL::Memory::Array::Array<typename T_Array::t_Element> m_Weights; // filled by specialized sub-classes
    public:
      WeightedFilter1D()
      { }
      void filter(const T_Array& source,T_Array& _dest) const
      {
        LIBSL_BEGIN
          _dest.allocate(source.size() / T_DownSampleFactor);
        ForIndex(x,_dest.size()) {
          typename T_Array::t_Element sum       = 0;
          typename T_Array::t_Element sumWeight = 0;
          int  s = int(m_Weights.size()) / 2;
          ForIndex(n,m_Weights.size()) {
            typename T_Array::t_Element weight = typename T_Array::t_Element(m_Weights[n]);
            int i                              = int(x * T_DownSampleFactor + n) - s + int((m_Weights.size()&1)?0:1);
            if (i >= 0 && i < (int)source.size()) {
              sum                    += source.at(uint(i)) * weight;
              sumWeight              += weight;
            }
          }
          _dest.at(x) = sum / sumWeight;
        }
        LIBSL_END
      }
    };

    //! Filters along 2 dimensions
    template <class T_Array2D,uint T_DownSampleFactor=1>
    class Filter2D
    {
    public:
      typedef T_Array2D t_Array;
      enum {e_DownSampleFactor = T_DownSampleFactor};
    public:
      virtual void filter(const T_Array2D& source,T_Array2D& _dest) const = 0;
    };

    //! Select a row of a 2D array as 1D array
    template <class T_Array2D>
    class RowSelect
    {
    private:
      T_Array2D& m_Array2D;
      uint       m_Row;
    public:
      typedef typename T_Array2D::t_Element t_Element;
    public:
      RowSelect(T_Array2D& array2d,uint row) : m_Array2D (array2d) 
      {
        m_Row = row;
      }
      void                           allocate(uint sz)   { }
      uint                           size() const        { return m_Array2D.xsize();     }
      typename T_Array2D::t_Element& at(uint i)          { return m_Array2D.at(i,m_Row); }
      typename T_Array2D::t_Element& operator[] (uint i) { return m_Array2D.at(i,m_Row); }
      const typename T_Array2D::t_Element& at(uint i)         const { return m_Array2D.at(i,m_Row); }
      const typename T_Array2D::t_Element& operator[](uint i) const { return m_Array2D.at(i,m_Row); }
    };

    //! Select a column of a 2D array as 1D array
    template <class T_Array2D>
    class ColumnSelect
    {
    private:
      T_Array2D& m_Array2D;
      uint       m_Col;
    public:
      typedef typename T_Array2D::t_Element t_Element;
    public:
      ColumnSelect(T_Array2D& array2d,uint col) : m_Array2D (array2d) 
      {
        m_Col = col;
      }
      void                           allocate(uint sz)   { }
      uint                           size() const        { return m_Array2D.ysize();     }
      typename T_Array2D::t_Element& at(uint i)          { return m_Array2D.at(m_Col,i); }
      typename T_Array2D::t_Element& operator[](uint i)  { return m_Array2D.at(m_Col,i); }
      const typename T_Array2D::t_Element& at(uint i)         const { return m_Array2D.at(m_Col,i); }
      const typename T_Array2D::t_Element& operator[](uint i) const { return m_Array2D.at(m_Col,i); }
    };

    //! SeparableFilter2D
    //! Uses a 1D filter on rows and then columns
    template < class T_Array2D, template <class,uint,uint> class T_Filter, uint T_Size, uint T_DownSampleFactor=1 >
    class SeparableFilter2D : public Filter2D<T_Array2D,T_DownSampleFactor>
    {
    public:
      typedef T_Array2D                                                  t_Array;
      typedef T_Filter<RowSelect<t_Array>   ,T_Size,T_DownSampleFactor>  t_RowFilter;
      typedef T_Filter<ColumnSelect<t_Array>,T_Size,T_DownSampleFactor>  t_ColFilter;
    public:
      void filter(const t_Array& source,t_Array& _dest) const
      {
        t_RowFilter rfilter;
        t_ColFilter cfilter;
        t_Array     tmp;
        tmp.allocate(source.xsize()/T_DownSampleFactor,source.ysize());
        // filter rows
        ForIndex(n,source.ysize()) {
          RowSelect<t_Array>  src_row(const_cast<t_Array&>(source),n);
          RowSelect<t_Array>  dst_row(tmp   ,n);
          rfilter.filter(src_row, dst_row);
        }
        // filter columns
        _dest.allocate(source.xsize()/T_DownSampleFactor,source.ysize()/T_DownSampleFactor);
        ForIndex(n,tmp.xsize()) {
          ColumnSelect<t_Array>  src_row(tmp  ,n);
          ColumnSelect<t_Array>  dst_row(_dest,n);
          cfilter.filter(src_row, dst_row);
        }
      }
    };

    // --------------------------------------------------
    // Filter implementations
    // --------------------------------------------------

    //! 1D BoxFilter
    template <class T_Array,uint T_Size,uint T_DownSampleFactor=1>
    class BoxFilter : public WeightedFilter1D<T_Array,T_DownSampleFactor>
    {
    protected:
      void fillInWeights()
      {
        this->m_Weights.allocate(T_Size);
        ForIndex(n,this->m_Weights.size()) {
          this->m_Weights[n] = typename T_Array::t_Element(1);
        }
      }
    public:
      typedef T_Array t_Array;
      enum {e_DownSampleFactor = T_DownSampleFactor};
    public:
      BoxFilter() : WeightedFilter1D<T_Array,T_DownSampleFactor>()
      {  fillInWeights();  }
    };

    //! 1D Gaussian filter
    template <class T_Array,uint T_Size,uint T_DownSampleFactor=1>
    class GaussianFilter : public WeightedFilter1D<T_Array,T_DownSampleFactor>
    {
    protected:

      double GaussFW1(double x) const
      {
        // Gaussian with standard deviation choosen to accomodate "full width at half height" = 1
        // see: http://mathworld.wolfram.com/FullWidthatHalfMaximum.html
        double w = 1.0; // full width at half height should be 1.
        double s = w / (2.0*sqrt(2.0f*log(2.0))); 
        double fract = 1.0 / (s * sqrt(2.0 * M_PI));
        return fract * exp( - x*x / (2.0*s*s) );
      }

      void fillInWeights()
      {
        int   sz;
        if ((T_Size&1) == 0) { sz = Math::max(4,T_Size); }
        else                 { sz = Math::max(3,T_Size); }
        this->m_Weights.allocate(sz);
        float freq = 2.0f / float(sz);
        ForIndex(n,this->m_Weights.size()) {
          this->m_Weights[n] = float(GaussFW1((n - float(sz)/2.0f + 0.5f) * freq / 1.1f)); // 10% safety distance
        }
      }

    public:
      typedef T_Array t_Array;
      enum {e_DownSampleFactor = T_DownSampleFactor};
    public:
      GaussianFilter() : WeightedFilter1D<T_Array,T_DownSampleFactor>()
      {  fillInWeights();  }
    };

    //! 2D Box filter
    template <class T_Array,uint T_Size,uint T_DownSampleFactor=1>
    class BoxFilter2D : public SeparableFilter2D<T_Array,BoxFilter,T_Size,T_DownSampleFactor>
    {
    public:
      BoxFilter2D() : SeparableFilter2D<T_Array,BoxFilter,T_Size,T_DownSampleFactor>()
      { }
    };

    //! 2D Gaussian filter
    template <class T_Array,uint T_Size,uint T_DownSampleFactor=1>
    class GaussianFilter2D : public SeparableFilter2D<T_Array,GaussianFilter,T_Size,T_DownSampleFactor>
    {
    public:
      GaussianFilter2D() : SeparableFilter2D<T_Array,GaussianFilter,T_Size,T_DownSampleFactor>()
      { }
    };

  } //end namespace Filter
} //end namespace LibSL
