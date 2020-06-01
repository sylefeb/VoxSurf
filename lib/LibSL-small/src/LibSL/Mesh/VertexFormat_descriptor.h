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

#pragma once

namespace LibSL {
  namespace Mesh {

    template <class VertexFormat> class FormatDescriptor
    {

    public:

      enum {idx_position=MVF_index_item<VertexFormat,MVF_BASE_POSITION>::value};
      enum {idx_normal=MVF_index_item<VertexFormat,MVF_BASE_NORMAL>::value};
      enum {idx_color0=MVF_index_item<VertexFormat,MVF_BASE_COLOR0>::value};
      enum {idx_color1=MVF_index_item<VertexFormat,MVF_BASE_COLOR1>::value};
      enum {idx_texcoord0=MVF_index_item<VertexFormat,MVF_BASE_TEXCOORD0>::value};
      enum {idx_texcoord1=MVF_index_item<VertexFormat,MVF_BASE_TEXCOORD1>::value};
      enum {idx_texcoord2=MVF_index_item<VertexFormat,MVF_BASE_TEXCOORD2>::value};
      enum {idx_texcoord3=MVF_index_item<VertexFormat,MVF_BASE_TEXCOORD3>::value};
      enum {idx_texcoord4=MVF_index_item<VertexFormat,MVF_BASE_TEXCOORD4>::value};
      enum {idx_texcoord5=MVF_index_item<VertexFormat,MVF_BASE_TEXCOORD5>::value};
      enum {idx_texcoord6=MVF_index_item<VertexFormat,MVF_BASE_TEXCOORD6>::value};
      enum {idx_texcoord7=MVF_index_item<VertexFormat,MVF_BASE_TEXCOORD7>::value};

      enum {has_position=(idx_position > -1)};
      enum {has_normal=(idx_normal > -1)};
      enum {has_color0=(idx_color0 > -1)};
      enum {has_color1=(idx_color1 > -1)};
      enum {has_texcoord0=(idx_texcoord0 > -1)};
      enum {has_texcoord1=(idx_texcoord1 > -1)};
      enum {has_texcoord2=(idx_texcoord2 > -1)};
      enum {has_texcoord3=(idx_texcoord3 > -1)};
      enum {has_texcoord4=(idx_texcoord4 > -1)};
      enum {has_texcoord5=(idx_texcoord5 > -1)};
      enum {has_texcoord6=(idx_texcoord6 > -1)};
      enum {has_texcoord7=(idx_texcoord7 > -1)};

      enum {size_of=MVF_sizeof<VertexFormat>::value};

      enum {offset_position =MVF_offset_item<VertexFormat,MVF_BASE_POSITION>::value};
      enum {offset_normal   =MVF_offset_item<VertexFormat,MVF_BASE_NORMAL>::value};
      enum {offset_color0   =MVF_offset_item<VertexFormat,MVF_BASE_COLOR0>::value};
      enum {offset_color1   =MVF_offset_item<VertexFormat,MVF_BASE_COLOR1>::value};
      enum {offset_texcoord0=MVF_offset_item<VertexFormat,MVF_BASE_TEXCOORD0>::value};
      enum {offset_texcoord1=MVF_offset_item<VertexFormat,MVF_BASE_TEXCOORD1>::value};
      enum {offset_texcoord2=MVF_offset_item<VertexFormat,MVF_BASE_TEXCOORD2>::value};
      enum {offset_texcoord3=MVF_offset_item<VertexFormat,MVF_BASE_TEXCOORD3>::value};
      enum {offset_texcoord4=MVF_offset_item<VertexFormat,MVF_BASE_TEXCOORD4>::value};
      enum {offset_texcoord5=MVF_offset_item<VertexFormat,MVF_BASE_TEXCOORD5>::value};
      enum {offset_texcoord6=MVF_offset_item<VertexFormat,MVF_BASE_TEXCOORD6>::value};
      enum {offset_texcoord7=MVF_offset_item<VertexFormat,MVF_BASE_TEXCOORD7>::value};

      // TODO FIXME dangerous use of -1 and unsigned int ; replace by int max

      typedef typename Loki::TL::TypeAtNonStrict<VertexFormat,static_cast<unsigned int>(idx_position),MVF_empty>::Result  type_position;
      typedef typename Loki::TL::TypeAtNonStrict<VertexFormat,static_cast<unsigned int>(idx_normal),MVF_empty>::Result    type_normal;
      typedef typename Loki::TL::TypeAtNonStrict<VertexFormat,static_cast<unsigned int>(idx_color0),MVF_empty>::Result    type_color0;
      typedef typename Loki::TL::TypeAtNonStrict<VertexFormat,static_cast<unsigned int>(idx_color1),MVF_empty>::Result    type_color1;
      typedef typename Loki::TL::TypeAtNonStrict<VertexFormat,static_cast<unsigned int>(idx_texcoord0),MVF_empty>::Result type_texcoord0;
      typedef typename Loki::TL::TypeAtNonStrict<VertexFormat,static_cast<unsigned int>(idx_texcoord1),MVF_empty>::Result type_texcoord1;
      typedef typename Loki::TL::TypeAtNonStrict<VertexFormat,static_cast<unsigned int>(idx_texcoord2),MVF_empty>::Result type_texcoord2;
      typedef typename Loki::TL::TypeAtNonStrict<VertexFormat,static_cast<unsigned int>(idx_texcoord3),MVF_empty>::Result type_texcoord3;
      typedef typename Loki::TL::TypeAtNonStrict<VertexFormat,static_cast<unsigned int>(idx_texcoord4),MVF_empty>::Result type_texcoord4;
      typedef typename Loki::TL::TypeAtNonStrict<VertexFormat,static_cast<unsigned int>(idx_texcoord5),MVF_empty>::Result type_texcoord5;
      typedef typename Loki::TL::TypeAtNonStrict<VertexFormat,static_cast<unsigned int>(idx_texcoord6),MVF_empty>::Result type_texcoord6;
      typedef typename Loki::TL::TypeAtNonStrict<VertexFormat,static_cast<unsigned int>(idx_texcoord7),MVF_empty>::Result type_texcoord7;

      //   xx   xxxx
      //  x x      x
      // xxxxxx  xx
      //    x   xxxxx

      static int size(int i)
      {
        static const int tbl[]={
          Loki::TL::TypeAtNonStrict<VertexFormat,0,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,1,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,2,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,3,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,4,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,5,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,6,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,7,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,8,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,9,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,10,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,11,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,12,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,13,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,14,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,15,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,16,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,17,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,18,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,19,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,20,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,21,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,22,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,23,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,24,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,25,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,26,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,27,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,28,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,29,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,30,MVF_empty>::Result::size_of,
          Loki::TL::TypeAtNonStrict<VertexFormat,31,MVF_empty>::Result::size_of
        };
        return (tbl[i]);
      }

      static int offset(int i)
      {
        static const int tbl[]={
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,0>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,1>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,2>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,3>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,4>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,5>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,6>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,7>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,8>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,9>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,10>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,11>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,12>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,13>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,14>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,15>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,16>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,17>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,18>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,19>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,20>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,21>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,22>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,23>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,24>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,25>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,26>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,27>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,28>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,29>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,30>::Result >::value,
          MVF_offset<VertexFormat,typename Loki::TL::TypeAtNonStrict<VertexFormat,31>::Result >::value
        };
        return (tbl[i]);
      }
    };

  } // namespace Mesh
} //namespace LibSL
