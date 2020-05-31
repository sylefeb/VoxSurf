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
#include "LibSL.precompiled.h"
// -------------------------------------------------

#include <LibSL/Errors/Errors.h>
#include <LibSL/Math/Math.h>
#include <LibSL/Memory/Array.h>

#include "VertexFormat.h"

#include <sstream>
#include <cstring>

using namespace std;

// -------------------------------------------------

#define NAMESPACE LibSL::Mesh

using namespace LibSL::Errors;
using namespace LibSL::Math;
using namespace LibSL::Memory::Array;

using namespace NAMESPACE;

// -------------------------------------------------

void NAMESPACE::MVF::addAttribute(const Attribute& a)
{
  if (!a.isValid()) {
    throw Fatal("MVF::addAttribute - adding a non-valid attribute!");
  }
  if (findAttributeByBinding(a.binding) != NULL) {
    stringstream sstr;
    sstr << a.binding;
    throw Fatal("MVF::addAttribute - binding '%s' already defined in MVF",sstr.str().c_str());
  }
  m_AttributesByBinding.insert( make_pair(a.binding,(int)m_Attributes.size()) );
  m_Attributes         .push_back(a);
}

// -------------------------------------------------

const MVF::Attribute *NAMESPACE::MVF::findAttributeByBinding(e_Binding b) const
{
  std::map<e_Binding,int>::const_iterator A = m_AttributesByBinding.find(b);
  if (A == m_AttributesByBinding.end()) {
    return NULL;
  } else {
    return &(m_Attributes[A->second]);
  }
}

// -------------------------------------------------

void NAMESPACE::MVF::convertVertexFrom(void *dst,const void *src,const MVF *srcmvf) const
{
  sl_assert(dst    != NULL);
  sl_assert(src    != NULL);
  sl_assert(srcmvf != NULL);
  ForIndex(i,m_Attributes.size()) {
    const Attribute *a = srcmvf->findAttributeByBinding( m_Attributes[i].binding );
    if (a != NULL) {
      // attribute is present in source mvf
      // -> copy data
      memcpy((uchar*)dst + m_Attributes[i].offset,(const uchar *)src + a->offset,min(m_Attributes[i].size_of,a->size_of) );
    }
  }
}

// -------------------------------------------------

v3f *NAMESPACE::MVF::pos3(void *data)
{
  const Attribute *a = findAttributeByBinding( Position );
  if (a == NULL) {
    throw Fatal("MVF::pos3 - no position in MVF");
  }
  if (a->numComponents != 3) {
    throw Fatal("MVF::pos3 - position has %d components, not 3",a->numComponents);
  }
  if (a->type != Float) {
    throw Fatal("MVF::pos3 - position type is not float");
  }
  return (v3f*)((const uchar *)data + a->offset);
}

// -------------------------------------------------

const v3f *NAMESPACE::MVF::pos3(const void *data) const
{
  const Attribute *a = findAttributeByBinding( Position );
  if (a == NULL) {
    throw Fatal("MVF::pos3 - no position in MVF");
  }
  if (a->numComponents != 3) {
    throw Fatal("MVF::pos3 - position has %d components, not 3",a->numComponents);
  }
  if (a->type != Float) {
    throw Fatal("MVF::pos3 - position type is not float");
  }
  return (v3f*)((const uchar *)data + a->offset);
}

// -------------------------------------------------

void *NAMESPACE::MVF::attr(void *data,MVF::e_Type type,MVF::e_Binding binding,uint nComp)
{
	const Attribute *a = findAttributeByBinding( binding );
	if (a == NULL) {
		throw Fatal("MVF::attr - attribute not found in MVF");
	}
    sl_assert(a->binding == binding);
	if (a->numComponents != nComp) {
		throw Fatal("MVF::attr - attribute has incorrect number of components: %i",a->numComponents);
	}
	if (a->type != type) {
		throw Fatal("MVF::attr - attribute type mismatch");
	}
	return (void*)((uchar *)data + a->offset);
}

// -------------------------------------------------

bool NAMESPACE::MVF::hasAttr(MVF::e_Type type,MVF::e_Binding binding,uint nComp)
{
	const Attribute *a = findAttributeByBinding( binding );
	if (a == NULL) {
		return false;
	}
	sl_assert(a->binding == binding);
	if (a->numComponents != nComp) {
		return false;
	}
	if (a->type != type) {
		return false;
	}
	return true;
}

// -------------------------------------------------

bool NAMESPACE::MVF::load(FILE *f)
{
  int n  = 0;
  int nr = (int)fread(&n,sizeof(int),1,f);
  if (nr == 0) {
    return false;
  }
  m_Attributes         .clear();
  m_AttributesByBinding.clear();
  Array<Attribute> attribs;
  attribs.allocate(n);
  nr = (int)fread(&(attribs[0]),sizeof(Attribute),n,f);
  ForIndex(a,attribs.size()) {
    addAttribute(attribs[a]);
  }
  return true;
}

// -------------------------------------------------

void NAMESPACE::MVF::save(FILE *f) const
{
  int n = (int)m_Attributes.size();
  fwrite(&n,sizeof(int),1,f);
  fwrite(&(m_Attributes[0]),sizeof(Attribute),n,f);
}

// -------------------------------------------------

std::ostream& operator<<(std::ostream& s,const LibSL::Mesh::MVF::e_Type& t)
{
  switch (t)
  {
  case MVF::Byte:   s << "Byte";   break;
  case MVF::Float:  s << "Float";  break;
  case MVF::Double: s << "Double"; break;
  case MVF::Int:    s << "Int";    break;
  default: sl_assert(false);
  }
  return s;
}

// -------------------------------------------------

std::ostream& operator<<(std::ostream& s,const LibSL::Mesh::MVF::e_Binding& b)
{
  switch (b)
  {
  case MVF::Position:  s << "Position";  break;
  case MVF::Normal:    s << "Normal";    break;
  case MVF::Color0:    s << "Color0";    break;
  case MVF::Color1:    s << "Color1";    break;
  case MVF::TexCoord0: s << "TexCoord0"; break;
  case MVF::TexCoord1: s << "TexCoord1"; break;
  case MVF::TexCoord2: s << "TexCoord2"; break;
  case MVF::TexCoord3: s << "TexCoord3"; break;
  case MVF::TexCoord4: s << "TexCoord4"; break;
  case MVF::TexCoord5: s << "TexCoord5"; break;
  case MVF::TexCoord6: s << "TexCoord6"; break;
  case MVF::TexCoord7: s << "TexCoord7"; break;
  default: sl_assert(false);
  }
  return s;
}

// -------------------------------------------------

std::ostream& operator<<(std::ostream& s,const LibSL::Mesh::MVF::Attribute& a)
{
  s << '[' << a.index
    << ": "
    << a.binding
    << ", " << int(a.numComponents)
    << 'x'
    << a.type
    << " szof: " << a.size_of
    << " addr: " << a.offset;
  return s;
}

// -------------------------------------------------

std::ostream& operator<<(std::ostream& s,const LibSL::Mesh::MVF& mvf)
{
  ForIndex(a,mvf.attributes().size()) {
    s << mvf.attributes()[a] << std::endl;
  }
  return s;
}

// -------------------------------------------------
