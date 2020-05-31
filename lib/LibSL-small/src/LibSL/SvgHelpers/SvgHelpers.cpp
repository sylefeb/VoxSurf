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
#include "LibSL.precompiled.h"
// ------------------------------------------------------

#include "SvgHelpers.h"
#include <LibSL/Errors/Errors.h>
#include <LibSL/Math/Matrix4x4.h>
#include <LibSL/Math/Vertex.h>
#include <LibSL/Geometry/AAB.h>

using namespace LibSL;
using namespace LibSL::Math;

// ------------------------------------------------------

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

// ------------------------------------------------------

#define NAMESPACE LibSL::SvgHelpers

// ------------------------------------------------------

NAMESPACE::Svg::Svg(
  const std::string &fname, 
  LibSL::Geometry::AAB<2> viewbox,
  m4x4f transform)
  : m_File(fname.c_str())
{
  m_Trsf = transform;
  m_StrokeColor = "#000000";
  m_FillColor = "#ffffff";
  m_StrokeWidth = 1.0f;
  m_File << "\
<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n\
<!-- Created with LibSL -->\n\
<svg\n\
   xmlns=\"http://www.w3.org/2000/svg\"\n\
   baseProfile=\"full\"\n\
   version=\"1.1\"\n";
   if (!viewbox.empty()) {
     m_File << "viewBox=\""
       << viewbox.minCorner()[0] << " "
       << viewbox.minCorner()[1] << " "
       << viewbox.extent()[0] << " "
       << viewbox.extent()[1]
       << "\"";
   }
  m_File << ">\n\
  <g id=\"layer1\">\n\
    ";
}

// ------------------------------------------------------

void NAMESPACE::Svg::setProperties(std::string strokeColor, std::string fillColor, float strokeWidth)
{
	m_StrokeColor = strokeColor;
	m_FillColor   = fillColor;
	m_StrokeWidth = strokeWidth;
}

// ------------------------------------------------------

void NAMESPACE::Svg::startPath()
{
  m_File << "<path \n";
  m_File << "style=\"fill:"<<m_FillColor<<";fill-rule:nonzero;stroke:" << m_StrokeColor << ";stroke-width:" << m_StrokeWidth << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n";
  m_File << "d=\"M ";
}

// ------------------------------------------------------

void NAMESPACE::Svg::addPoint(float x,float y)
{
  v3f t = m_Trsf.mulPoint(v3f(x, y, 0.0f));
  m_File << t[0] << ',' << t[1] << ' ';
}

// ------------------------------------------------------

void NAMESPACE::Svg::endPath(bool open)
{
  if (!open) {
    m_File << "Z ";
  }
  m_File << "\"\n";
  m_File << "/>\n";
}

// ------------------------------------------------------

void NAMESPACE::Svg::startPolygon()
{
  m_File << "<path \n";
  m_File << "style=\"fill:" << m_FillColor << ";fill-rule:nonzero;stroke:" << m_StrokeColor << ";stroke-width:" << m_StrokeWidth << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n";
  m_File << "d=\"";
}

void NAMESPACE::Svg::startPolygonPath()
{
  m_File << "M ";
}

// ------------------------------------------------------

void NAMESPACE::Svg::endPolygonPath()
{
  m_File << "Z ";
}

// ------------------------------------------------------

void NAMESPACE::Svg::endPolygon()
{
  m_File << "\"\n";
  m_File << "/>\n";
}

// ------------------------------------------------------

void NAMESPACE::Svg::addCircle(float x, float y,float r)
{
  v3f t = m_Trsf.mulPoint(v3f(x, y, 0.0f));
  m_File << "<circle \
    style=\"fill:"<<m_FillColor<<";fill-rule:evenodd;stroke:"<<m_StrokeColor<<";stroke-width:"<<m_StrokeWidth<<"px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n\
    cx=\"" << t[0] << "\"\n\
    cy=\"" << t[1] << "\"\n\
    r=\"" << r << "\" />\n";
}

void NAMESPACE::Svg::addText(float x, float y, const char *txt)
{
  v3f t = m_Trsf.mulPoint(v3f(x, y, 0.0f));
  m_File << "<text \
    style=\"fill:"<<m_FillColor<<";\"\n\
    x=\"" << t[0] << "\"\n\
    y=\"" << t[1] << "\"\n\
        >" << txt << "</text>";
}

// ------------------------------------------------------

NAMESPACE::Svg::~Svg()
{
  m_File << "\
</g>\n\
</svg>\n\
";
  //m_File.close();
}

// ------------------------------------------------------
