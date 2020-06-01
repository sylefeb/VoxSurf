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

#include "StlHelpers.h"
#include <LibSL/Errors/Errors.h>

// ------------------------------------------------------

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

// ------------------------------------------------------

#define NAMESPACE LibSL::StlHelpers

// ------------------------------------------------------

string NAMESPACE::loadFileIntoString(const char *file)
{
  ifstream infile(file, ios::binary );
  if (!infile) {
    throw LibSL::Errors::Fatal("[loadFileIntoString] - file '%s' not found",file);
  }
  ostringstream strstream;
  while (infile) { // TODO: improve efficienty
    ifstream::int_type c = infile.get();
    if (c != EOF)
      strstream << char(c);
    else
      break;
  } 
  return strstream.str();
}

// ------------------------------------------------------

std::string NAMESPACE::extractFileName(const std::string& path)
{
  // search for last '\\' or '/'
  size_t pos;
  pos = path.rfind("\\");
  if (pos == string::npos) {
    pos   = path.rfind("/");
    if (pos == string::npos) {
      return path;
    }
  }
  string fname = path.substr(pos+1,path.length()-pos);
  return fname;
}

// ------------------------------------------------------

std::string NAMESPACE::extractExtension(const std::string& path)
{
  // search for last '.'
  size_t pos   = path.rfind(".");
  if (pos == string::npos) {
    return path;
  }
  string newfname = path.substr(pos+1,path.length()-pos);
  return newfname;
}

// ------------------------------------------------------

std::string NAMESPACE::extractPath(const std::string& path)
{
  // search for last '\\' or '/'
  size_t pos;
  pos = path.rfind("\\");
  if (pos == string::npos) {
    pos   = path.rfind("/");
    if (pos == string::npos) {
      return path;
    }
  }
  string dname = path.substr(0,pos);
  return dname;
}

// ------------------------------------------------------

std::string NAMESPACE::removeExtensionFromFileName(const std::string& fname)
{
  // search for last '.'
  size_t pos   = fname.rfind(".");
  if (pos == string::npos) {
    return fname;
  }
  string newfname = fname.substr(0,pos);
  return newfname;
}

// ------------------------------------------------------

std::string NAMESPACE::replaceBy(const std::string& match,const std::string& replacement,const std::string& str)
{
  std::string tmp = str;
  while (1) {
    size_t pos = tmp.find(match);
    if (pos == string::npos) {
      return tmp;
    }
    // replace
    tmp.replace(pos,match.length(),replacement);
  }
  return tmp;
}

// ------------------------------------------------------
