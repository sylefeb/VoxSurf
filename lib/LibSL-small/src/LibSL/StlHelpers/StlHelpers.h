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
// LibSL::StlHelpers
// ------------------------------------------------------
//
// STL syntax helpers
// 
// ------------------------------------------------------
// Sylvain Lefebvre - 2006-02-23
// ------------------------------------------------------

#pragma once

#include <LibSL/LibSL.common.h>

#include <string>

namespace LibSL {
	namespace StlHelpers {

    LIBSL_DLL std::string loadFileIntoString         (const char *);
    LIBSL_DLL std::string extractFileName            (const std::string& path);
    LIBSL_DLL std::string extractExtension           (const std::string& path);
    LIBSL_DLL std::string extractPath                (const std::string& path);
    LIBSL_DLL std::string removeExtensionFromFileName(const std::string& fname);
    LIBSL_DLL std::string replaceBy                  (const std::string& match,const std::string& replacement,const std::string& str);

    /// fixes transform + tolower issue (defined both as a one arg and two arg function)
    template<class charT> charT toLower(charT c) {
	return tolower(c); // explicitely call one argument version of tolower
    }
  } //namespace LibSL::StlHelpers
} //namespace LibSL

// ------------------------------------------------------

#define ForIterator(T,C,I)      for (T::iterator I=C.begin();I!=C.end();I++)
#define ForConstIterator(T,C,I) for (T::const_iterator I=C.begin();I!=C.end();I++)

#define ForIterator_typename(T,C,I)      for (typename T::iterator I=C.begin();I!=C.end();I++)
#define ForConstIterator_typename(T,C,I) for (typename T::const_iterator I=C.begin();I!=C.end();I++)

// ------------------------------------------------------
