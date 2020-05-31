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
// LibSL::BasicParser
// ------------------------------------------------------
//
// Basic parser
// 
// ------------------------------------------------------
// Sylvain Lefebvre - 2008-05-10
// ------------------------------------------------------
/*

NOTES

2012-11-06 WARNING  Changed all methods to keep the stopper by default ('keep_stopper=true').
           This may lead to trouble in prior parsing codes.

2012-11-06 TODO/FIXME readInt / readFloat do not have the same behaviour regarding stoppers.
           (readInt stops if a non-digit char is found, which is better and required by some parsers)

*/
#pragma once

#include <LibSL/Memory/Array.h>
#include <LibSL/System/Types.h>
#include <LibSL/CppHelpers/CppHelpers.h>

#include <fstream>
#include <cstdio>
#include <cstring>

#ifdef _MSC_VER
#pragma warning ( disable : 4996 )
#endif

// ------------------------------------------------------

#define IS_SPACE(c) ((c) == ' ' || (c) == '\t' || (c) == '\r')
#define IS_EOL(c)   ((c) == '\n')

// ------------------------------------------------------

namespace LibSL {
  namespace BasicParser {

    //! file stream for basic parser
    class FileStream
    {
    protected:

      FILE *m_File;
      bool  m_EOF;

    public:

      FileStream(const char *fname) 
      {
        m_EOF  = false;
        m_File = NULL;
				fopen_s(&m_File, fname, "rb");
        if (m_File == NULL) {
          throw LibSL::Errors::Fatal("LibSL::CppHelpers::FileStream - cannot open '%s'",fname);
        }
      }

      void close()
      {
        fclose(m_File);
      }

      int getc()
      {
        int c = fgetc(m_File);
        if (c == EOF) {
          m_EOF = true;
        }
        return c;
      }

      void ungetc(int c)
      {
        ::ungetc(c,m_File);
        if (c == EOF) {
          m_EOF = false;
        }
      }

      bool eof() const { return (m_EOF); }
    };

    //! buffer stream for basic parser
    class BufferStream
    {
    protected:
      
      const char *m_Buffer;
      uint        m_Size;
      int         m_Pos;
      bool        m_EOF;
    
    public:

      BufferStream(const char *buffer,uint size) 
      {
        sl_assert(buffer != NULL);
        sl_assert(size > 0);
        m_Buffer = buffer;
        m_Size   = size;
        m_Pos    = 0;
        m_EOF    = false;
      }

      void close()
      {
        
      }

      int getc()
      {
        if (m_Pos >= int(m_Size)) {
          m_EOF = true;
          return EOF;
        }
        char c = m_Buffer[m_Pos];
        m_Pos  = m_Pos + 1;
        return int(c);
      }

      void ungetc(int c)
      {
        if (m_Pos > 0) {
          m_Pos = m_Pos - 1;
        }
        if (c == EOF) {
          m_EOF = false;
        }
      }

      bool eof() const { return (m_EOF); }
    };

    //! basic parser
    template<class T_Stream>
    class Parser
    {
    protected:

      bool             m_EOLAsSpace;
      T_Stream&        m_Stream;
	  
	  static const int StringBufferSize = 1024;
	  LibSL::Memory::Array::Array<char> m_StringBuffer;

    public:

      Parser(T_Stream& stream,bool eol_as_space = true)
        : m_Stream(stream), m_StringBuffer(StringBufferSize)
      {
        m_EOLAsSpace = eol_as_space;
      }

      ~Parser()
      {
        m_Stream.close();
      }

      bool        reachChar(char toBeReached)
      {
        if (m_Stream.eof()) {
          return (false);
        }
        int c;
        do {
          c = m_Stream.getc();
          if (m_Stream.eof()) {
            return (false);
          }
        } while (char(c) != toBeReached);
        return (true);
      }

      int         reachOneOf(const char *delims)
      {
        if (m_Stream.eof()) {
          return (EOF);
        }
        int c;
        do {
          c = m_Stream.getc();
          if (m_Stream.eof()) {
            return (EOF);
          }
        } while (strchr(delims,char(c)) == NULL);
        return (c);
      }

      bool        skipSpaces()
      {
        if (m_Stream.eof()) {
          return (false);
        }
        int c;
        do {
          c = m_Stream.getc();
          if (c == EOF) {
            return (false);
          }
        } while (IS_SPACE(c) || (m_EOLAsSpace && IS_EOL(c)));
        m_Stream.ungetc(c);
        return (true);
      }

      int         readChar(bool advance = true)
      {
        if (!skipSpaces()) {
          return (EOF);
        }
        int c = m_Stream.getc();
        if (!advance) {
          m_Stream.ungetc(c);
        }
        return (c);
      }

      char *readString(const char *eos = NULL, const char *accepted = NULL)
      {
        int bufpos = 0;
        m_StringBuffer[bufpos] = '\0';
        int         c;
        if (!skipSpaces()) {
          return (m_StringBuffer.raw());
        }
        do {
          sl_assert(bufpos < StringBufferSize);
          c = m_Stream.getc();
          // end here ?
          if (m_Stream.eof()) {
            // yes: eof
            m_StringBuffer[bufpos++] = '\0';
            return (m_StringBuffer.raw());
          } else {
            if (eos == NULL) {
              if (accepted == NULL) {
                // yes: space
                if (IS_SPACE(c) || IS_EOL(c)) {
                  m_Stream.ungetc(c);
                  m_StringBuffer[bufpos++] = '\0';
                  return (m_StringBuffer.raw());
                }
              } else {
                if (strchr(accepted, c) == NULL) {
                  // not one of the accepted char
                  m_Stream.ungetc(c);
                  m_StringBuffer[bufpos++] = '\0';
                  return (m_StringBuffer.raw());
                }
              }
            } else {
              // yes: user given stopper
              if (strchr(eos,c) != NULL) {
                m_Stream.ungetc(c);
                m_StringBuffer[bufpos++] = '\0';
                return (m_StringBuffer.raw());
              }
            }
          }
          // no: append
          m_StringBuffer[bufpos++] = c;
        } while ( ! m_Stream.eof() );
        sl_assert(false); 
        return (NULL);
      }

      char *readUntil(const char *eos)
      {
        int bufpos = 0;
        m_StringBuffer[bufpos] = '\0';
        int         c;
        do {
          sl_assert(bufpos < StringBufferSize);
          c = m_Stream.getc();
          // end here ?
          if (m_Stream.eof()) {
            // yes: eof
            m_StringBuffer[bufpos++] = '\0';
            return (m_StringBuffer.raw());
          } else {
            // yes: user given stopper
            if (strchr(eos,c) != NULL) {
              m_Stream.ungetc(c);
              m_StringBuffer[bufpos++] = '\0';
              return (m_StringBuffer.raw());
            }
          }
          // no: append
          m_StringBuffer[bufpos++] = c;
        } while ( ! m_Stream.eof() );
        sl_assert(false); 
        return (NULL);
      }

      bool reachString(const char *str,const char *eos=" ")
      {
        while (!m_Stream.eof()) {
          const char *read = readString(eos);
          if (!strcmp(read,str)) {
            return (true);
          }
        }
        return (false);
      }

      float       readFloat()
      {
        const char *str = readString(NULL,"0123456789.-+e");
        return float(atof(str));
      }

      int         readInt()
      {
        bool first = true;
        int  n     = 0;
        int  s     = 1;
        int  c;
        if (!skipSpaces()) {
          //return (0);
		      sl_assert(false);
        }
        do {
          c = m_Stream.getc();
          // end here ?
          if (m_Stream.eof()) {
            // yes: eof
            return (s*n);
          } else if (first && (c == '-' || c == '+')) {
            if (c == '-') {
              s = -1;
            }
          } else if (IS_SPACE(c) || IS_EOL(c)) {
		      	m_Stream.ungetc(c);
            return (s*n);
          } else if (!isdigit(c)) {
            m_Stream.ungetc(c);
            return (s*n);
          } else {
            // no: add
            sl_assert(isdigit(c) != 0);
            n     = n * 10 + int(c - '0');
          }
          first = false;
        } while (!m_Stream.eof());
        sl_assert(false);
        return (0);
      }

      bool eof() const { return m_Stream.eof(); }

      char *trim(char *str,const char *charList)
      {
		char *start = str;
        int len = (int)strlen(str);
		// left trim
		ForIndex(i,len) {
			if (strchr(charList,str[i]) != NULL ) {
				start ++;
			} else {
				break;
			}
		}
		// right trim
		int i = len - 1;
		while (i >= 0) {
			if (strchr(charList,str[i]) != NULL ) {
				str[i] = '\0';
			}
			i --;
		}
        return start;
      }

    };

  } //namespace LibSL::BasicParser
} //namespace LibSL

// ------------------------------------------------------
