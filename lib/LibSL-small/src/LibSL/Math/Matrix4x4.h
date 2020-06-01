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
// LibSL::Math
// ------------------------------------------------------
//
// Matrix 4x4
//
// Layout is row major
// Access with at(column,row)
//
//  WARNING: transpose before using in OpenGL
//
// ------------------------------------------------------
// Sylvain Lefebvre - 2006-06-13
// ------------------------------------------------------

#pragma once

#include <LibSL/Math/Math.h>
#include <LibSL/Math/Tuple.h>
#include <LibSL/CppHelpers/CppHelpers.h>
#include <LibSL/Errors/Errors.h>

namespace LibSL {
  namespace Math {

#define LibSL__Math__M(IJ) const T_Type& m##IJ=(*this)[_##IJ]; // helper to map matrix entries to mij variables

#if defined(_WIN32) || defined(_WIN64)
#define LIBSL_PUBLIC_CLASS __declspec(dllexport)
#else
#define LIBSL_PUBLIC_CLASS /*__attribute__((visibility("default")))*/
#endif

    template <typename T_Type> class Quaternion;

    //! Matrix 4x4 class
    template <typename T_Type>
		class LIBSL_PUBLIC_CLASS Matrix4x4 : public Tuple<T_Type, 16>
    {
    public:

      enum { // column/row
        _00 = 0,
        _10 = 1,
        _20 = 2,
        _30 = 3,
        _01 = 4,
        _11 = 5,
        _21 = 6,
        _31 = 7,
        _02 = 8,
        _12 = 9,
        _22 =10,
        _32 =11,
        _03 =12,
        _13 =13,
        _23 =14,
        _33 =15};

        //! Undefined matrix
        Matrix4x4()
        {
        }

        //! Matrix from table
        explicit Matrix4x4(const T_Type *p)
        {
          ForIndex(n,16) {
            (*this)[n] = p[n];
          }
        }

        //! Matrix from Tuple
        Matrix4x4(const Tuple<T_Type,16>& t)
        {
          ForIndex(n,16) {
            (*this)[n]=t[n];
          }
        }

        //! Matrix from Tuple
        Matrix4x4(const Matrix4x4& t) : Tuple<T_Type,16>(t)
        {
          ForIndex(n,16) {
            (*this)[n]=t[n];
          }
        }

        //! Matrix from quaternion, scaling and translation
        Matrix4x4(const Quaternion<T_Type>& q,const Tuple<T_Type,3>& s,const Tuple<T_Type,3>& t)
        {
          (*this)=q.toMatrix();
          ForIndex(i,3) {
            at(i,i)*=s[i];
          }
          ForIndex(i,3) {
            at(3,i)=t[i];
          }
        }

        //! Matrix from values
        Matrix4x4(
          T_Type m00,T_Type m10,T_Type m20,T_Type m30,
          T_Type m01,T_Type m11,T_Type m21,T_Type m31,
          T_Type m02,T_Type m12,T_Type m22,T_Type m32,
          T_Type m03,T_Type m13,T_Type m23,T_Type m33)
        {
          (*this)[_00]=m00; (*this)[_10]=m10; (*this)[_20]=m20; (*this)[_30]=m30;
          (*this)[_01]=m01; (*this)[_11]=m11; (*this)[_21]=m21; (*this)[_31]=m31;
          (*this)[_02]=m02; (*this)[_12]=m12; (*this)[_22]=m22; (*this)[_32]=m32;
          (*this)[_03]=m03; (*this)[_13]=m13; (*this)[_23]=m23; (*this)[_33]=m33;
        }

        //! Matrix from vectors
        Matrix4x4(Tuple<float,3> p,Tuple<float,3> u,Tuple<float,3> v,Tuple<float,3> w)
        {
          (*this)[_00]=u[0]; (*this)[_10]=v[0]; (*this)[_20]=w[0]; (*this)[_30]=p[0];
          (*this)[_01]=u[1]; (*this)[_11]=v[1]; (*this)[_21]=w[1]; (*this)[_31]=p[1];
          (*this)[_02]=u[2]; (*this)[_12]=v[2]; (*this)[_22]=w[2]; (*this)[_32]=p[2];
          (*this)[_03]=0;    (*this)[_13]=0;    (*this)[_23]=0;    (*this)[_33]=1;
        }

        //! Set to identity matrix
        void eqIdentity()
        {
          Tuple<T_Type,16>::fill(T_Type(0));
          ForIndex(i,4) {
            at(i,i)=T_Type(1);
          }
        }

        static Matrix4x4 identity();

        //! Matrix determinant
        T_Type det() const
        {
          LibSL__Math__M(00); LibSL__Math__M(10); LibSL__Math__M(20); LibSL__Math__M(30);
          LibSL__Math__M(01); LibSL__Math__M(11); LibSL__Math__M(21); LibSL__Math__M(31);
          LibSL__Math__M(02); LibSL__Math__M(12); LibSL__Math__M(22); LibSL__Math__M(32);
          LibSL__Math__M(03); LibSL__Math__M(13); LibSL__Math__M(23); LibSL__Math__M(33);
          return (
            m00*m33*m11*m22 - m00*m23*m11*m32 + m00*m23*m12*m31 + m00*m13*m21*m32
            - m00*m33*m12*m21 - m00*m13*m31*m22 + m23*m02*m30*m11 + m11*m03*m20*m32
            - m11*m33*m02*m20 - m03*m30*m11*m22 - m13*m30*m02*m21 - m13*m01*m20*m32
            - m03*m10*m21*m32 - m23*m12*m01*m30 + m33*m02*m10*m21 + m33*m12*m01*m20
            - m33*m22*m01*m10 + m03*m21*m12*m30 + m13*m31*m02*m20 + m13*m01*m30*m22
            + m23*m32*m01*m10 + m03*m10*m31*m22 - m03*m20*m12*m31 - m23*m02*m10*m31);
        }

        //! Computes inverse matrix
        Matrix4x4 inverse() const
        {
          LibSL__Math__M(00); LibSL__Math__M(10); LibSL__Math__M(20); LibSL__Math__M(30);
          LibSL__Math__M(01); LibSL__Math__M(11); LibSL__Math__M(21); LibSL__Math__M(31);
          LibSL__Math__M(02); LibSL__Math__M(12); LibSL__Math__M(22); LibSL__Math__M(32);
          LibSL__Math__M(03); LibSL__Math__M(13); LibSL__Math__M(23); LibSL__Math__M(33);
          //float d = det();
          T_Type d = det();
          if (LibSL::Math::abs(d) < 1e-18f) {
            throw Errors::Fatal("Matrix4x4::inverse - nul determinant");
          }
          Matrix4x4 i;

          i[_00] =   ( m33*m11*m22 - m23*m11*m32 + m23*m12*m31 + m13*m21*m32 - m33*m12*m21  - m13*m31*m22);
          i[_10] = - ( m13*m20*m32 + m23*m12*m30 - m33*m12*m20 + m33*m22*m10 - m13*m30*m22 - m23*m32*m10);
          i[_20] =   ( m23*m30*m11 - m11*m33*m20 - m13*m30*m21 + m33*m10*m21 + m13*m31*m20 - m23*m10*m31);
          i[_30] = - (-m11*m20*m32 + m30*m11*m22 + m10*m21*m32 - m21*m12*m30 - m10*m31*m22 + m20*m12*m31);

          i[_01] = - ( m33*m22*m01 - m23*m32*m01 + m03*m21*m32 - m33*m02*m21 + m23*m02*m31 - m03*m31*m22);
          i[_11] =   ( m00*m33*m22 - m00*m23*m32 + m23*m02*m30 + m03*m20*m32 - m33*m02*m20 - m03*m30*m22);
          i[_21] = - ( m00*m33*m21 - m00*m23*m31 - m03*m21*m30 + m03*m20*m31 - m33*m01*m20 + m23*m01*m30);
          i[_31] =   ( m00*m21*m32 - m00*m31*m22 - m30*m02*m21 - m01*m20*m32 + m31*m02*m20 + m01*m30*m22);

          i[_02] =   ( m03*m11*m32 - m33*m02*m11 - m13*m01*m32 + m33*m12*m01 + m13*m02*m31 - m03*m31*m12);
          i[_12] = - (-m00*m13*m32 + m00*m33*m12 - m33*m02*m10 + m03*m10*m32 - m03*m12*m30 + m13*m30*m02);
          i[_22] =   ( m33*m00*m11 - m33*m01*m10 - m03*m30*m11 - m13*m00*m31 + m13*m01*m30 + m03*m10*m31);
          i[_32] = - ( m32*m00*m11 - m32*m01*m10 - m02*m30*m11 - m12*m00*m31 + m12*m01*m30 + m02*m10*m31);

          i[_03] = - (-m23*m02*m11 + m03*m11*m22 + m23*m12*m01 - m03*m21*m12 - m13*m01*m22 + m13*m02*m21);
          i[_13] =   ( m00*m23*m12 - m00*m13*m22 + m13*m02*m20 + m03*m10*m22 - m03*m20*m12 - m23*m02*m10);
          i[_23] = - ( m23*m00*m11 - m23*m01*m10 - m03*m20*m11 - m13*m00*m21 + m13*m01*m20 + m03*m10*m21);
          i[_33] =   ( m22*m00*m11 - m22*m01*m10 - m02*m20*m11 - m12*m00*m21 + m12*m01*m20 + m02*m10*m21);

          i = i / d;

          return (i);
        }

        Matrix4x4 transpose() const
        {
          Matrix4x4 r;
          ForIndex(i,4) {
            ForIndex(j,4) {
              r.at(i,j) = at(j,i);
            }
          }
          return (r);
        }

        const T_Type& at(uint c,uint r) const
        {
          return (*this)[c+r*4];
        }

        T_Type& at(uint c,uint r)
        {
          return (*this)[c+r*4];
        }

        Matrix4x4 operator *(const Matrix4x4& m) const
        {
          Matrix4x4 r;
          r.fill(T_Type(0));
          ForIndex (i,4) {
            ForIndex (j,4) {
              ForIndex (k,4) {
                r.at(j,i)+=at(k,i)*m.at(j,k);
              }
            }
          }
          return (r);
        }

        Tuple<T_Type,4> operator *(const Tuple<T_Type,4>& v) const
        {
          Tuple<T_Type,4> r;
          r.fill(T_Type(0));
          ForIndex (i,4) {
            ForIndex (k,4) {
              r[i]+=at(k,i)*v[k];
            }
          }
          return (r);
        }

        Tuple<T_Type,4> mul(const Tuple<T_Type,4>& v) const
        {
          Tuple<T_Type,4> r;
          r.fill(T_Type(0));
          ForIndex (i,4) {
            ForIndex (k,4) {
              r[i]+=at(k,i)*v[k];
            }
          }
          return (r);
        }

        Tuple<T_Type,3> mulPoint(const Tuple<T_Type,3>& p) const
        {
          Tuple<T_Type,3> r;
          r.fill(T_Type(0));
          ForIndex (i,3) {
            ForIndex (k,3) {
              r[i]+=at(k,i)*p[k];
            }
            r[i]+=at(3,i);
          }
          return (r);
        }

        Tuple<T_Type,3> mulVector(const Tuple<T_Type,3>& v) const
        {
          Tuple<T_Type,3> r;
          r.fill(T_Type(0));
          ForIndex (i,3) {
            ForIndex (k,3) {
              r[i]+=at(k,i)*v[k];
            }
          }
          return (r);
        }

    };


    // static methods

    template< typename T_Type >
    Matrix4x4<T_Type> Matrix4x4<T_Type>::identity()
    {
      Matrix4x4<T_Type> m;
      m.eqIdentity();
      return m;
    }

    // typedefs

    typedef Matrix4x4<float>  m4x4f;
    typedef Matrix4x4<double> m4x4d;

    // translation matrix
    template <typename T_Type>
    Matrix4x4<T_Type> translationMatrix(const Tuple<T_Type,3>& t)
    {
      Matrix4x4<T_Type> m; m.eqIdentity();
      ForIndex(i,3) {
        m.at(3,i)=t[i];
      }
      return (m);
    }

    // scale matrix
    template <typename T_Type>
    Matrix4x4<T_Type> scaleMatrix(const Tuple<T_Type,3>& s)
    {
      Matrix4x4<T_Type> m; m.eqIdentity();
      ForIndex(i,3) {
        m.at(i,i)=s[i];
      }
      return (m);
    }

    // remove scaling from matrix
    template <typename T_Type>
    Matrix4x4<T_Type> removeScaling(const Matrix4x4<T_Type>& s)
    {
      Tuple<T_Type,3> v[3];
      ForIndex(n,3) {
        ForIndex(i,3) {
          v[n][i] = s.at(i,n);
        }
        v[n] = normalize_safe(v[n]);
      }
      Matrix4x4<T_Type> m; m = s;
      ForIndex(n,3) {
        ForIndex(i,3) {
           m.at(i,n) = v[n][i];
        }
      }
      return (m);
    }

    // mirror matrix
    template <typename T_Type>
    Matrix4x4<T_Type> mirrorMatrix(const Tuple<T_Type,3>& n)
    {
      return Matrix4x4<T_Type>(
        1-2*n[0]*n[0],      n[0]*n[1],          n[0]*n[2],         0,
        n[1]*n[0],      1-2*n[1]*n[1],          n[1]*n[2],         0,
        n[2]*n[0],          n[2]*n[1],      1-2*n[2]*n[2],         0,
        0,                  0,                  0,                 1);
    }

    // lookat matrix (RH)
    template <typename T_Type>
    Matrix4x4<T_Type> lookatMatrix(
      const Tuple<T_Type,3>& eye,
      const Tuple<T_Type,3>& at,
      const Tuple<T_Type,3>& up)
    {
      Tuple<T_Type,3> zaxis = normalize(eye - at); // revert for LH
      Tuple<T_Type,3> xaxis = normalize(cross(up,zaxis));
      Tuple<T_Type,3> yaxis = cross(zaxis, xaxis);
      return Matrix4x4<T_Type>(
        xaxis[0],           xaxis[1],           xaxis[2],          -dot(xaxis, eye),
        yaxis[0],           yaxis[1],           yaxis[2],          -dot(yaxis, eye),
        zaxis[0],           zaxis[1],           zaxis[2],          -dot(zaxis, eye),
        0,                  0,                  0,                  1);
    }

    // perspective matrix from frustum (RH)
    // maps (x,y,z) to [-1,1]x[-1,1]x[**0**,1] (Direct3D convention)
    template <typename T_Type>
    Matrix4x4<T_Type> perspectiveMatrixD3D(
      T_Type l,T_Type r,
      T_Type b,T_Type t,
      T_Type zn,T_Type zf)
    {
      return Matrix4x4<T_Type>(
        2*zn/(r-l),          0, (l+r)/(r-l),             0,
        0,          2*zn/(t-b), (t+b)/(t-b),             0,
        0,                   0,  zf/(zn-zf), zn*zf/(zn-zf),
        0,                   0,           -1,            0);
    }

    // perspective matrix from FOV (RH)
    // maps (x,y,z) to [-1,1]x[-1,1]x[**0**,1] (Direct3D convention)
    template <typename T_Type>
    Matrix4x4<T_Type> perspectiveMatrixD3D(
      T_Type fov,T_Type aspect,
      T_Type zn, T_Type zf)
    {
      T_Type yScale = LibSL::Math::cot(fov/2.0f);
      T_Type xScale = yScale / aspect;
      return Matrix4x4<T_Type>(
        xScale,    0,          0,             0,
        0,    yScale,          0,             0,
        0,         0, zf/(zn-zf), zn*zf/(zn-zf),
        0,         0,         -1,             0);
    }


    // perspective matrix from frustum (RH)
    // maps (x,y,z) to [-1,1]x[-1,1]x[**-1**,1] (OpenGL convention)
    template <typename T_Type>
    Matrix4x4<T_Type> perspectiveMatrixGL(
      T_Type l,T_Type r,
      T_Type b,T_Type t,
      T_Type zn,T_Type zf)
    {
      return Matrix4x4<T_Type>(
        2*zn/(r-l),          0, (l+r)/(r-l),             0,
        0,          2*zn/(t-b), (t+b)/(t-b),             0,
        0,                   0, (zn+zf)/(zn-zf), 2*zn*zf/(zn-zf),
        0,                   0,           -1,            0);
    }

    // perspective matrix from FOV (RH)
    // maps (x,y,z) to [-1,1]x[-1,1]x[**-1**,1] (OpenGL convention)
    template <typename T_Type>
    Matrix4x4<T_Type> perspectiveMatrixGL(
      T_Type fov,T_Type aspect,
      T_Type zn, T_Type zf)
    {
      T_Type yScale = LibSL::Math::cot(fov/2.0f);
      T_Type xScale = yScale / aspect;
      return Matrix4x4<T_Type>(
        xScale,    0,          0,             0,
        0,    yScale,          0,             0,
        0,         0, (zn+zf)/(zn-zf), 2*zn*zf/(zn-zf),
        0,         0,         -1,             0);
    }

    // ortho matrix
    // maps (x,y,z) to [-1,1]x[-1,1]x[**0**,1] (Direct3D convention)
    template <typename T_Type>
    Matrix4x4<T_Type> orthoMatrixD3D(
      T_Type l,T_Type r,
      T_Type b,T_Type t,
      T_Type zn,T_Type zf)
    {
      float tx = - (r+l) / (r-l);
      float ty = - (t+b) / (t-b);
      float tz =   zn / (zn-zf);
      return Matrix4x4<T_Type>(
        2.0f / (r-l),            0,             0,          tx,
                   0, 2.0f / (t-b),             0,          ty,
                   0,            0,  1.0f/(zn-zf),          tz,
                   0,            0,             0,           1);
    }

    // ortho matrix
    // maps (x,y,z) to [-1,1]x[-1,1]x[**-1**,1] (OpenGL convention)
    template <typename T_Type>
    Matrix4x4<T_Type> orthoMatrixGL(
      T_Type l,T_Type r,
      T_Type b,T_Type t,
      T_Type zn,T_Type zf)
    {
      float tx = - (r+l) / (r-l);
      float ty = - (t+b) / (t-b);
      float tz = - (zf+zn) / (zf-zn);
      return Matrix4x4<T_Type>(
        2.0f / (r-l),            0,             0,          tx,
                   0, 2.0f / (t-b),             0,          ty,
                   0,            0, -2.0f/(zf-zn),          tz,
                   0,            0,             0,           1);
    }

    // pick matrix
    template <typename T_Type>
    Matrix4x4<T_Type> pickMatrix(
      T_Type x, T_Type y,
      T_Type w, T_Type h,
      T_Type viewport[4])
    {
      float sx, sy;
      float tx, ty;
      sx = viewport[2] / w;
      sy = viewport[3] / h;
      tx = (viewport[2] + 2.0 * (viewport[0] - x)) / w;
      ty = (viewport[3] + 2.0 * (viewport[1] - y)) / h;
      return Matrix4x4<T_Type>(
        sx, 0, 0, tx,
        0, sy, 0, ty,
        0,  0, 1,  0,
        0,  0, 0,  1);
    }

    // print matrix
    template <typename T_Type>
    std::ostream& operator<<(std::ostream& s,const Matrix4x4<T_Type>& m)
    {
      ForIndex(r,4) {
        s << '[';
        ForIndex(c,4) {
          s << m.at(c,r);
          if (c < 3) s << ',';
        }
        s << ']';
        if (r < 3) s << std::endl;
      }
      return (s);
    }

#ifdef OPENGL
template <typename T_Type>
Matrix4x4<T_Type> perspectiveMatrix(T_Type fov,T_Type aspect, T_Type zn, T_Type zf)          { return perspectiveMatrixGL(fov,aspect,zn,zf); }
template <typename T_Type>
Matrix4x4<T_Type> perspectiveMatrix(T_Type l,T_Type r,T_Type b,T_Type t,T_Type zn,T_Type zf) { return perspectiveMatrixGL(l,r,b,t,zn,zf); }
template <typename T_Type>
Matrix4x4<T_Type> orthoMatrix(T_Type l,T_Type r,T_Type b,T_Type t,T_Type zn,T_Type zf)       { return orthoMatrixGL(l,r,b,t,zn,zf); }
#endif

#ifdef DIRECT3D
template <typename T_Type>
Matrix4x4<T_Type> perspectiveMatrix(T_Type fov,T_Type aspect, T_Type zn, T_Type zf)          { return perspectiveMatrixD3D(fov,aspect,zn,zf); }
template <typename T_Type>
Matrix4x4<T_Type> perspectiveMatrix(T_Type l,T_Type r,T_Type b,T_Type t,T_Type zn,T_Type zf) { return perspectiveMatrixD3D(l,r,b,t,zn,zf); }
template <typename T_Type>
Matrix4x4<T_Type> orthoMatrix(T_Type l,T_Type r,T_Type b,T_Type t,T_Type zn,T_Type zf)       { return orthoMatrixD3D(l,r,b,t,zn,zf); }
#endif

  } //namespace LibSL::Math
} //namespace LibSL

