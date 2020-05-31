.. Copyright (c) 2017, Wolf Vollprecht, Johan Mabille and Sylvain Corlay

   Distributed under the terms of the BSD 3-Clause License.

   The full license is in the file LICENSE, distributed with this software.

.. image:: xtensor-blas.svg
   :alt: xtensor-blas

xtensor bindings for BLAS and LAPACK

Introduction
------------

``xtensor-blas`` is an extension to the xtensor library, offering bindings to BLAS and LAPACK libraries through cxxblas and cxxlapack from the FLENS_ project.

``xtensor-blas`` currently provides non-broadcasting ``dot``, ``norm`` (1- and 2-norm for vectors), ``inverse``, ``solve``, ``eig``, ``cross``, ``det``, ``slogdet``, ``matrix_rank``, ``inv``, ``cholesky``, ``qr``, ``svd`` in the ``xt::linalg`` namespace (check the corresponding ``xlinalg.hpp`` header for the function signatures). The functions, and signatures, are trying to be 1-to-1 equivalent to NumPy. Low-level functions to interface with BLAS or LAPACK with xtensor containers are also offered in the blas and lapack namespace.

``xtensor`` and ``xtensor-blas`` require a modern C++ compiler supporting C++14. The following C++ compilers are supported:

- On Windows platforms, Visual C++ 2015 Update 2, or more recent
- On Unix platforms, gcc 4.9 or a recent version of Clang

Licensing
---------

We use a shared copyright model that enables all contributors to maintain the
copyright on their contributions.

This software is licensed under the BSD-3-Clause license. See the LICENSE file for details.

.. toctree::
   :caption: INSTALLATION
   :maxdepth: 2

   installation
   performance

.. toctree::
   :caption: USAGE
   :maxdepth: 2

   usage
   reference

.. toctree::
   :caption: DEVELOPER ZONE

   build-options
   releasing

.. _FLENS: https://github.com/michael-lehn/FLENS
