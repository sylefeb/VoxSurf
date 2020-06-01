.. Copyright (c) 2017, Johan Mabille and Sylvain Corlay

   Distributed under the terms of the BSD 3-Clause License.

   The full license is in the file LICENSE, distributed with this software.


.. raw:: html

   <style>
   .rst-content .section>img {
       width: 30px;
       margin-bottom: 0;
       margin-top: 0;
       margin-right: 15px;
       margin-left: 15px;
       float: left;
   }
   </style>

Installation
============

Although ``xtl`` is a header-only library, we provide standardized means to install it, with package managers or with cmake.

Besides the xtl headers, all these methods place the ``cmake`` project configuration file in the right location so that third-party projects can use cmake's ``find_package`` to locate xtl headers.

.. image:: conda.svg

Using the conda package
-----------------------

A package for xtl is available on the conda package manager.

.. code::

    conda install -c conda-forge xtl 

.. image:: spack.svg

Using the Spack package
-----------------------

A package for xtl is available on the Spack package manager.

.. code::

    spack install xtl
    spack load xtl

.. image:: cmake.svg

From source with cmake
----------------------

You can also install ``xtl`` from source with cmake. On Unix platforms, from the source directory:

.. code::

    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/path/to/prefix ..
    make install

On Windows platforms, from the source directory:

.. code::

    mkdir build
    cd build
    cmake -G "NMake Makefiles" -DCMAKE_INSTALL_PREFIX=/path/to/prefix ..
    nmake
    nmake install

See the section of the documentation on :doc:`build-options`, for more details on how to cmake options.
