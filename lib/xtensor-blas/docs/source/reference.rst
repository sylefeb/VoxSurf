.. Copyright (c) 2017, Wolf Vollprecht, Johan Mabille and Sylvain Corlay

   Distributed under the terms of the BSD 3-Clause License.

   The full license is in the file LICENSE, distributed with this software.


Reference
=========

Defined in ``xtensor-blas/xlinalg.hpp``

The functions here are closely modeled after NumPy's linalg package.

Matrix, vector and tensor products
--------------------------

.. doxygenfunction:: xt::linalg::dot
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::vdot
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::outer
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::matrix_power
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::kron
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::tensordot(const xexpression<T>&, const xexpression<O>&, std::size_t)
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::tensordot(const xexpression<T>&, const xexpression<O>&, const std::vector<std::size_t>&, const std::vector<std::size_t>&)
    :project: xtensor-blas

Decompositions
--------------

.. doxygenfunction:: xt::linalg::cholesky
    :project: xtensor-blas

.. doxygenenum:: xt::linalg::qrmode
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::qr
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::svd
    :project: xtensor-blas

Matrix eigenvalues
------------------

.. doxygenfunction:: xt::linalg::eig
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::eigvals
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::eigh
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::eigvalsh
    :project: xtensor-blas


Norms and other numbers
-----------------------

.. doxygenenum:: xt::linalg::normorder
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::norm(const xexpression<E>&, int)
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::norm(const xexpression<E>&, normorder)
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::norm(const xexpression<E>&)
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::det(const xexpression<E>&)
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::slogdet(const xexpression<E>&)
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::matrix_rank
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::trace
    :project: xtensor-blas

Solving equations and inverting matrices
----------------------------------------

.. doxygenfunction:: xt::linalg::solve
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::lstsq
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::inv
    :project: xtensor-blas

.. doxygenfunction:: xt::linalg::pinv
    :project: xtensor-blas

Other
-----

.. doxygenfunction:: xt::linalg::cross
    :project: xtensor-blas
