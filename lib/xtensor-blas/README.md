# ![xtensor](docs/source/xtensor-blas.svg)

[![Travis](https://travis-ci.org/xtensor-stack/xtensor-blas.svg?branch=master)](https://travis-ci.org/xtensor-stack/xtensor-blas)
[![Appveyor](https://ci.appveyor.com/api/projects/status/i5c8u3q0uksx0m06?svg=true)](https://ci.appveyor.com/project/xtensor-stack/xtensor-blas)
[![Azure](https://dev.azure.com/xtensor-stack/xtensor-stack/_apis/build/status/xtensor-stack.xtensor-blas?branchName=master)](https://dev.azure.com/xtensor-stack/xtensor-stack/_build/latest?definitionId=5&branchName=master)
[![Documentation](http://readthedocs.org/projects/xtensor-blas/badge/?version=latest)](https://xtensor-blas.readthedocs.io/en/latest/?badge=latest)
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/xtensor-stack/xtensor/stable?filepath=notebooks%2Fxtensor.ipynb)
[![Join the Gitter Chat](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/QuantStack/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

## Introduction

`xtensor-blas` is an extension to the xtensor library, offering bindings to BLAS and LAPACK libraries through cxxblas and cxxlapack from the [FLENS](https://github.com/michael-lehn/FLENS) project.

`xtensor-blas` currently provides non-broadcasting `dot`, `norm` (1- and 2-norm for vectors), `inverse`, `solve`,
`eig`, `cross`, `det`, `slogdet`, `matrix_rank`, `inv`, `cholesky`, `qr`, `svd` in the `xt::linalg` namespace (check the corresponding `xlinalg.hpp` header for the function signatures). The functions, and signatures, are trying to be 1-to-1 equivalent to NumPy.
Low-level functions to interface with BLAS or LAPACK with xtensor containers are also offered in the `blas` and `lapack` namespace.

`xtensor` and `xtensor-blas` require a modern C++ compiler supporting C++14. The following C++ compilers are supported:

 - On Windows platforms, Visual C++ 2015 Update 2, or more recent
 - On Unix platforms, gcc 4.9 or a recent version of Clang

## Installation

xtensor-blas is a header-only library. We provide a package for the conda package manager.

```
conda install -c conda-forge xtensor-blas
```

which will also install the core `xtensor` package.

Or you can directly install it from the sources:

```
cmake -D CMAKE_INSTALL_PREFIX=your_install_prefix
make install
```

To build the tests or actually use `xtensor-blas`, you will need binaries for

 - `openblas`
 - `lapack`

which are also available on conda-forge.

## Trying it online

You can play with `xtensor` interactively in a Jupyter notebook right now! Just click on the binder link below:

[![Binder](binder-logo.svg)](https://mybinder.org/v2/gh/xtensor-stack/xtensor/stable?filepath=notebooks%2Fxtensor.ipynb)

The C++ support in Jupyter is powered by the [xeus-cling](https://github.com/xtensor-stack/xeus-cling) C++ kernel. Together with xeus-cling, xtensor enables a similar workflow to that of NumPy with the IPython Jupyter kernel.

## Documentation

For more information on using `xtensor`, check out the reference documentation

http://xtensor-blas.readthedocs.io/

## Dependency on `xtensor`

`xtensor-blas` depends on the `xtensor` package

| `xtensor-blas`  | `xtensor` |
|-----------------|-----------|
| master          |  ^0.21.4  |
| 0.17.2          |  ^0.21.4  |
| 0.17.1          |  ^0.21.2  |
| 0.17.0          |  ^0.21.1  |
| 0.16.1          |  ^0.20.4  |
| 0.16.0          |  ^0.20.0  |
| 0.15.2          |  ^0.19.0  |
| 0.15.1          |  ^0.19.0  |
| 0.15.0          |  ^0.19.0  |
| 0.14.0          |  ^0.18.0  |
| 0.13.x          |  ^0.17.0  |
| 0.12.x          |  ^0.17.0  |

## License

We use a shared copyright model that enables all contributors to maintain the
copyright on their contributions.

This software is licensed under the BSD-3-Clause license. See the [LICENSE](LICENSE) file for details.
