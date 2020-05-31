/***************************************************************************
* Copyright (c) 2016, Johan Mabille, Sylvain Corlay and Wolf Vollprecht    *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/

#ifndef BENCHMARK_BLAS_HPP
#define BENCHMARK_BLAS_HPP

#include <benchmark/benchmark.h>

#include "xtensor/xnoalias.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xarray.hpp"

#include "xtensor-blas/xlinalg.hpp"

namespace xt
{
    namespace benchmark_dot
    {

        /****************************
         * Benchmark initialization *
         ****************************/

        template <class V>
        inline void init_benchmark_data(V& lhs, V& rhs, std::size_t size0, std::size_t size1)
        {
            using T = typename V::value_type;
            for (std::size_t i = 0; i < size0; ++i)
            {
                for (std::size_t j = 0; j < size1; ++j)
                {
                    lhs(i, j) = T(0.5) * T(j) / T(j + 1) + std::sqrt(T(i)) * T(9.) / T(size1);
                    rhs(i, j) = T(10.2) / T(i + 2) + T(0.25) * T(j);
                }
            }
        }

        template <class V>
        inline void init_xtensor_benchmark(V& lhs, V& rhs,
                                           std::size_t size0, size_t size1)
        {
            lhs.reshape({ size0, size1 });
            rhs.reshape({ size0, size1 });
            init_benchmark_data(lhs, rhs, size0, size1);
        }

        template <class E>
        inline auto benchmark_dot(benchmark::State& state)
        {
            using size_type = typename E::size_type;

            E x, y;
            init_xtensor_benchmark(x, y, state.range(0), state.range(0));

            while (state.KeepRunning())
            {
                auto res = xt::linalg::dot(x, y);
                benchmark::DoNotOptimize(res.data());
            }
        }

        template <class E>
        inline auto benchmark_transpose_dot(benchmark::State& state)
        {
            using size_type = typename E::size_type;

            E x, y;
            init_xtensor_benchmark(x, y, state.range(0), state.range(0));

            while (state.KeepRunning())
            {
                auto res = xt::linalg::dot(xt::transpose(x), y);
                benchmark::DoNotOptimize(res.data());
            }
        }

        template <class E>
        inline auto benchmark_transpose_with_assign_dot(benchmark::State& state)
        {
            using size_type = typename E::size_type;

            E x, y;
            init_xtensor_benchmark(x, y, state.range(0), state.range(0));

            while (state.KeepRunning())
            {
                xtensor<double, 2> t = xt::transpose(x);
                auto res = xt::linalg::dot(t, y);
                benchmark::DoNotOptimize(res.data());
            }
        }

        BENCHMARK_TEMPLATE(benchmark_dot, xt::xtensor<double, 2>)->Range(32, 32<<3);
        BENCHMARK_TEMPLATE(benchmark_transpose_dot, xt::xtensor<double, 2>)->Range(32, 32<<3);
        BENCHMARK_TEMPLATE(benchmark_transpose_with_assign_dot, xt::xtensor<double, 2>)->Range(32, 32<<3);
    }
}

#endif
