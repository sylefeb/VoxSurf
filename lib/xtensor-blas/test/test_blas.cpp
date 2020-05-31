/***************************************************************************
* Copyright (c) Wolf Vollprecht, Johan Mabille and Sylvain Corlay          *
* Copyright (c) QuantStack                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/

#include "gtest/gtest.h"
#include "xtensor/xarray.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xrandom.hpp"

#include "xtensor-blas/xblas.hpp"
#include "xtensor-blas/xlinalg.hpp"

namespace xt
{
    TEST(xblas, matrix_times_vector)
    {
        xarray<double> m1 = {{1, 2, 3},
                             {4, 5, 6}};
        xarray<double> b = {1, 2, 3};

        auto res = linalg::dot(m1, b);
        xarray<double> expected = {14, 32};
        EXPECT_EQ(expected, res);

        xarray<double> next_row = {{7, 8, 9}};
        auto res2 = linalg::dot(concatenate(xtuple(m1, next_row)), b);
        xarray<double> expected2 = {14, 32, 50};
        EXPECT_EQ(expected2, res2);
    }

    TEST(xblas, dot_2d)
    {
        xarray<double> a = {{1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}};
        xarray<double> b = {5, 4, 3, 2, 1};
        xarray<double> c = {1, 2};
        auto res = linalg::dot(a, b);
        xarray<double> expected = {35, 35};

        auto res_ca = linalg::dot(c, a);
        xarray<double> expected_ca = {3, 6, 9, 12, 15};

        EXPECT_EQ(expected, res);
        EXPECT_EQ(expected_ca, res_ca);
    }

    TEST(xblas, matrix_matrix)
    {
        xarray<double> a = arange(3 * 3);
        a.reshape({3, 3});
        xarray<double> b = arange(5 * 3);
        b.reshape({3, 5});

        auto ab = linalg::dot(a, b);
        xarray<double> ab_expected ={{ 25,  28,  31,  34,  37},
                                     { 70,  82,  94, 106, 118},
                                     {115, 136, 157, 178, 199}};
        EXPECT_EQ(ab_expected, ab);
    }

    TEST(xblas, view_dot)
    {
        xarray<double> a = {1, 2, 3, 4, 5};
        xarray<double> b = {5, 4, 3, 2, 1};
        auto res = linalg::dot(a, b);

        xarray<double> expected = {35};
        EXPECT_EQ(expected, res);

        xarray<double> m1
          {{1, 2, 3},
           {4, 5, 6}};

        xarray<double> c = {1, 2};
        auto res2 = xt::linalg::dot(xt::view(m1, xt::all(), 1), c);
        xarray<double> expected2 = {12};
        EXPECT_EQ(expected2, res2);
    }

    TEST(xblas, norm)
    {
        auto a = linalg::norm(xt::arange<double>(15), 1);
        auto b = linalg::norm(xt::arange<double>(15), 2);
        xarray<double> c = {6, 4, 2, 1};
        auto res = linalg::norm(c);

        EXPECT_EQ(a, 105.0);
        EXPECT_NEAR(b, 31.859064644147981, 1e-6);
        EXPECT_NEAR(res, 7.5498344352707498, 1e-6);
    }

    TEST(xblas, normFloat)
    {
        auto a = linalg::norm(xt::arange<float>(15), 1);
        auto b = linalg::norm(xt::arange<float>(15), 2);
        xarray<float> c = {6, 4, 2, 1};
        auto res = linalg::norm(c);

        EXPECT_EQ(a, 105.0);
        EXPECT_NEAR(b, 31.859064644147981, 1e-6);
        EXPECT_NEAR(res, 7.5498344352707498, 1e-6);
    }

    TEST(xblas, outer)
    {
        xarray<double> a = {1, 1, 1};

        xarray<double> b = arange(0, 3);

        xarray<double> expected = {{0, 1, 2},
                                   {0, 1, 2},
                                   {0, 1, 2}};

        auto t = linalg::outer(a, b);
        auto t2 = linalg::outer(a, xt::arange<double>(0, 3));
        auto t3 = linalg::outer(xt::ones<double>({3}), xt::arange<double>(0, 3));

        EXPECT_TRUE(all(equal(expected, t)));
        EXPECT_TRUE(all(equal(expected, t2)));
        EXPECT_TRUE(all(equal(expected, t3)));
    }

    TEST(xblas, outer_random)
    {
        xt::random::seed(123);
        xt::xarray<double> expected = xt::random::randn<double>({5});
        xt::random::seed(123);
        auto x = xt::random::randn<double>({5});
        auto weights = xt::xarray<double> ({1});  // should perform identity

        auto result = linalg::outer(x, weights);

        // shapes are different
        for (std::size_t i = 0; i < 5; ++i)
        {
            EXPECT_EQ(result.data()[i], expected.data()[i]);
        }
    }

    TEST(xblas, nan_result)
    {
        xt::xarray<double> X = {{1, 2, 3},
                                {1, 2, 3}};

        auto M = xt::xarray<double>::from_shape({3, 3});
        M(0, 0) = std::numeric_limits<double>::quiet_NaN();
        M(0, 1) = std::numeric_limits<double>::quiet_NaN();
        xt::blas::gemm(X, X, M, true, false, 1.0, 0.0);
        for (std::size_t i = 0; i < M.size(); ++i)
        {
            EXPECT_FALSE(std::isnan(M.storage()[i]));
        }
    }

    TEST(xblas, gemm_transpose)
    {
        xt::xarray<double> X = {{1, 2, 3},
                                {1, 2, 3}};

        auto M = xt::xarray<double>::from_shape({3, 3});
        auto O = xt::xarray<double>::from_shape({2, 2});

        xt::blas::gemm(X, X, M, true, false, 1.0, 0.0);
        xt::blas::gemm(X, X, O, false, true, 1.0, 0.0);

        xt::xarray<double> expM = {{ 2,  4,  6},
                                   { 4,  8, 12},
                                   { 6, 12, 18}};

        xt::xarray<double> expO = {{ 14, 14},
                                   { 14, 14}};

        EXPECT_TRUE(all(equal(expM, M)));
        EXPECT_TRUE(all(equal(expO, O)));
    }

    TEST(xblas, gemv_transpose)
    {
        xt::xarray<double> X = {{1, 2, 3},
                                {1, 2, 3}};
        xt::xarray<double> v = {1, 2};
        auto R = xt::xarray<double>::from_shape({3});

        xt::blas::gemv(X, v, R, true, 1, 0);

        xt::xarray<double> expR = {3, 6, 9};

        EXPECT_TRUE(all(equal(expR, R)));
    }
}
