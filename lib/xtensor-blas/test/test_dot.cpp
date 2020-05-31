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
#include "xtensor/xstrided_view.hpp"

#include "xtensor-blas/xlinalg.hpp"

namespace xt
{
    TEST(xdot, matrix_times_vector)
    {
        xarray<float> a = xt::ones<float>({1, 4});
        xarray<float> b = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {1, 1, 1}};

        xarray<float> e1 = {{13, 16, 19}};

        auto r1 = linalg::dot(a, b);
        EXPECT_EQ(e1, r1);

        xarray<float> c = xt::ones<float>({3, 1});

        auto r2 = linalg::dot(b, c);
        xarray<float> e2 = {{6, 15, 24, 3}};
        e2.reshape({4, 1});
        EXPECT_EQ(e2, r2);

        EXPECT_THROW(linalg::dot(b, a), std::runtime_error);
        EXPECT_THROW(linalg::dot(c, b), std::runtime_error);
    }

    TEST(xdot, matrix_transpose_times_column)
    {
        xarray<double, layout_type::row_major> a = xt::ones<double>({2, 4});
        xarray<double, layout_type::row_major> b = xt::ones<double>({2, 1});
        auto r1 = linalg::dot(xt::transpose(a), b);
        EXPECT_TRUE(all(equal(r1, 2.0)));
    }

    TEST(xdot, matrix_transpose_times_column_cm)
    {
        xarray<double, layout_type::column_major> a = xt::ones<double>({2, 4});
        xarray<double, layout_type::column_major> b = xt::ones<double>({2, 1});

        auto r1 = linalg::dot(xt::transpose(a), b);
        EXPECT_TRUE(all(equal(r1, 2.0)));
    }

    TEST(xdot, square_matrix_times_vector)
    {
        xarray<float> a = {{1, 1, 1}};
        xarray<float> b = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

        auto r1 = linalg::dot(a, b);

        xarray<float> e1 = {{12, 15, 18}};
        EXPECT_EQ(r1, e1);

        auto r2 = linalg::dot(b, xt::transpose(a));
        xarray<float> e2 = xarray<float>::from_shape({3, 1});
        e2(0, 0) = 6.f;
        e2(1, 0) = 15.f;
        e2(2, 0) = 24.f;
        EXPECT_EQ(r2, e2);

        EXPECT_THROW(linalg::dot(b, a), std::runtime_error);
    }

    TEST(xdot, vector_times_vector)
    {
        xarray<float> a = xt::ones<float>({1, 3});
        xarray<float> b = xt::ones<float>({3, 1});

        auto r1 = linalg::dot(a, b);

        xarray<float> e1 = xarray<float>::from_shape({1, 1});
        e1(0, 0) = 3;

        EXPECT_EQ(e1, r1);

        auto r2 = linalg::dot(b, a);
        xarray<float> e2 = xt::ones<float>({3, 3});
        EXPECT_EQ(e2, r2);

        auto r3 = linalg::dot(b, e1);
        EXPECT_EQ(b * 3.f, r3);
    }

    TEST(xdot, matrix_times_1d)
    {
        xarray<float> a = xt::ones<float>({5, 3});
        xarray<float> b = xt::ones<float>({5});
        xarray<float> c = xt::ones<float>({3});
        auto r1 = linalg::dot(xt::transpose(a), b);

        EXPECT_TRUE(all(equal(r1, 5.0)));

        auto r2 = linalg::dot(c, xt::transpose(a));
        EXPECT_TRUE(all(equal(r2, 3.0)));

        auto r3 = linalg::dot(a, c);
        EXPECT_TRUE(all(equal(r3, 3.0)));

        auto r4 = linalg::dot(c, xt::ones<float>({3, 5}));
        EXPECT_TRUE(all(equal(r3, 3.0)));
    }

    TEST(xdot, A_times_A_T)
    {
        xarray<float> a = xt::ones<float>({5, 3});

        auto r1 = linalg::dot(a, xt::transpose(a));
        EXPECT_TRUE(all(equal(r1, 3.0)));

        auto r2 = linalg::dot(xt::transpose(a), a);
        EXPECT_TRUE(all(equal(r2, 5.0)));
    }

    TEST(xdot, matrix_times_vector_cm)
    {
        xarray<float, layout_type::column_major> a = xt::ones<float>({1, 4});
        xarray<float, layout_type::column_major> b = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {1, 1, 1}};

        xarray<float, layout_type::column_major> e1 = {{13, 16, 19}};

        auto r1 = linalg::dot(a, b);
        EXPECT_EQ(e1, r1);

        xarray<float, layout_type::column_major> c = xt::ones<float>({3, 1});

        auto r2 = linalg::dot(b, c);
        xarray<float, layout_type::column_major> e2 = {{6, 15, 24, 3}};
        e2.reshape({4, 1});
        EXPECT_EQ(e2, r2);

        EXPECT_THROW(linalg::dot(b, a), std::runtime_error);
        EXPECT_THROW(linalg::dot(c, b), std::runtime_error);
    }

    TEST(xdot, square_matrix_times_vector_cm)
    {
        xarray<float, layout_type::column_major> a = {{1, 1, 1}};
        xarray<float, layout_type::column_major> b = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

        auto r1 = linalg::dot(a, b);

        xarray<float, layout_type::column_major> e1 = {{12, 15, 18}};
        EXPECT_EQ(r1, e1);

        auto r2 = linalg::dot(b, xt::transpose(a));
        xarray<float, layout_type::column_major> e2 = xarray<float, layout_type::column_major>::from_shape({3, 1});
        e2(0, 0) = 6.f;
        e2(1, 0) = 15.f;
        e2(2, 0) = 24.f;
        EXPECT_EQ(r2, e2);

        EXPECT_THROW(linalg::dot(b, a), std::runtime_error);
    }

    TEST(xdot, vector_times_vector_cm)
    {
        xarray<float, layout_type::column_major> a = xt::ones<float>({1, 3});
        xarray<float, layout_type::column_major> b = xt::ones<float>({3, 1});

        auto r1 = linalg::dot(a, b);

        xarray<float, layout_type::column_major> e1 = xarray<float, layout_type::column_major>::from_shape({1, 1});
        e1(0, 0) = 3;

        EXPECT_EQ(e1, r1);

        auto r2 = linalg::dot(b, a);
        xarray<float, layout_type::column_major> e2 = xt::ones<float>({3, 3});
        EXPECT_EQ(e2, r2);

        auto r3 = linalg::dot(b, e1);
        EXPECT_EQ(b * 3.f, r3);
    }

    TEST(xdot, on_view)
    {
        xt::xarray<int> a = xt::reshape_view(xt::arange<int>(10 * 10 * 10), {10, 10, 10});
        xt::xarray<int> b = xt::reshape_view(xt::arange<int>(10 * 10 * 10), {10, 10, 10});

        auto res = xt::linalg::dot(view(a, 0, 0), view(b, 0, 0));
        auto res1 = xt::linalg::dot(view(a, 0, range(0, 3)), transpose(view(b, 0, range(0, 3))));

        EXPECT_EQ(res(0), 285.);
        EXPECT_EQ(res1(0, 0), 285.);
        EXPECT_EQ(res1(1, 2), 3635.);

        EXPECT_EQ(res1.dimension(), 2u);
        EXPECT_EQ(res1.shape()[0], 3u);
        EXPECT_EQ(res1.shape()[1], 3u);

        auto res2 = xt::linalg::dot(strided_view(a, {0, 0}), strided_view(b, {0, 0}));
        auto res3 = xt::linalg::dot(strided_view(a, {0, range(0, 3)}), transpose(strided_view(b, {0, range(0, 3)})));
        EXPECT_EQ(res2(0), 285.);
        EXPECT_EQ(res3(0, 0), 285.);
        EXPECT_EQ(res3(1, 2), 3635.);

        EXPECT_EQ(res3.dimension(), 2u);
        EXPECT_EQ(res3.shape()[0], 3u);
        EXPECT_EQ(res3.shape()[1], 3u);
    }

}
