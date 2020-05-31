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
#include "xtensor/xcomplex.hpp"
#include "xtensor/xrandom.hpp"
#include "xtensor/xio.hpp"
#include "xtensor-blas/xblas.hpp"
#include "xtensor-blas/xlapack.hpp"
#include "xtensor-blas/xlinalg.hpp"

using namespace std::complex_literals;

namespace xt
{
    TEST(xlinalg, matrixpower)
    {
        xarray<double> t1arg_0 = {{0, 1, 2},
                                  {3, 4, 5},
                                  {6, 7, 8}};

        auto t1res = xt::linalg::matrix_power(t1arg_0, 2);
        xarray<double> t1expected = {{ 15, 18, 21},
                                     { 42, 54, 66},
                                     { 69, 90,111}};
        EXPECT_TRUE(allclose(t1res, t1expected));

        auto t2res = xt::linalg::matrix_power(t1arg_0, 5);
        xarray<double> t2expected = {{ 32400,  41796,  51192},
                                     { 99468, 128304, 157140},
                                     {166536, 214812, 263088}};
        EXPECT_TRUE(allclose(t2res, t2expected));

        auto t3res = xt::linalg::matrix_power(t1arg_0, 41);
        xarray<double> t3expected = {{  1.06199622e+45,  1.36986674e+45,  1.67773727e+45},
                                     {  3.26000325e+45,  4.20507151e+45,  5.15013977e+45},
                                     {  5.45801029e+45,  7.04027628e+45,  8.62254226e+45}};
        EXPECT_TRUE(allclose(t3res, t3expected));

        xarray<double> t4arg_0 = {{-2., 1., 3.},
                                  { 3., 2., 1.},
                                  { 1., 2., 5.}};

        auto t4res = xt::linalg::matrix_power(t4arg_0, -2);
        xarray<double> t4expected = {{ 0.09259259,-0.09259259, 0.01851852},
                                     { 0.35185185, 0.64814815,-0.46296296},
                                     {-0.2037037 ,-0.2962963 , 0.25925926}};
        EXPECT_TRUE(allclose(t4res, t4expected));

        auto t5res = xt::linalg::matrix_power(t4arg_0, -13);
        xarray<double> t5expected = {{-0.02119919,-0.02993041, 0.02400524},
                                     { 0.15202629, 0.21469317,-0.17217602},
                                     {-0.0726041 ,-0.10253451, 0.08222825}};
        EXPECT_TRUE(allclose(t5res, t5expected));
    }

    TEST(xlinalg, det)
    {
        xarray<double> a = {{1,2}, {3,4}};
        double da = linalg::det(a);
        EXPECT_EQ(da, -2.0);
        xarray<double> b = {{0, 1,2}, {3,4, 5}, {6,7,8}};
        double db = linalg::det(b);
        EXPECT_EQ(db, 0.0);
        xarray<double> c = {{12, 1,2}, {3,4, 5}, {6,7,8}};
        double dc = linalg::det(c);
        EXPECT_NEAR(dc, -36, 1e-06);

        xarray<std::complex<double>> arg_0 = {{ 0.95368636+0.32324664i, 0.49936872+0.22164004i, 0.30452434+0.78922905i},
                                              { 0.84118920+0.59652768i, 0.42052057+0.97211559i, 0.19916742+0.83068058i},
                                              { 0.67065616+0.56830636i, 0.00268706+0.29410473i, 0.69147455+0.7052149i }};
        auto res = linalg::det(arg_0);

        auto res_i = std::imag(res);
        auto res_r = std::real(res);
        EXPECT_NEAR(0.4201495908415372, res_i, 1e-06);
        EXPECT_NEAR(-0.07633013993862534, res_r, 1e-06);
    }

    TEST(xlinalg, slogdet)
    {
        xarray<std::complex<double>> arg_0 = {{ 0.13373658+0.43025551i, 0.42593478+0.17539337i, 0.18840853+0.24669458i},
                                              { 0.82800224+0.11797823i, 0.40310379+0.14037109i, 0.88204561+0.96870283i},
                                              { 0.35427657+0.1233739i , 0.22740960+0.94019582i, 0.05410180+0.86462543i}};
        auto resc = linalg::slogdet(arg_0);
        auto sc = std::get<0>(resc);
        auto sl = std::real(std::get<1>(resc));
        auto scr = std::real(sc);
        auto sci = std::imag(sc);

        EXPECT_NEAR(-0.8818794751329891, sl, 1e-06);
        EXPECT_NEAR(0.8473375077176295, scr, 1e-06);
        EXPECT_NEAR(0.5310547504870624, sci, 1e-06);

        xarray<double> arg_b = {{ 0.20009016, 0.33997118, 0.74433611},
                                { 0.52721448, 0.2449798 , 0.49085606},
                                { 0.49757477, 0.97304175, 0.05011255}};
        auto res = linalg::slogdet(arg_b);
        double expected_0 = 1.0;
        double expected_1 = -1.3017524147193602;
        auto sres = std::get<0>(res);
        auto lres = std::get<1>(res);
        EXPECT_EQ(expected_0, sres);
        EXPECT_NEAR(expected_1, lres, 1e-06);
    }

    TEST(xlinalg, svd)
    {
        xarray<double> arg_0 = {{0,1,2},
                                {3,4,5},
                                {6,7,8}};

        auto res = linalg::svd(arg_0);

        xarray<double, layout_type::column_major> expected_0 = {{-0.13511895, 0.90281571, 0.40824829},
                                                                {-0.49633514, 0.29493179,-0.81649658},
                                                                {-0.85755134,-0.31295213, 0.40824829}};
        xarray<double, layout_type::column_major> expected_1 = {  1.42267074e+01,  1.26522599e+00,  5.89938022e-16};
        xarray<double, layout_type::column_major> expected_2 = {{-0.4663281 ,-0.57099079,-0.67565348},
                                                                {-0.78477477,-0.08545673, 0.61386131},
                                                                {-0.40824829, 0.81649658,-0.40824829}};

        EXPECT_TRUE(allclose(std::get<0>(res), expected_0));
        EXPECT_TRUE(allclose(std::get<1>(res), expected_1));
        EXPECT_TRUE(allclose(std::get<2>(res), expected_2));
    }

    TEST(xlinalg, svd_horizontal_vertical)
    {
        xarray<double> a = xt::ones<double>({3, 1});
        xarray<double> b = xt::ones<double>({1, 3});
        xarray<double> u, s, vt;

        std::tie(u, s, vt) = linalg::svd(a, false);
        EXPECT_TRUE(allclose(a, xt::linalg::dot(u * s, vt)));

        std::tie(u, s, vt) = linalg::svd(b, false);
        EXPECT_TRUE(allclose(b, xt::linalg::dot(u * s, vt)));
    }

    TEST(xlinalg, matrix_rank)
    {
        xarray<double> eall = eye<double>(4);
        int a = linalg::matrix_rank(eall);
        EXPECT_EQ(4, a);

        xarray<double> b = eye<double>(4);
        b(1, 1) = 0;
        int rb = linalg::matrix_rank(b);
        EXPECT_EQ(3, rb);
        xarray<double> ones_arr = ones<double>({4, 4});
        int ro = linalg::matrix_rank(ones_arr);
        EXPECT_EQ(1, ro);
        xarray<double> zarr = zeros<double>({4, 4});
        int rz = linalg::matrix_rank(zarr);
        EXPECT_EQ(0, rz);
    }

    TEST(xlinalg, eigh)
    {
        xarray<double> arg_0 = {{ -761. , -208. , -582. },
                                { -208. , -623. ,-1605.5},
                                { -582. ,-1605.5, -476. }};
        auto res = xt::linalg::eigh(arg_0);
        xarray<double> expected_0 = {-2351.3290686 , -609.79206435, 1101.12113295};
        xarray<double, layout_type::column_major> expected_1 = {{-0.33220683,-0.93041946,-0.15478453},
                                                                {-0.66309119, 0.34708777,-0.66320446},
                                                                {-0.67078216, 0.11768479, 0.73225787}};
        auto vals = std::get<0>(res);
        auto vecs = std::get<1>(res);
        EXPECT_TRUE(allclose(expected_0, vals));
        EXPECT_TRUE(allclose(expected_1, vecs));

        auto vals_2 = xt::linalg::eigvalsh(arg_0);
        EXPECT_TRUE(allclose(expected_0, vals_2));

        xarray<std::complex<double>> complarg_0 = {{ 1.+0.i,-0.-2.i},
                                                   { 0.+2.i, 5.+0.i}};
        auto complres = xt::linalg::eigh(complarg_0);

        xarray<double> complexpected_0 = { 0.17157288, 5.82842712};
        auto cmvals = std::get<0>(complres);
        auto cmvecs = std::get<1>(complres);
        EXPECT_TRUE(allclose(complexpected_0, cmvals));
        xarray<std::complex<double>, layout_type::column_major> complexpected_1 = {{-0.92387953+0.i        ,-0.38268343+0.i        },
                                                                                   { 0.00000000+0.38268343i, 0.00000000-0.92387953i}}; 
        EXPECT_TRUE(allclose(imag(complexpected_1), imag(cmvecs)));
        EXPECT_TRUE(allclose(real(complexpected_1), real(cmvecs)));

        auto cmvals2 = xt::linalg::eigvalsh(complarg_0);
        EXPECT_TRUE(allclose(complexpected_0, cmvals2));
    }

    TEST(xlinalg, pinv)
    {
        xarray<double> arg_0 = {{ 1.47351391, 0.94686323, 0.92236842,-1.44141916,-1.53123963,-0.36949144},
                                {-0.76686921,-0.01087083,-1.11100036,-0.59745592,-0.99849726, 0.45296729},
                                {-0.35274989,-1.27760231, 1.50092545,-2.7243503 ,-0.79326768,-1.00826405},
                                { 0.05763039, 1.04069983,-0.502178  ,-1.01776144, 0.6496664 ,-0.2374513 },
                                {-1.45517735, 0.42523508, 0.41400096, 0.87164292, 1.87754145, 0.16358461},
                                { 1.07487297,-0.26417364, 1.82998799, 0.97985789,-0.74820612,-0.75097366},
                                { 0.91375249, 1.14211989,-0.23055478,-0.48264987,-0.4591723 , 0.83185472},
                                { 0.05318152,-0.30014836, 1.68456715, 0.07388112, 0.0607432 ,-0.51529535},
                                {-1.36227295,-0.12015569,-0.45599178,-1.07135129,-0.27405687, 0.50177945}};
        auto res = xt::linalg::pinv(arg_0);
        xarray<double> expected = {{-0.12524671,-0.41299325, 0.09108576,-0.07346514,-0.29603324,-0.1702256 ,
                      0.28503799,-0.08979346,-0.29415286},
                    { 0.38896886, 0.28355746,-0.26406943, 0.36828839, 0.29271933, 0.19368219,
                     -0.11433648, 0.05102866, 0.11050527},
                    { 0.05803293,-0.11353613, 0.08736754,-0.20744326, 0.21647828, 0.11201996,
                      0.26187311, 0.24550066, 0.13766844},
                    { 0.01037663, 0.15710257,-0.24647364,-0.06611398, 0.0482986 , 0.21835662,
                     -0.18763349, 0.01747553,-0.02066564},
                    {-0.23681835,-0.44982904, 0.14998634, 0.07608489, 0.04751429,-0.27436861,
                      0.18626331,-0.01520331,-0.18144752},
                    {-0.25580995,-0.33606213, 0.10424938,-0.64100544, 0.01084618,-0.28963573,
                      0.89629794, 0.12987378, 0.22451995}};
        EXPECT_TRUE(allclose(expected, res));

        xarray<std::complex<double>> cmpl_arg_0 = {{-0.32865615+1.56868725i, 0.28804396+0.52266479i},
                                                   {-1.29703842+0.34647524i,-2.14982936+0.31425111i},
                                                   {-0.69224750-1.36725801i, 2.22948403+1.4612309i }};
        auto cmpl_res = xt::linalg::pinv(cmpl_arg_0);
        xarray<std::complex<double>> cmpl_expected = {{-0.06272312-0.24840107i,-0.20530381-0.00548715i,-0.14179276+0.16337684i},
                         { 0.05975312-0.0502577i ,-0.17431091-0.05525696i, 0.16047967-0.14140846i}};
        EXPECT_TRUE(allclose(real(cmpl_expected), real(cmpl_res)));
        EXPECT_TRUE(allclose(imag(cmpl_expected), imag(cmpl_res)));
    }

    TEST(xlinalg, pinv_small)
    {
        xt::xtensor<float, 2> d1 {{1.f, 2.f}};
        auto r1 = xt::linalg::pinv(d1);
        xt::xtensor<float, 2> e1 = {{ 0.2f}, { 0.4f}};
        EXPECT_TRUE(allclose(r1, e1));

        xt::xtensor<float, 2> d2 {{ 1.f}};
        auto r2 = xt::linalg::pinv(d2);
        EXPECT_EQ(r2(0), 1.f);
    }

    TEST(xlinalg, mat_norm)
    {
        xarray<double> arg_0 = {{ 0.06817001, 0.50274712,-0.36802027,-0.93123204},
                                {-0.5990272 ,-0.67439921,-0.09397038,-1.55915724},
                                { 2.22694395, 0.59099048,-0.43162172, 0.19410077},
                                { 0.41859591, 1.68555153, 1.82660739, 1.24427635}};
        auto res1 = xt::linalg::norm(arg_0, 1);
        auto res2 = xt::linalg::norm(arg_0, linalg::normorder::frob);
        auto res3 = xt::linalg::norm(arg_0, linalg::normorder::inf);
        auto res4 = xt::linalg::norm(arg_0, linalg::normorder::neg_inf);
        auto res5 = xt::linalg::norm(arg_0, linalg::normorder::nuc);
        auto res6 = xt::linalg::norm(arg_0, 2);
        double exp1 = 3.92876639061;
        double exp2 = 4.23639347394;
        double exp3 = 5.17503118283;
        double exp4 = 1.87016943835;
        double exp5 = 7.42677006218;
        double exp6 = 3.29152325862;

        EXPECT_NEAR(exp1, res1, 1e-06);
        EXPECT_NEAR(exp2, res2, 1e-06);
        EXPECT_NEAR(exp3, res3, 1e-06);
        EXPECT_NEAR(exp4, res4, 1e-06);
        EXPECT_NEAR(exp5, res5, 1e-06);
        EXPECT_NEAR(exp6, res6, 1e-06);

        xarray<std::complex<double>> cmplarg_0 = {{ 0.40101756+0.71233018i, 0.62731701+0.42786349i, 0.32415089+0.2977805i },
                                                  { 0.24475928+0.49208478i, 0.69475518+0.74029639i, 0.59390240+0.35772892i},
                                                  { 0.63179202+0.41720995i, 0.44025718+0.65472131i, 0.08372648+0.37380143i}};
        auto cmplres1 = xt::linalg::norm(cmplarg_0, 1);
        auto cmplres2 = xt::linalg::norm(cmplarg_0, linalg::normorder::frob);
        auto cmplres3 = xt::linalg::norm(cmplarg_0, linalg::normorder::inf);
        auto cmplres4 = xt::linalg::norm(cmplarg_0, linalg::normorder::neg_inf);
        auto cmplres5 = xt::linalg::norm(cmplarg_0, linalg::normorder::nuc);
        auto cmplres6 = xt::linalg::norm(cmplarg_0, 2);

        double cmplexp1 = 2.56356133004;
        double cmplexp2 = 2.14347558031;
        double cmplexp3 = 2.25815855456;
        double cmplexp4 = 1.92915797164;
        double cmplexp5 = 2.77947580342;
        double cmplexp6 = 2.0683368289;

        EXPECT_NEAR(cmplexp1, cmplres1, 1e-06);
        EXPECT_NEAR(cmplexp2, cmplres2, 1e-06);
        EXPECT_NEAR(cmplexp3, cmplres3, 1e-06);
        EXPECT_NEAR(cmplexp4, cmplres4, 1e-06);
        EXPECT_NEAR(cmplexp5, cmplres5, 1e-06);
        EXPECT_NEAR(cmplexp6, cmplres6, 1e-06);
    }

    TEST(xlinalg, vec_norm)
    {
        xarray<double> arg_0 = { 0.23451288, 0.98799529, 0.76599595, 0.77700444, 0.02798196};

        EXPECT_NEAR(2.79349050582, xt::linalg::norm(arg_0, 1), 1e-6);
        EXPECT_NEAR(1.49077149771, xt::linalg::norm(arg_0, 2), 1e-6);
        EXPECT_NEAR(1.23766843269, xt::linalg::norm(arg_0, 3), 1e-6);
        EXPECT_NEAR(1.13587319901, xt::linalg::norm(arg_0, 4), 1e-6);
        EXPECT_NEAR(5.0, xt::linalg::norm(arg_0, 0), 1e-6);
        EXPECT_NEAR(0.0229325662443, xt::linalg::norm(arg_0, -1), 1e-6);
        EXPECT_NEAR(0.0277379546324, xt::linalg::norm(arg_0, -2), 1e-6);
        EXPECT_NEAR(0.987995286517, xt::linalg::norm(arg_0, linalg::normorder::inf), 1e-6);
        EXPECT_NEAR(0.0279819550429, xt::linalg::norm(arg_0, linalg::normorder::neg_inf), 1e-6);

        xarray<std::complex<double>> arg_1 = { 0.23451288+0.77700444i, 0.98799529+0.02798196i, 0.76599595+0.17390652i};
        EXPECT_NEAR(2.58550383197, xt::linalg::norm(arg_1, 1), 1e-06);
        EXPECT_NEAR(1.50088078633, xt::linalg::norm(arg_1, 2), 1e-06);
        EXPECT_NEAR(1.25673399279, xt::linalg::norm(arg_1, 3), 1e-06);
        EXPECT_NEAR(1.15326879931, xt::linalg::norm(arg_1, 4), 1e-06);
        EXPECT_NEAR(3.0, xt::linalg::norm(arg_1, 0), 1e-06);
        EXPECT_NEAR(0.284338433895, xt::linalg::norm(arg_1, -1), 1e-06);
        EXPECT_NEAR(0.490145522524, xt::linalg::norm(arg_1, -2), 1e-06);
        EXPECT_NEAR(0.98839145888, xt::linalg::norm(arg_1, linalg::normorder::inf), 1e-06);
        EXPECT_NEAR(0.785489192861, xt::linalg::norm(arg_1, linalg::normorder::neg_inf), 1e-06);
    }

    TEST(xlinalg, vdot)
    {
        xarray<double> arg_0 = { 0.23451288, 0.98799529, 0.76599595, 0.77700444, 0.02798196};
        xarray<double> arg_1 = { 0.17390652, 0.15408224, 0.07708648, 0.8898657 , 0.7503787 };
        auto res = xt::linalg::vdot(arg_0, arg_1);
        EXPECT_NEAR(0.964490439715, res, 1e-06);

        xarray<std::complex<double>> carg_0 = { 0.23451288+0.17390652i, 0.98799529+0.15408224i, 0.76599595+0.07708648i,
                                                0.77700444+0.8898657i , 0.02798196+0.7503787i };
        xarray<std::complex<double>> carg_1 = { 0.17390652+0.23451288i, 0.15408224+0.98799529i, 0.07708648+0.76599595i,
                                                0.88986570+0.77700444i, 0.75037870+0.02798196i};
        auto res_c = xt::linalg::vdot(carg_0, carg_1);

        EXPECT_NEAR(1.9289808794290355, std::real(res_c), 1e-06);
        EXPECT_NEAR(0.8075433553117102, std::imag(res_c), 1e-06);
    }

    TEST(xlinalg, kron)
    {
        xarray<int> arg_0 = {{2,1,8},
                             {3,5,0},
                             {2,6,2},
                             {4,4,6}};

        xarray<int> arg_1 = {{3,0,6,4,7},
                             {6,7,1,5,7}};

        auto res = xt::linalg::kron(arg_0, arg_1);

        xarray<int> expected = {{ 6, 0,12, 8,14, 3, 0, 6, 4, 7,24, 0,48,32,56},
                                {12,14, 2,10,14, 6, 7, 1, 5, 7,48,56, 8,40,56},
                                { 9, 0,18,12,21,15, 0,30,20,35, 0, 0, 0, 0, 0},
                                {18,21, 3,15,21,30,35, 5,25,35, 0, 0, 0, 0, 0},
                                { 6, 0,12, 8,14,18, 0,36,24,42, 6, 0,12, 8,14},
                                {12,14, 2,10,14,36,42, 6,30,42,12,14, 2,10,14},
                                {12, 0,24,16,28,12, 0,24,16,28,18, 0,36,24,42},
                                {24,28, 4,20,28,24,28, 4,20,28,36,42, 6,30,42}};

        EXPECT_EQ(expected, res);
    }

    TEST(xlinalg, cholesky)
    {
        xarray<double> arg_0 = {{  4, 12,-16},
                                { 12, 37,-43},
                                {-16,-43, 98}};

        auto res = xt::linalg::cholesky(arg_0);
        xarray<double> expected = {{ 2., 0., 0.},
                                   { 6., 1., 0.},
                                   {-8., 5., 3.}};
        EXPECT_EQ(expected, res);

        xarray<std::complex<double>> cmplarg_0 = {{ 1.+0.i,-0.-2.i},
                                                  { 0.+2.i, 5.+0.i}};
        auto cmplres = xt::linalg::cholesky(cmplarg_0);
        xarray<std::complex<double>> cmplexpected = {{ 1.+0.i, 0.+0.i},
                                                     { 0.+2.i, 1.+0.i}};
        EXPECT_EQ(cmplexpected, cmplres);
    }

    TEST(xlinalg, qr)
    {
        xarray<double, layout_type::column_major> a = xt::random::rand<double>({9, 6});
        auto res = xt::linalg::qr(a);
        xarray<double> q = std::get<0>(res);
        xarray<double> r = std::get<1>(res);
        auto resf = xt::linalg::qr(a, linalg::qrmode::complete);
        auto resr = xt::linalg::qr(a, linalg::qrmode::r);
        xarray<double> qf = std::get<0>(resf);
        xarray<double> rf = std::get<1>(resf);

        auto neara = xt::linalg::dot(q, r);
        EXPECT_TRUE(allclose(neara, a));
        auto nearaf = xt::linalg::dot(qf, rf);
        EXPECT_TRUE(allclose(nearaf, a));

        EXPECT_EQ(std::get<1>(resr), xt::view(rf, xt::range(0, 6), xt::all()));
        EXPECT_EQ(std::get<0>(resr).size(), 0u);
        EXPECT_EQ(std::get<0>(resr).dimension(), 1u);

        xarray<double, layout_type::column_major> erawR = {{-1.00444014e+01,  0.00000000e+00,  6.74440143e-01, 2.24813381e-01},
                                                           {-9.58743044e+00, -1.25730337e+01, -6.22814365e-03, 3.37562246e-01},
                                                           {-1.29027101e+01, -7.34080303e+00, -4.07831856e+00, -5.76331089e-01}};

        xarray<double, layout_type::column_major> eTau = { 1.32854123, 1.79535299, 1.50132395};

        xarray<double, layout_type::column_major> AA = {{ 3.3,  1.,  2.},
                                                        { 0. , 10.,  8.},
                                                        { 9. ,  7., 12.},
                                                        { 3. , 10.,  5.}};

        auto resraw = xt::linalg::qr(AA, linalg::qrmode::raw);
        auto tau = std::get<1>(resraw);
        auto rawR = std::get<0>(resraw);

        EXPECT_TRUE(allclose(tau, eTau));
        EXPECT_TRUE(allclose(erawR, rawR));
    }

    TEST(xlinalg, lstsq)
    {
        xarray<double> arg_0 = {{ 0., 1.},
                                { 1., 1.},
                                { 2., 1.},
                                { 3., 1.}};

        xarray<double> arg_1 = {{-1., 0.2, 0.9, 2.1}, 
                                { 2., 3. , 2. , 1. }};
        arg_1 = transpose(arg_1);
        auto res = xt::linalg::lstsq(arg_0, arg_1);

        xarray<double, layout_type::column_major> el_0 = {{ 1.  ,-0.4 },
                                                          {-0.95, 2.6 }};
        xarray<double> el_1 = { 0.05, 1.2 };
        int el_2 = 2;
        xarray<double> el_3 = { 4.10003045, 1.09075677};

    
        EXPECT_TRUE(allclose(el_0, std::get<0>(res)));
        EXPECT_TRUE(allclose(el_1, std::get<1>(res)));
        EXPECT_EQ(el_2, std::get<2>(res));
        EXPECT_TRUE(allclose(el_3, std::get<3>(res)));

        xarray<std::complex<double>> carg_0 = {{ 0., 1.},
                                               { 1. - 3i, 1.},
                                               { 2., 1.},
                                               { 3., 1.}};
        xarray<std::complex<double>> carg_1 = {{-1. , 0.2+4i, 0.9, 2.1-1i}, {2,3i,2,1}};
        carg_1 = transpose(carg_1);
        auto cres = xt::linalg::lstsq(carg_0, carg_1);

        xarray<std::complex<double>, layout_type::column_major> cel_0 = {{-0.40425532-0.38723404i,-0.61702128-0.44680851i},
                                                                         { 1.44680851+1.02765957i, 2.51063830+0.95744681i}};
        xarray<double> cel_1 = { 16.11787234,  2.68085106};
        int cel_2 = 2;
        xarray<double> cel_3 = { 5.01295356, 1.36758789};

        EXPECT_TRUE(allclose(imag(cel_0), imag(std::get<0>(cres))));
        EXPECT_TRUE(allclose(real(cel_0), real(std::get<0>(cres))));
        EXPECT_TRUE(allclose(cel_1, std::get<1>(cres)));
        EXPECT_EQ(cel_2, std::get<2>(cres));
        EXPECT_TRUE(allclose(cel_3, std::get<3>(cres)));
    }

    TEST(xlinalg, trace)
    {
        auto e1 = eye<double>(10);
        xarray<double> e2 = eye<double>(5);

        auto t1 = linalg::trace(e1);
        auto t11 = linalg::trace(e1, 1);
        auto t1n1 = linalg::trace(e1, -1);
        EXPECT_EQ(10, t1());
        EXPECT_EQ(0, t11());
        EXPECT_EQ(0, t1n1());

        auto t2 = linalg::trace(e2);
        auto t22 = linalg::trace(e2, 1);
        EXPECT_EQ(5, t2());
        EXPECT_EQ(0, t22());

        xarray<double> ar = xt::arange(9);
        ar.reshape({3, 3});

        auto ar1 = linalg::trace(ar);
        auto ar2 = linalg::trace(ar, 1);
        auto ar3 = linalg::trace(ar, -1);

        EXPECT_EQ(12, ar1());
        EXPECT_EQ(6,  ar2());
        EXPECT_EQ(10, ar3());
    }

    TEST(xlinalg, dots)
    {
        xarray<double> arg_0 = {{{ 0, 1, 2},
                                 { 3, 4, 5},
                                 { 6, 7, 8}},

                                {{ 9,10,11},
                                 {12,13,14},
                                 {15,16,17}}};

        xarray<double> arg_1 = {{{ 0, 1, 2},
                                 { 3, 4, 5},
                                 { 6, 7, 8}},

                                 {{ 9,10,11},
                                  {12,13,14},
                                  {15,16,17}},

                                 {{18,19,20},
                                  {21,22,23},
                                  {24,25,26}}};

        auto res1 = xt::linalg::dot(arg_0, arg_1);
        xarray<double> expected1 = {{{{  15,  18,  21},
                                      {  42,  45,  48},
                                      {  69,  72,  75}},

                                     {{  42,  54,  66},
                                      { 150, 162, 174},
                                      { 258, 270, 282}},

                                     {{  69,  90, 111},
                                      { 258, 279, 300},
                                      { 447, 468, 489}}},
                                
                                    {{{  96, 126, 156},
                                      { 366, 396, 426},
                                      { 636, 666, 696}},

                                     {{ 123, 162, 201},
                                      { 474, 513, 552},
                                      { 825, 864, 903}},
                                
                                     {{ 150, 198, 246},
                                      { 582, 630, 678},
                                      {1014,1062,1110}}}};

        EXPECT_TRUE(allclose(expected1, res1));

        auto res2 = xt::linalg::dot(arg_1, arg_0);
        xarray<double> expected2 = {{{{  15,  18,  21},
                                      {  42,  45,  48}},

                                     {{  42,  54,  66},
                                      { 150, 162, 174}},

                                     {{  69,  90, 111},
                                      { 258, 279, 300}}},

                                    {{{  96, 126, 156},
                                      { 366, 396, 426}},

                                     {{ 123, 162, 201},
                                      { 474, 513, 552}},

                                     {{ 150, 198, 246},
                                      { 582, 630, 678}}},

                                    {{{ 177, 234, 291},
                                      { 690, 747, 804}},

                                     {{ 204, 270, 336},
                                      { 798, 864, 930}},

                                     {{ 231, 306, 381},
                                      { 906, 981,1056}}}};

        EXPECT_TRUE(allclose(expected2, res2));

        xarray<double> arg_2 = {0, 1, 2};
        auto res3 = xt::linalg::dot(arg_0, arg_2);

        xarray<double> expected3 = {{ 5, 14, 23},
                                    {32, 41, 50}};

        EXPECT_TRUE(allclose(expected3, res3));

        auto res4 = xt::linalg::dot(arg_2, arg_0);

        xarray<double> expected4 = {{15, 18, 21},
                                    {42, 45, 48}};

        EXPECT_TRUE(allclose(expected4, res4));
    }

    TEST(xlinalg, negative_strides)
    {
        xt::xarray<double> A = {{2, 3}, {5, 7}, {11, 13}};

        auto A1 = xt::view(A, xt::range(0, 3), 0);
        auto A2 = xt::view(A, xt::range(-1, -4, -1), 1);

        auto res = xt::linalg::dot(A1, A2);
        EXPECT_EQ(res(), 94);
    }

    TEST(xlinalg, asserts)
    {
        EXPECT_THROW(xt::linalg::eigh(xt::ones<double>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::eig(xt::ones<double>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::solve(xt::ones<double>({3, 1}), xt::ones<double>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::inv(xt::ones<double>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::eigvals(xt::ones<double>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::eigvalsh(xt::ones<double>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::det(xt::ones<double>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::slogdet(xt::ones<double>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::cholesky(xt::ones<double>({3, 1})), std::runtime_error);

        EXPECT_THROW(xt::linalg::eigh(xt::ones<std::complex<double>>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::eig(xt::ones<std::complex<double>>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::solve(xt::ones<std::complex<double>>({3, 1}), xt::ones<std::complex<double>>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::inv(xt::ones<std::complex<double>>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::eigvals(xt::ones<std::complex<double>>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::eigvalsh(xt::ones<std::complex<double>>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::det(xt::ones<std::complex<double>>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::slogdet(xt::ones<std::complex<double>>({3, 1})), std::runtime_error);
        EXPECT_THROW(xt::linalg::cholesky(xt::ones<std::complex<double>>({3, 1})), std::runtime_error);
    }

}
