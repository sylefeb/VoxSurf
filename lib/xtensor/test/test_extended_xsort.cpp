/***************************************************************************
* Copyright (c) Johan Mabille, Sylvain Corlay and Wolf Vollprecht          *
* Copyright (c) QuantStack                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/

// This file is generated from test/files/cppy_source/test_extended_xsort.cppy by preprocess.py!

#include "gtest/gtest.h"

#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xsort.hpp>

namespace xt
{
    using namespace xt::placeholders;

    template <class T>
    bool check_partition_equal(const T& a1, const T& a2, std::size_t kth)
    {
        auto p = a1[kth];
        EXPECT_EQ(p, a2[kth]);

        for (std::size_t i = 0; i < kth; ++i)
        {
            EXPECT_TRUE(a1[i] < p);
            EXPECT_TRUE(a2[i] < p);
        }

        for (std::size_t i = kth + 1; i < a1.size(); ++i)
        {
            EXPECT_TRUE(p < a1[i]);
            EXPECT_TRUE(p < a2[i]);
        }
        return true;
    }

    template <class X, class Y, class Z>
    bool check_argpartition_equal(const X& data, const Y& a1, const Z& a2, std::size_t kth)
    {
        auto p = static_cast<std::size_t>(a1[kth]);
        EXPECT_EQ(p, std::size_t(a2[kth]));
        auto el = data[static_cast<std::size_t>(a1[kth])];
        for (std::size_t i = 0; i < kth; ++i)
        {
            EXPECT_TRUE(data[static_cast<std::size_t>(a1[i])] < el);
            EXPECT_TRUE(data[static_cast<std::size_t>(a2[i])] < el);
        }

        for (std::size_t i = kth + std::size_t(1); i < a1.size(); ++i)
        {
            EXPECT_TRUE(el < data[static_cast<std::size_t>(a1[i])]);
            EXPECT_TRUE(el < data[static_cast<std::size_t>(a2[i])]);
        }
        return true;
    }

    /*py
    a = np.random.randint(0, 1000, size=(20,))
    */
    TEST(xtest_extended, partition)
    {
        // py_a
        xarray<long> py_a = {102,435,860,270,106, 71,700, 20,614,121,466,214,330,458, 87,372, 99,871,
                             663,130};
        
        // py_p5 = np.partition(a, 5)
        xarray<long> py_p5 = { 20, 71, 87, 99,102,106,121,700,614,435,466,214,330,458,270,372,860,871,
                              663,130};
        // py_p0 = np.partition(a, 0)
        xarray<long> py_p0 = { 20,435,860,270,106, 71,700,102,614,121,466,214,330,458, 87,372, 99,871,
                              663,130};
        // py_p13 = np.partition(a, 13)
        xarray<long> py_p13 = { 20,102, 99, 87,106, 71,121,270,130,435,372,214,330,458,614,466,860,871,
                               663,700};
        // py_p19 = np.partition(a, 19)
        xarray<long> py_p19 = { 20,102, 99, 87,106, 71,121,270,130,435,372,214,330,458,663,614,466,700,
                               860,871};
        
        // py_a5 = np.argpartition(a, 5)
        xarray<long> py_a5 = { 7, 5,14,16, 0, 4, 9, 6, 8, 1,10,11,12,13, 3,15, 2,17,18,19};
        // py_a0 = np.argpartition(a, 0)
        xarray<long> py_a0 = { 7, 1, 2, 3, 4, 5, 6, 0, 8, 9,10,11,12,13,14,15,16,17,18,19};
        // py_a13 = np.argpartition(a, 13)
        xarray<long> py_a13 = { 7, 0,16,14, 4, 5, 9, 3,19, 1,15,11,12,13, 8,10, 2,17,18, 6};
        // py_a19 = np.argpartition(a, 19)
        xarray<long> py_a19 = { 7, 0,16,14, 4, 5, 9, 3,19, 1,15,11,12,13,18, 8,10, 6, 2,17};

        auto part_a0 = xt::partition(py_a, 0);

        check_partition_equal(py_p0, part_a0, 0);
        check_partition_equal(py_p5, xt::partition(py_a, 5), 5);
        check_partition_equal(py_p13, xt::partition(py_a, 13), 13);
        check_partition_equal(py_p19, xt::partition(py_a, 19), 19);

        auto parta_a0 = xt::argpartition(py_a, 0);
        check_argpartition_equal(py_a, py_a0, parta_a0, 0);
        check_argpartition_equal(py_a, py_a5, xt::argpartition(py_a, 5), 5);
        check_argpartition_equal(py_a, py_a13, xt::argpartition(py_a, 13), 13);
        check_argpartition_equal(py_a, py_a19, xt::argpartition(py_a, 19), 19);

        // py_median = np.median(a)
        long py_median = 300;
        EXPECT_EQ(py_median, xt::median(py_a));
    }

    /*py
    a = np.random.randint(0, 20, size=(20,))
    */
    TEST(xtest_extended, multi_partition)
    {
        // py_a
        xarray<long> py_a = { 1,11, 5, 1, 0,11,11,16, 9,15,14,14,18,11,19, 2, 4,18, 6, 8};
        
        // py_p0 = np.partition(a, (4, 5, 6))
        xarray<long> py_p0 = { 1, 1, 0, 2, 4, 5, 6, 8, 9,11,14,14,18,11,19,16,11,18,11,15};
        // py_p1 = np.partition(a, (2, 7, 12))
        xarray<long> py_p1 = { 0, 1, 1, 2, 4, 5, 6, 8, 9,11,11,11,11,15,19,16,14,18,18,14};

        auto part_p0 = xt::partition(py_a, {4, 5, 6});
        auto part_p1 = xt::partition(py_a, {2, 7, 12});

        EXPECT_EQ(part_p0(4), py_p0(4));
        EXPECT_EQ(part_p0(5), py_p0(5));
        EXPECT_EQ(part_p0(6), py_p0(6));

        EXPECT_EQ(part_p1(2), py_p1(2));
        EXPECT_EQ(part_p1(7), py_p1(7));
        EXPECT_EQ(part_p1(12), py_p1(12));

        // py_a0 = np.argpartition(a, (4, 5, 6))
        xarray<long> py_a0 = { 0, 3, 4,15,16, 2,18,19, 8, 1,10,11,12,13,14, 7, 6,17, 5, 9};
        // py_a1 = np.argpartition(a, (2, 7, 12))
        xarray<long> py_a1 = { 4, 3, 0,15,16, 2,18,19, 8,13, 1, 6, 5, 9,14, 7,10,17,12,11};

        auto part_a0 = xt::argpartition(py_a, {4, 5, 6});
        auto part_a1 = xt::argpartition(py_a, {2, 7, 12});

        EXPECT_EQ(py_a[part_a0(4)], py_a[static_cast<std::size_t>(py_a0(4))]);
        EXPECT_EQ(py_a[part_a0(5)], py_a[static_cast<std::size_t>(py_a0(5))]);
        EXPECT_EQ(py_a[part_a0(6)], py_a[static_cast<std::size_t>(py_a0(6))]);

        EXPECT_EQ(py_a[part_a1(2)], py_a[static_cast<std::size_t>(py_a1(2))]);
        EXPECT_EQ(py_a[part_a1(7)], py_a[static_cast<std::size_t>(py_a1(7))]);
        EXPECT_EQ(py_a[part_a1(12)], py_a[static_cast<std::size_t>(py_a1(12))]);
    }

    /*py
    a = np.random.rand(5, 5, 5)
    */
    TEST(xtest_extended, axis_median)
    {
        // py_a
        xarray<double> py_a = {{{0.0650515929852795,0.9488855372533332,0.9656320330745594,
                                 0.8083973481164611,0.3046137691733707},
                                {0.0976721140063839,0.6842330265121569,0.4401524937396013,
                                 0.1220382348447788,0.4951769101112702},
                                {0.0343885211152184,0.9093204020787821,0.2587799816000169,
                                 0.662522284353982 ,0.311711076089411 },
                                {0.5200680211778108,0.5467102793432796,0.184854455525527 ,
                                 0.9695846277645586,0.7751328233611146},
                                {0.9394989415641891,0.8948273504276488,0.5978999788110851,
                                 0.9218742350231168,0.0884925020519195}},
                              
                               {{0.1959828624191452,0.0452272889105381,0.3253303307632643,
                                 0.388677289689482 ,0.2713490317738959},
                                {0.8287375091519293,0.3567533266935893,0.2809345096873808,
                                 0.5426960831582485,0.1409242249747626},
                                {0.8021969807540397,0.0745506436797708,0.9868869366005173,
                                 0.7722447692966574,0.1987156815341724},
                                {0.0055221171236024,0.8154614284548342,0.7068573438476171,
                                 0.7290071680409873,0.7712703466859457},
                                {0.0740446517340904,0.3584657285442726,0.1158690595251297,
                                 0.8631034258755935,0.6232981268275579}},
                              
                               {{0.3308980248526492,0.0635583502860236,0.3109823217156622,
                                 0.325183322026747 ,0.7296061783380641},
                                {0.6375574713552131,0.8872127425763265,0.4722149251619493,
                                 0.1195942459383017,0.713244787222995 },
                                {0.7607850486168974,0.5612771975694962,0.770967179954561 ,
                                 0.4937955963643907,0.5227328293819941},
                                {0.4275410183585496,0.0254191267440952,0.1078914269933045,
                                 0.0314291856867343,0.6364104112637804},
                                {0.3143559810763267,0.5085706911647028,0.907566473926093 ,
                                 0.2492922291488749,0.4103829230356297}},
                              
                               {{0.7555511385430487,0.2287981654916225,0.076979909828793 ,
                                 0.289751452913768 ,0.1612212872540044},
                                {0.9296976523425731,0.808120379564417 ,0.6334037565104235,
                                 0.8714605901877177,0.8036720768991145},
                                {0.1865700588860358,0.8925589984899778,0.5393422419156507,
                                 0.8074401551640625,0.8960912999234932},
                                {0.3180034749718639,0.1100519245276768,0.2279351625419417,
                                 0.4271077886262563,0.8180147659224931},
                                {0.8607305832563434,0.0069521305311907,0.5107473025775657,
                                 0.417411003148779 ,0.2221078104707302}},
                              
                               {{0.1198653673336828,0.337615171403628 ,0.9429097039125192,
                                 0.3232029320207552,0.5187906217433661},
                                {0.7030189588951778,0.363629602379294 ,0.9717820827209607,
                                 0.9624472949421112,0.2517822958253642},
                                {0.4972485058923855,0.3008783098167697,0.2848404943774676,
                                 0.0368869473545328,0.6095643339798968},
                                {0.5026790232288615,0.0514787512499894,0.2786464642366114,
                                 0.9082658859666537,0.2395618906669724},
                                {0.1448948720912231,0.489452760277563 ,0.9856504541106007,
                                 0.2420552715115004,0.6721355474058786}}};
        // py_m = np.median(a)
        double py_m = 0.489452760277563;
        
        // py_m0 = np.median(a, 0)
        xarray<double> py_m0 = {{0.1959828624191452,0.2287981654916225,0.3253303307632643,
                                 0.325183322026747 ,0.3046137691733707},
                                {0.7030189588951778,0.6842330265121569,0.4722149251619493,
                                 0.5426960831582485,0.4951769101112702},
                                {0.4972485058923855,0.5612771975694962,0.5393422419156507,
                                 0.662522284353982 ,0.5227328293819941},
                                {0.4275410183585496,0.1100519245276768,0.2279351625419417,
                                 0.7290071680409873,0.7712703466859457},
                                {0.3143559810763267,0.489452760277563 ,0.5978999788110851,
                                 0.417411003148779 ,0.4103829230356297}};
        // py_m1 = np.median(a, 1)
        xarray<double> py_m1 = {{0.0976721140063839,0.8948273504276488,0.4401524937396013,
                                 0.8083973481164611,0.311711076089411 },
                                {0.1959828624191452,0.3567533266935893,0.3253303307632643,
                                 0.7290071680409873,0.2713490317738959},
                                {0.4275410183585496,0.5085706911647028,0.4722149251619493,
                                 0.2492922291488749,0.6364104112637804},
                                {0.7555511385430487,0.2287981654916225,0.5107473025775657,
                                 0.4271077886262563,0.8036720768991145},
                                {0.4972485058923855,0.337615171403628 ,0.9429097039125192,
                                 0.3232029320207552,0.5187906217433661}};
        // py_m2 = np.median(a, 2)
        xarray<double> py_m2 = {{0.8083973481164611,0.4401524937396013,0.311711076089411 ,
                                 0.5467102793432796,0.8948273504276488},
                                {0.2713490317738959,0.3567533266935893,0.7722447692966574,
                                 0.7290071680409873,0.3584657285442726},
                                {0.325183322026747 ,0.6375574713552131,0.5612771975694962,
                                 0.1078914269933045,0.4103829230356297},
                                {0.2287981654916225,0.808120379564417 ,0.8074401551640625,
                                 0.3180034749718639,0.417411003148779 },
                                {0.337615171403628 ,0.7030189588951778,0.3008783098167697,
                                 0.2786464642366114,0.489452760277563 }};

        EXPECT_EQ(py_m, xt::median(py_a));
        EXPECT_EQ(py_m0, xt::median(py_a, 0));
        EXPECT_EQ(py_m1, xt::median(py_a, 1));
        EXPECT_EQ(py_m2, xt::median(py_a, 2));
    }

    /*py
    a = np.random.permutation(np.arange(5 * 5 * 5)).reshape(5, 5, 5)
    */
    TEST(xtest_extended, axis_partition)
    {
        // py_a
        xarray<long> py_a = {{{ 46, 73,109, 61, 71},
                              { 33, 82,  9, 86, 90},
                              { 29, 14,114, 83, 50},
                              {  4, 65, 44, 27, 12},
                              { 63,  5, 57, 25, 84}},
                            
                             {{ 21, 94, 26, 75, 43},
                              { 77,121,124, 99,106},
                              {123, 20, 13, 42,  3},
                              { 64,102, 40, 80, 62},
                              { 78, 60, 45, 17, 28}},
                            
                             {{ 52, 54, 24,107, 34},
                              { 67, 91,120,101,  7},
                              { 48, 87, 93, 98,105},
                              { 53,119,110,104, 76},
                              {116,  8, 31, 36,118}},
                            
                             {{ 11, 47, 88,  1, 59},
                              {112,  6, 41, 30, 22},
                              { 49, 56,113,  0, 81},
                              { 38, 39, 32, 51, 70},
                              { 95,103,100, 66, 89}},
                            
                             {{ 18, 35,108, 19,  2},
                              { 92, 79,117, 58, 72},
                              { 96, 15,111, 10, 85},
                              { 69, 97,115, 68, 23},
                              { 37, 16,122, 55, 74}}};

        // py_p0 = np.partition(a, 2, 0)
        xarray<long> py_p0 = {{{ 11, 35, 24,  1,  2},
                               { 33,  6,  9, 30,  7},
                               { 29, 14, 13,  0,  3},
                               {  4, 39, 32, 27, 12},
                               { 37,  5, 31, 17, 28}},
                             
                              {{ 18, 47, 26, 19, 34},
                               { 67, 79, 41, 58, 22},
                               { 48, 15, 93, 10, 50},
                               { 38, 65, 40, 51, 23},
                               { 63,  8, 45, 25, 74}},
                             
                              {{ 21, 54, 88, 61, 43},
                               { 77, 82,117, 86, 72},
                               { 49, 20,111, 42, 81},
                               { 53, 97, 44, 68, 62},
                               { 78, 16, 57, 36, 84}},
                             
                              {{ 46, 94,109,107, 59},
                               {112, 91,124,101,106},
                               {123, 56,113, 83,105},
                               { 64,102,110, 80, 70},
                               { 95,103,100, 66, 89}},
                             
                              {{ 52, 73,108, 75, 71},
                               { 92,121,120, 99, 90},
                               { 96, 87,114, 98, 85},
                               { 69,119,115,104, 76},
                               {116, 60,122, 55,118}}};
        // py_p1 = np.partition(a, 4, 1)
        xarray<long> py_p1 = {{{  4,  5, 44, 27, 12},
                               { 29, 14, 57, 25, 50},
                               { 33, 65,  9, 61, 71},
                               { 46, 73,109, 83, 84},
                               { 63, 82,114, 86, 90}},
                             
                              {{ 64, 20, 13, 17,  3},
                               { 21, 60, 26, 42, 28},
                               { 77, 94, 40, 75, 43},
                               { 78,102, 45, 80, 62},
                               {123,121,124, 99,106}},
                             
                              {{ 48,  8, 24, 36, 76},
                               { 52, 54, 31, 98, 34},
                               { 53, 87, 93,101,  7},
                               { 67, 91,110,104,105},
                               {116,119,120,107,118}},
                             
                              {{ 38, 39, 32,  0, 70},
                               { 11, 47, 88,  1, 59},
                               { 49,  6, 41, 30, 22},
                               { 95, 56,100, 51, 81},
                               {112,103,113, 66, 89}},
                             
                              {{ 18, 15,108, 10, 23},
                               { 37, 16,111, 19,  2},
                               { 69, 35,115, 55, 72},
                               { 92, 79,117, 58, 74},
                               { 96, 97,122, 68, 85}}};
        // py_p2 = np.partition(a, 3, 2)
        xarray<long> py_p2 = {{{ 61, 46, 71, 73,109},
                               {  9, 33, 82, 86, 90},
                               { 14, 29, 50, 83,114},
                               {  4, 12, 27, 44, 65},
                               { 25, 57,  5, 63, 84}},
                             
                              {{ 21, 26, 43, 75, 94},
                               { 99, 77,106,121,124},
                               {  3, 13, 20, 42,123},
                               { 40, 62, 64, 80,102},
                               { 17, 28, 45, 60, 78}},
                             
                              {{ 24, 34, 52, 54,107},
                               {  7, 67, 91,101,120},
                               { 87, 48, 93, 98,105},
                               { 53, 76,104,110,119},
                               { 36, 31,  8,116,118}},
                             
                              {{  1, 11, 47, 59, 88},
                               { 30, 22,  6, 41,112},
                               {  0, 49, 56, 81,113},
                               { 32, 38, 39, 51, 70},
                               { 66, 89, 95,100,103}},
                             
                              {{  2, 18, 19, 35,108},
                               { 58, 72, 79, 92,117},
                               { 10, 85, 15, 96,111},
                               { 68, 23, 69, 97,115},
                               { 55, 37, 16, 74,122}}};

        auto p0 = xt::partition(py_a, 2, 0);
        auto p1 = xt::partition(py_a, 4, 1);
        auto p2 = xt::partition(py_a, 3, 2);

        EXPECT_EQ(xt::view(py_p0, 2, all(), all()), xt::view(p0, 2, all(), all()));
        EXPECT_EQ(xt::view(py_p1, all(), 4, all()), xt::view(p1, all(), 4, all()));
        EXPECT_EQ(xt::view(py_p2, all(), all(), 3), xt::view(p2, all(), all(), 3));

        // py_a0 = np.argpartition(a, 2, 0)
        xarray<long> py_a0 = {{{3,4,2,3,4},
                               {0,3,0,3,2},
                               {0,0,1,3,1},
                               {0,3,3,0,0},
                               {4,0,2,1,1}},
                             
                              {{4,3,1,4,2},
                               {2,4,3,4,3},
                               {2,4,2,4,0},
                               {3,0,1,3,4},
                               {0,2,1,0,4}},
                             
                              {{1,2,3,0,1},
                               {1,0,4,0,4},
                               {3,1,4,1,3},
                               {2,4,0,4,1},
                               {1,4,0,2,0}},
                             
                              {{0,1,0,2,3},
                               {3,2,1,2,1},
                               {1,3,3,0,2},
                               {1,1,2,1,3},
                               {3,3,3,3,3}},
                             
                              {{2,0,4,1,0},
                               {4,1,2,1,0},
                               {4,2,0,2,4},
                               {4,2,4,2,2},
                               {2,1,4,4,2}}};
        // py_a1 = np.argpartition(a, 4, 1)
        xarray<long> py_a1 = {{{3,4,3,3,3},
                               {2,2,4,4,2},
                               {1,3,1,0,0},
                               {0,0,0,2,4},
                               {4,1,2,1,1}},
                             
                              {{3,2,2,4,2},
                               {0,4,0,2,4},
                               {1,0,3,0,0},
                               {4,3,4,3,3},
                               {2,1,1,1,1}},
                             
                              {{2,4,0,4,3},
                               {0,0,4,2,0},
                               {3,2,2,1,1},
                               {1,1,3,3,2},
                               {4,3,1,0,4}},
                             
                              {{3,3,3,2,3},
                               {0,0,0,0,0},
                               {2,1,1,1,1},
                               {4,2,4,3,2},
                               {1,4,2,4,4}},
                             
                              {{0,2,0,2,3},
                               {4,4,2,0,0},
                               {3,0,3,4,1},
                               {1,1,1,1,4},
                               {2,3,4,3,2}}};
        // py_a2 = np.argpartition(a, 3, 2)
        xarray<long> py_a2 = {{{3,0,4,1,2},
                               {2,0,1,3,4},
                               {1,0,4,3,2},
                               {0,4,3,2,1},
                               {3,2,1,0,4}},
                             
                              {{0,2,4,3,1},
                               {3,0,4,1,2},
                               {4,2,1,3,0},
                               {2,4,0,3,1},
                               {3,4,2,1,0}},
                             
                              {{2,4,0,1,3},
                               {4,0,1,3,2},
                               {1,0,2,3,4},
                               {0,4,3,2,1},
                               {3,2,1,0,4}},
                             
                              {{3,0,1,4,2},
                               {3,4,1,2,0},
                               {3,0,1,4,2},
                               {2,0,1,3,4},
                               {3,4,0,2,1}},
                             
                              {{4,0,3,1,2},
                               {3,4,1,0,2},
                               {3,4,1,0,2},
                               {3,4,0,1,2},
                               {3,0,1,4,2}}};

        auto a0 = xt::argpartition(py_a, 2, 0);
        auto a1 = xt::argpartition(py_a, 4, 1);
        auto a2 = xt::argpartition(py_a, 3, 2);

        EXPECT_EQ(xt::cast<std::size_t>(xt::view(py_a0, 2, all(), all())), xt::view(a0, 2, all(), all()));
        EXPECT_EQ(xt::cast<std::size_t>(xt::view(py_a1, all(), 4, all())), xt::view(a1, all(), 4, all()));
        EXPECT_EQ(xt::cast<std::size_t>(xt::view(py_a2, all(), all(), 3)), xt::view(a2, all(), all(), 3));
    }

    /*py
    a = np.random.permutation(np.arange(5 * 5 * 5)).reshape(5, 5, 5)
    */
    TEST(xtest_extended, multi_k_axis_partition)
    {
        // py_a
        xarray<long> py_a = {{{ 55, 76, 33,103, 60},
                              {  2, 10,  9,100, 67},
                              { 75, 71, 40, 80,108},
                              { 78, 37, 89, 54, 74},
                              {  0, 90, 11,106, 73}},
                            
                             {{ 20, 28, 72, 30, 63},
                              { 81, 27, 97, 88, 18},
                              { 39, 77, 68, 34, 29},
                              { 70,114,  3, 83,  4},
                              {101, 43, 69, 87,119}},
                            
                             {{ 93,123,122, 82, 51},
                              { 48, 58, 61, 99, 79},
                              { 50, 21,109,117, 95},
                              { 91, 14, 13, 15,  7},
                              { 64, 19, 44, 35, 56}},
                            
                             {{  6,111, 12,110,102},
                              { 49, 25, 41, 38, 47},
                              { 42,  8,105, 16,  1},
                              { 26, 65, 22, 85, 46},
                              {107,118,115,120, 57}},
                            
                             {{ 53,113, 24, 17, 66},
                              { 32, 86, 31, 84, 62},
                              { 96, 59,121, 94, 52},
                              {124,112, 92, 23, 36},
                              {104, 98,  5,116, 45}}};

        // py_p0 = np.partition(a, (1, 2), 0)
        xarray<long> py_p0 = {{{  6, 28, 12, 17, 51},
                               {  2, 10,  9, 38, 18},
                               { 39,  8, 40, 16,  1},
                               { 26, 14,  3, 15,  4},
                               {  0, 19,  5, 35, 45}},
                             
                              {{ 20, 76, 24, 30, 60},
                               { 32, 25, 31, 84, 47},
                               { 42, 21, 68, 34, 29},
                               { 70, 37, 13, 23,  7},
                               { 64, 43, 11, 87, 56}},
                             
                              {{ 53,111, 33, 82, 63},
                               { 48, 27, 41, 88, 62},
                               { 50, 59,105, 80, 52},
                               { 78, 65, 22, 54, 36},
                               {101, 90, 44,106, 57}},
                             
                              {{ 55,123,122,110,102},
                               { 49, 58, 61,100, 67},
                               { 75, 71,109,117,108},
                               { 91,114, 89, 85, 46},
                               {107,118,115,120,119}},
                             
                              {{ 93,113, 72,103, 66},
                               { 81, 86, 97, 99, 79},
                               { 96, 77,121, 94, 95},
                               {124,112, 92, 83, 74},
                               {104, 98, 69,116, 73}}};
        // py_p1 = np.partition(a, (1, 4), 1)
        xarray<long> py_p1 = {{{  0, 10,  9, 54, 60},
                               {  2, 37, 11, 80, 67},
                               { 55, 71, 33,100, 73},
                               { 75, 76, 40,103, 74},
                               { 78, 90, 89,106,108}},
                             
                              {{ 20, 27,  3, 30,  4},
                               { 39, 28, 68, 34, 18},
                               { 70, 43, 69, 83, 29},
                               { 81, 77, 72, 87, 63},
                               {101,114, 97, 88,119}},
                             
                              {{ 48, 14, 13, 15,  7},
                               { 50, 19, 44, 35, 51},
                               { 64, 21, 61, 82, 56},
                               { 91, 58,109, 99, 79},
                               { 93,123,122,117, 95}},
                             
                              {{  6,  8, 12, 16,  1},
                               { 26, 25, 22, 38, 46},
                               { 42, 65, 41, 85, 47},
                               { 49,111,105,110, 57},
                               {107,118,115,120,102}},
                             
                              {{ 32, 59,  5, 17, 36},
                               { 53, 86, 24, 23, 45},
                               { 96, 98, 31, 84, 52},
                               {104,112, 92, 94, 62},
                               {124,113,121,116, 66}}};
        // py_p2 = np.partition(a, (1, 3), 2)
        xarray<long> py_p2 = {{{ 33, 55, 60, 76,103},
                               {  2,  9, 10, 67,100},
                               { 40, 71, 75, 80,108},
                               { 37, 54, 74, 78, 89},
                               {  0, 11, 73, 90,106}},
                             
                              {{ 20, 28, 30, 63, 72},
                               { 18, 27, 81, 88, 97},
                               { 29, 34, 39, 68, 77},
                               {  3,  4, 70, 83,114},
                               { 43, 69, 87,101,119}},
                             
                              {{ 51, 82, 93,122,123},
                               { 48, 58, 61, 79, 99},
                               { 21, 50, 95,109,117},
                               {  7, 13, 14, 15, 91},
                               { 19, 35, 44, 56, 64}},
                             
                              {{  6, 12,102,110,111},
                               { 25, 38, 41, 47, 49},
                               {  1,  8, 16, 42,105},
                               { 22, 26, 46, 65, 85},
                               { 57,107,115,118,120}},
                             
                              {{ 17, 24, 53, 66,113},
                               { 31, 32, 62, 84, 86},
                               { 52, 59, 94, 96,121},
                               { 23, 36, 92,112,124},
                               {  5, 45, 98,104,116}}};

        auto p0 = xt::partition(py_a, {1, 2}, 0);
        auto p1 = xt::partition(py_a, {1, 4}, 1);
        auto p2 = xt::partition(py_a, {1, 3}, 2);

        EXPECT_EQ(xt::view(py_p0, 2, all(), all()), xt::view(p0, 2, all(), all()));
        EXPECT_EQ(xt::view(py_p1, all(), 4, all()), xt::view(p1, all(), 4, all()));
        EXPECT_EQ(xt::view(py_p2, all(), all(), 3), xt::view(p2, all(), all(), 3));

        EXPECT_EQ(xt::view(py_p0, 1, all(), all()), xt::view(p0, 1, all(), all()));
        EXPECT_EQ(xt::view(py_p1, all(), 1, all()), xt::view(p1, all(), 1, all()));
        EXPECT_EQ(xt::view(py_p2, all(), all(), 1), xt::view(p2, all(), all(), 1));

        // py_a0 = np.argpartition(a, (1, 2), 0)
        xarray<long> py_a0 = {{{3,1,3,4,2},
                               {0,0,0,3,1},
                               {1,3,0,3,3},
                               {3,2,1,2,1},
                               {0,2,4,2,4}},
                             
                              {{1,0,4,1,0},
                               {4,3,4,4,3},
                               {3,2,1,1,1},
                               {1,0,2,4,2},
                               {2,1,0,1,2}},
                             
                              {{4,3,0,2,1},
                               {2,1,3,1,4},
                               {2,4,3,0,4},
                               {0,3,3,0,4},
                               {1,0,2,0,3}},
                             
                              {{0,2,2,3,3},
                               {3,2,2,0,0},
                               {0,0,2,2,0},
                               {2,1,0,3,3},
                               {3,3,3,3,1}},
                             
                              {{2,4,1,0,4},
                               {1,4,1,2,2},
                               {4,1,4,4,2},
                               {4,4,4,1,0},
                               {4,4,1,4,0}}};
        // py_a1 = np.argpartition(a, (1, 4), 1)
        xarray<long> py_a1 = {{{4,1,1,3,0},
                               {1,3,4,2,1},
                               {0,2,0,1,4},
                               {2,0,2,0,3},
                               {3,4,3,4,2}},
                             
                              {{0,1,3,0,3},
                               {2,0,2,2,1},
                               {3,4,4,3,2},
                               {1,2,0,4,0},
                               {4,3,1,1,4}},
                             
                              {{1,3,3,3,3},
                               {2,4,4,4,0},
                               {4,2,1,0,4},
                               {3,1,2,1,1},
                               {0,0,0,2,2}},
                             
                              {{0,2,0,2,2},
                               {3,1,3,1,3},
                               {2,3,1,3,1},
                               {1,0,2,0,4},
                               {4,4,4,4,0}},
                             
                              {{1,2,4,0,3},
                               {0,1,0,3,4},
                               {2,4,1,1,2},
                               {4,3,3,2,1},
                               {3,0,2,4,0}}};
        // py_a2 = np.argpartition(a, (1, 3), 2)
        xarray<long> py_a2 = {{{2,0,4,1,3},
                               {0,2,1,4,3},
                               {2,1,0,3,4},
                               {1,3,4,0,2},
                               {0,2,4,1,3}},
                             
                              {{0,1,3,4,2},
                               {4,1,0,3,2},
                               {4,3,0,2,1},
                               {2,4,0,3,1},
                               {1,2,3,0,4}},
                             
                              {{4,3,0,2,1},
                               {0,1,2,4,3},
                               {1,0,4,2,3},
                               {4,2,1,3,0},
                               {1,3,2,4,0}},
                             
                              {{0,2,4,3,1},
                               {1,3,2,4,0},
                               {4,1,3,0,2},
                               {2,0,4,1,3},
                               {4,0,2,1,3}},
                             
                              {{3,2,0,4,1},
                               {2,0,4,3,1},
                               {4,1,3,0,2},
                               {3,4,2,1,0},
                               {2,4,1,0,3}}};

        auto a0 = xt::argpartition(py_a, {1, 2}, 0);
        auto a1 = xt::argpartition(py_a, {1, 4}, 1);
        auto a2 = xt::argpartition(py_a, {1, 3}, 2);

        EXPECT_EQ(xt::cast<std::size_t>(xt::view(py_a0, 2, all(), all())), xt::view(a0, 2, all(), all()));
        EXPECT_EQ(xt::cast<std::size_t>(xt::view(py_a1, all(), 4, all())), xt::view(a1, all(), 4, all()));
        EXPECT_EQ(xt::cast<std::size_t>(xt::view(py_a2, all(), all(), 3)), xt::view(a2, all(), all(), 3));

        EXPECT_EQ(xt::cast<std::size_t>(xt::view(py_a0, 1, all(), all())), xt::view(a0, 1, all(), all()));
        EXPECT_EQ(xt::cast<std::size_t>(xt::view(py_a1, all(), 1, all())), xt::view(a1, all(), 1, all()));
        EXPECT_EQ(xt::cast<std::size_t>(xt::view(py_a2, all(), all(), 1)), xt::view(a2, all(), all(), 1));
    }
}

