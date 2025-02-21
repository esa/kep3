// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <heyoka/taylor.hpp>
#include <vector>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/ta/pontryagin_cartesian.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using heyoka::taylor_adaptive;
using heyoka::taylor_outcome;

using kep3::ta::get_ta_pc;
using kep3::ta::get_ta_pc_cache_dim;
using kep3::ta::get_ta_pc_var;
using kep3::ta::get_ta_pc_var_cache_dim;

// This needs to be the first test as cache dimension will be assumed to be zero here.
TEST_CASE("caches")
{
    {
        // MASS
        //  The non variational one.
        REQUIRE(get_ta_pc_cache_dim() == 0u);
        auto ta_cached = get_ta_pc(1e-16, kep3::optimality_type::MASS);
        REQUIRE(get_ta_pc_cache_dim() == 1u);
        ta_cached = get_ta_pc(1e-16, kep3::optimality_type::MASS);
        REQUIRE(get_ta_pc_cache_dim() == 1u);
        ta_cached = get_ta_pc(1e-8, kep3::optimality_type::MASS);
        REQUIRE(get_ta_pc_cache_dim() == 2u);

        // The variational integrator.
        REQUIRE(get_ta_pc_var_cache_dim() == 0u);
        ta_cached = get_ta_pc_var(1e-16, kep3::optimality_type::MASS);
        REQUIRE(get_ta_pc_var_cache_dim() == 1u);
        ta_cached = get_ta_pc_var(1e-16, kep3::optimality_type::MASS);
        REQUIRE(get_ta_pc_var_cache_dim() == 1u);
        ta_cached = get_ta_pc_var(1e-8, kep3::optimality_type::MASS);
        REQUIRE(get_ta_pc_var_cache_dim() == 2u);
    }
    {
        // TIME
        //  The non variational one. (no cache is not empty)
        REQUIRE(get_ta_pc_cache_dim() == 2u);
        auto ta_cached = get_ta_pc(1e-16, kep3::optimality_type::TIME);
        REQUIRE(get_ta_pc_cache_dim() == 3u);
        ta_cached = get_ta_pc(1e-16, kep3::optimality_type::TIME);
        REQUIRE(get_ta_pc_cache_dim() == 3u);
        ta_cached = get_ta_pc(1e-8, kep3::optimality_type::TIME);
        REQUIRE(get_ta_pc_cache_dim() == 4u);

        // The variational integrator.
        REQUIRE(get_ta_pc_var_cache_dim() == 2u);
        ta_cached = get_ta_pc_var(1e-16, kep3::optimality_type::TIME);
        REQUIRE(get_ta_pc_var_cache_dim() == 3u);
        ta_cached = get_ta_pc_var(1e-16, kep3::optimality_type::TIME);
        REQUIRE(get_ta_pc_var_cache_dim() == 3u);
        ta_cached = get_ta_pc_var(1e-8, kep3::optimality_type::TIME);
        REQUIRE(get_ta_pc_var_cache_dim() == 4u);
    }
}

TEST_CASE("dynamics_mass")
{
    auto ta_cached = get_ta_pc(1e-16, kep3::optimality_type::MASS);
    REQUIRE(ta_cached.is_variational() == false);
    REQUIRE(ta_cached.get_dim() == 14);
    REQUIRE(ta_cached.get_pars().size() == 5); // [mu, c1, c2, eps, l0]

    {
        // We test a meaningless case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        std::vector<double> ic = {1., 0., 0., 0., 1., 0., 10., 1., 1., 1., 1., 1., 1., 1.};
        std::vector<double> pars = {1., 0.01, 1., 0.5, 1.};
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        auto out = ta.propagate_until(1.2345);
        std::vector<double> const ground_truth
            = {3.29796945478044e-01,  9.43846387947806e-01,  -6.98808299448645e-05, -9.44346244116373e-01,
               3.29727154212398e-01,  -1.65088454784805e-05, 9.99636411645986e+00,  -2.93037901256056e-01,
               1.61477939441094e-01,  1.27398659084316e+00,  9.47167118849066e-01,  1.85641413513971e-02,
               -6.14020501980864e-01, 9.99960774095174e-01};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-13);
    }
}

TEST_CASE("variational_dynamics_mass")
{
    auto ta_cached = get_ta_pc_var(1e-16, kep3::optimality_type::MASS);
    REQUIRE(ta_cached.is_variational() == true);
    REQUIRE(ta_cached.get_dim() == 14 + 14 * 8);
    REQUIRE(ta_cached.get_pars().size() == 5);
    REQUIRE(ta_cached.get_vorder() == 1);
    REQUIRE(ta_cached.get_vargs().size() == 8);

    {
        // We test a meaningless case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        std::vector<double> ic = {1., 0., 0., 0., 1., 0., 10., 1., 1., 1., 1., 1., 1., 1.};
        std::vector<double> pars = {1., 0.01, 1., 0.5, 1.};
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        auto out = ta.propagate_until(1.2345);
        std::vector<double> const ground_truth
            = {3.2979694547804350e-01,  9.4384638794780595e-01,  -6.9880829944868176e-05, -9.4434624411637313e-01,
               3.2972715421239757e-01,  -1.6508845478493283e-05, 9.9963641164598567e+00,  -2.9303790125605605e-01,
               1.6147793944109437e-01,  1.2739865908431554e+00,  9.4716711884906613e-01,  1.8564141351397081e-02,
               -6.1402050198086411e-01, 9.9996077409517425e-01,  2.4807113499712627e-05,  -1.8351141891890819e-05,
               -1.2818470510930369e-05, -8.5862816895937083e-05, 3.8310761104761820e-05,  5.3798677075455241e-05,
               -9.8439713220596816e-07, 1.1002747510345603e-06,  -1.8116386426256103e-05, 7.8320073528235891e-05,
               -1.3039975755904371e-05, 3.9559384976377910e-05,  -1.3397644584935795e-04, 4.7166350895785629e-05,
               -6.8702394931855638e-07, 7.7402258043757117e-07,  -7.7410298349099223e-06, -7.8270913104314583e-06,
               7.1927187950894463e-05,  3.6928270174194258e-05,  2.8940566765850425e-05,  -1.2227435927319579e-04,
               -3.5056346905390064e-07, 3.9701899665195141e-07,  4.3025732061979052e-05,  -1.1307741359621161e-05,
               5.9755824209336541e-06,  -1.0323690710348835e-04, 4.6713450168555002e-06,  6.0671057095680562e-05,
               -1.8208045330794651e-06, 2.0217364007402280e-06,  -1.0051867107766068e-05, 2.2223999075462209e-04,
               -2.3387904095795131e-05, -7.9450488983894245e-06, -2.6559987642840456e-04, 8.4581885146911274e-05,
               -1.3346761093572137e-06, 1.4974967381791003e-06,  3.3217922956180182e-05,  2.3278281398961398e-06,
               1.6043921032697509e-04,  -3.1117753363833245e-05, 5.8632246493092551e-06,  -1.7074415278991729e-04,
               -8.2718677759161214e-08, 9.6438759149039954e-08,  1.1803628114620423e-06,  3.9436827275116442e-07,
               -3.7470751549576763e-08, -2.2707173669061929e-06, -8.7535295896060392e-07, -3.5815361574013834e-07,
               -1.8238267043205607e-05, 2.0205230652148911e-05,  1.4421278777367776e+00,  9.2031553535555877e-01,
               -4.7732704896099501e-05, -1.5430545138592504e+00, -1.1124361873101991e+00, 5.7099713365790765e-05,
               -7.5016034151555402e-08, 9.4828622184191523e-08,  1.2653497327417409e+00,  1.6020737674955388e+00,
               -8.9798078791050680e-05, -2.3845927734959336e+00, -3.2140431397897290e-01, 1.4112944674914900e-04,
               -1.4871350129381810e-06, 1.6824457767373286e-06,  -7.6364027374917128e-05, -6.4405992370654094e-05,
               3.2977276917101150e-01,  1.2010714145562803e-04,  2.6834074244918730e-05,  9.4420766161851133e-01,
               5.0733780373013986e-08,  -6.1876102741464402e-08, -1.5764794684511851e+00, -3.7389609717733641e-01,
               2.2462477322421727e-05,  2.2649940942886935e+00,  6.3257564501327501e-01,  -4.9483316938361257e-05,
               2.3041257478972935e-07,  -2.6439734026273247e-07, -4.4901074031759292e-01, -1.2889175884820021e+00,
               2.5840943825693690e-05,  9.7755277936607832e-01,  7.7896421348181277e-01,  -5.0311367430719056e-05,
               3.9652291042017644e-07,  -4.4879620481642567e-07, 2.2259338708594027e-05,  1.6905367834749667e-05,
               -9.4392585514698957e-01, -3.8682769883439057e-05, -1.0143261673006455e-05, 3.2991499204769009e-01,
               -1.5221853007620489e-07, 1.7466197853629203e-07,  2.3551763027498985e-05,  7.8688806302431488e-06,
               -7.4760229283583263e-07, -4.5308065920295875e-05, -1.7466474729907059e-05, -7.1468319961342534e-06,
               9.9999980316480919e-01,  2.1926164655457741e-07};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-13);
    }
}

TEST_CASE("dynamics_time")
{
    auto ta_cached = get_ta_pc(1e-16, kep3::optimality_type::TIME);
    REQUIRE(ta_cached.is_variational() == false);
    REQUIRE(ta_cached.get_dim() == 14);
    REQUIRE(ta_cached.get_pars().size() == 3); // [mu, c1, c2, eps, l0]

    {
        // We test a meaningless case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        std::vector<double> ic = {1., 0., 0., 0., 1., 0., 10., 1., 1., 1., 1., 1., 1., 1.};
        std::vector<double> pars = {1., 0.01, 1.};
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        auto out = ta.propagate_until(1.2345);
        std::vector<double> const ground_truth
            = {3.29326772571214e-01,  9.43518279694245e-01,  -2.37243691077672e-04, -9.45216052816242e-01,
               3.29089569819722e-01,  -5.58634286200389e-05, 9.98765500000000e+00,  -2.93073464383497e-01,
               1.60767429870367e-01,  1.27401079872040e+00,  9.47277110771498e-01,  1.87535673100874e-02,
               -6.14093196354726e-01, 9.99866710251057e-01};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-13);
    }
}

TEST_CASE("variational_dynamics_time")
{
    auto ta_cached = get_ta_pc_var(1e-16, kep3::optimality_type::TIME);
    REQUIRE(ta_cached.is_variational() == true);
    REQUIRE(ta_cached.get_dim() == 14 + 14 * 8);
    REQUIRE(ta_cached.get_pars().size() == 3);
    REQUIRE(ta_cached.get_vorder() == 1);
    REQUIRE(ta_cached.get_vargs().size() == 8);

    {
        // We test a meaningless case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        std::vector<double> ic = {1., 0., 0., 0., 1., 0., 10., 1., 1., 1., 1., 1., 1., 1.};
        std::vector<double> pars = {1., 0.01, 1.};
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        auto out = ta.propagate_until(1.2345);
        std::vector<double> const ground_truth
            = {3.29326772571214e-01,  9.43518279694245e-01,  -2.37243691077678e-04, -9.45216052816242e-01,
               3.29089569819723e-01,  -5.58634286200446e-05, 9.98765500000000e+00,  -2.93073464383497e-01,
               1.60767429870367e-01,  1.27401079872040e+00,  9.47277110771498e-01,  1.87535673100874e-02,
               -6.14093196354726e-01, 9.99866710251056e-01,  8.41133971750316e-05,  -6.23805339365836e-05,
               -4.35576346536962e-05, -2.91230168218841e-04, 1.30236152736551e-04,  1.82818786897539e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  -6.15912904833284e-05, 2.66032245545687e-04,
               -4.43225236984674e-05, 1.34507269522455e-04,  -4.54930989076605e-04, 1.60305288190260e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  -2.62904019617088e-05, -2.65922554162050e-05,
               2.44305289111464e-04,  1.25451303614909e-04,  9.83231934912956e-05,  -4.15197128839756e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  1.45747427272735e-04,  -3.84863004523546e-05,
               2.03218553225493e-05,  -3.49785695279272e-04, 1.60038710912614e-05,  2.06198842045081e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  -3.42763969439676e-05, 7.55063066434171e-04,
               -7.95411628295219e-05, -2.67105446050315e-05, -9.02168112273193e-04, 2.87633150217543e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  1.12967971828034e-04,  7.93937073268113e-06,
               5.44943821249505e-04,  -1.05916110655104e-04, 1.98653792978839e-05,  -5.79800432453000e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,
               0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,
               0.00000000000000e+00,  0.00000000000000e+00,  1.44196561049032e+00,  9.20877930716147e-01,
               -1.62136878165076e-04, -1.54295185938881e+00, -1.11299681577074e+00, 1.93806447760147e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  1.26630644933279e+00,  1.60313599146158e+00,
               -3.05364032176103e-04, -2.38639680762912e+00, -3.22452855839271e-01, 4.80016576560558e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  -2.59565205932807e-04, -2.18973016400636e-04,
               3.29244288341037e-01,  4.08237064294930e-04,  9.12183025923601e-05,  9.44745593234806e-01,
               0.00000000000000e+00,  0.00000000000000e+00,  -1.57652935813679e+00, -3.74079622523328e-01,
               7.63176133867316e-05,  2.26516174045298e+00,  6.32816141123617e-01,  -1.68107758370038e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  -4.49254170142732e-01, -1.28907236682424e+00,
               8.78481189862216e-05,  9.78050639895916e-01,  7.79112671238256e-01,  -1.71054976102050e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  7.56303455884195e-05,  5.74555043756629e-05,
               -9.43788098497564e-01, -1.31402239732434e-04, -3.44548398801908e-05, 3.29727673372486e-01,
               0.00000000000000e+00,  0.00000000000000e+00,  8.00334644311539e-05,  2.67404597714927e-05,
               -2.55263719046308e-06, -1.53938341730639e-04, -5.93271855549472e-05, -2.42455086700587e-05,
               1.00000000000000e+00,  0.00000000000000e+00};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-12);
    }
}