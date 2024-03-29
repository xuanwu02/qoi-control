#ifndef _MDR_GE_QOI_WIZARD_HPP
#define _MDR_GE_QOI_WIZARD_HPP

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>
#include <cmath>
#include <numeric>
// #include "QoIsWizard/ge_LpSolver.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp"
#include <xtensor/xindex_view.hpp>

namespace MDR {

template<class Type>
Type xt_ptp(const xt::xarray<Type>& arr){
    return xt::amax(arr)() - xt::amin(arr)();
}

template<class Type>
double rel_L_inf(const xt::xarray<Type>& f_ori, const xt::xarray<Type>& f_dec){
    double range_f_ori = xt_ptp(f_ori);
    double error = xt::amax(xt::abs(f_ori - f_dec))();
    return error / range_f_ori;
}

template<class Type>
Type vector_range(const std::vector<Type>& vec){
    auto [min_it, max_it] = std::minmax_element(vec.begin(), vec.end());
    Type vec_range = *max_it - *min_it;
    return vec_range;
}

template<class Type>
xt::xarray<Type> ptr2xarray(const Type * data_ptr, size_t length){
    xt::xarray<double> data_array = xt::adapt(data_ptr, length, xt::no_ownership(), std::vector<std::size_t>{length});
    return data_array;
}

class GE_Var_ebList{
public:
    double eb_pres;
    double eb_den;
    double eb_vx;
    double eb_vy;
    double eb_vz;

    std::vector<double> vec;
    
    GE_Var_ebList() = default;

    GE_Var_ebList(double e_p, double e_d,
                  double e_x, double e_y, double e_z)
    : eb_pres(e_p), eb_den(e_d), eb_vx(e_x), eb_vy(e_y), eb_vz(e_z){
        vec = {e_p, e_d, e_x, e_y, e_z};
    }

    void update_all(double e_p, double e_d,
                  double e_x, double e_y, double e_z){
        eb_pres = e_p;
        eb_den = e_d;
        eb_vx = e_x;
        eb_vy = e_y;
        eb_vz = e_z;
        vec = {e_p, e_d, e_x, e_y, e_z};
    }
};

struct GE_QoI_errList{
    double VTOT;
    double T;
    double C;
    double M;
    double mu;
    double PT;

    GE_QoI_errList(double err_V, double err_T, double err_C, double err_M, double err_mu, double err_PT)
    : VTOT(err_V), T(err_T), C(err_C), M(err_M), mu(err_mu), PT(err_PT){}
};

template<class Type>
class GEQoI_gen{
public:

    const double gamma = 1.4;
    const double R = 287.1;
    const double mu_r = 1.716e-5;
    const double Tr = 273.15;
    const double S = 110.4;
    const int num_variables = 5;

    size_t num_elements;

    xt::xarray<Type> pressure, density, velocityX, velocityY, velocityZ;
    xt::xarray<Type> VTOT, T, C, M, mu, PT;

    GEQoI_gen(const std::vector<Type>& pres, const std::vector<Type>& den, 
              const std::vector<Type>& vx, const std::vector<Type>& vy, 
              const std::vector<Type>& vz, size_t num)
    : num_elements(num), pressure(xt::xarray<Type>(xt::adapt(pres, {num}))), density(xt::xarray<Type>(xt::adapt(den, {num}))),
      velocityX(xt::xarray<Type>(xt::adapt(vx, {num}))), velocityY(xt::xarray<Type>(xt::adapt(vy, {num}))), velocityZ(xt::xarray<Type>(xt::adapt(vz, {num}))){
        update_QoIs(pressure, density, velocityX, velocityY, velocityZ, VTOT, T, C, M, mu, PT);
    }

    void update_QoIs(const xt::xarray<Type>& pres, const xt::xarray<Type>& den, 
                     const xt::xarray<Type>& vx, const xt::xarray<Type>& vy, 
                     const xt::xarray<Type>& vz,
                     xt::xarray<Type>& VTOT_ref, xt::xarray<Type>& T_ref, 
                     xt::xarray<Type>& C_ref, xt::xarray<Type>& M_ref, 
                     xt::xarray<Type>& mu_ref, xt::xarray<Type>& PT_ref){
        VTOT_ref = compute_VTOT(vx, vy, vz);
        T_ref = compute_T(pres, den);
        C_ref = compute_C(pres, den);
        M_ref = compute_M(VTOT_ref, C_ref);
        mu_ref = compute_mu(pres, den);
        PT_ref = compute_PT(pres, M_ref);
    }

    void update_State(const std::vector<Type>& pres, const std::vector<Type>& den, 
                        const std::vector<Type>& vx, const std::vector<Type>& vy, 
                        const std::vector<Type>& vz, size_t num){
        num_elements = num;
        pressure = xt::adapt(pres, {num});
        density = xt::adapt(den, {num});
        velocityX = xt::adapt(vx, {num});
        velocityY = xt::adapt(vy, {num});
        velocityZ = xt::adapt(vz, {num});
        update_QoIs(pressure, density, velocityX, velocityY, velocityZ, VTOT, T, C, M, mu, PT);
    }

    GE_QoI_errList compute_QoI_error(const std::vector<Type>& dec_pres, const std::vector<Type>& dec_den, 
                                     const std::vector<Type>& dec_vx, const std::vector<Type>& dec_vy, 
                                     const std::vector<Type>& dec_vz){
        xt::xarray<Type> dec_pressure, dec_density, dec_velocityX, dec_velocityY, dec_velocityZ;
        dec_pressure = xt::adapt(dec_pres, {num_elements});
        dec_density = xt::adapt(dec_den, {num_elements});
        dec_velocityX = xt::adapt(dec_vx, {num_elements});
        dec_velocityY = xt::adapt(dec_vy, {num_elements});
        dec_velocityZ = xt::adapt(dec_vz, {num_elements});

        xt::xarray<Type> dec_VTOT, dec_T, dec_C, dec_M, dec_mu, dec_PT;
        update_QoIs(dec_pressure, dec_density, dec_velocityX, dec_velocityY, dec_velocityZ, dec_VTOT, dec_T, dec_C, dec_M, dec_mu, dec_PT);

        return GE_QoI_errList(
            rel_L_inf(VTOT, dec_VTOT),
            rel_L_inf(T, dec_T),
            rel_L_inf(C, dec_C),
            rel_L_inf(M, dec_M),
            rel_L_inf(mu, dec_mu),
            rel_L_inf(PT, dec_PT)
        );
    }

    GE_QoI_errList compute_QoI_error(const std::vector<const Type*>& dec_vars){
        assert(dec_vars.size() == num_variables);
        std::vector<xt::xarray<Type>> dec_arrays(num_variables);
        for(int i=0; i< num_variables; i++){
            dec_arrays[i] = ptr2xarray(dec_vars[i], num_elements);
        }

        xt::xarray<Type> dec_VTOT, dec_T, dec_C, dec_M, dec_mu, dec_PT;
        update_QoIs(dec_arrays[0], dec_arrays[1], dec_arrays[2], dec_arrays[3], dec_arrays[4], dec_VTOT, dec_T, dec_C, dec_M, dec_mu, dec_PT);

        return GE_QoI_errList(
            rel_L_inf(VTOT, dec_VTOT),
            rel_L_inf(T, dec_T),
            rel_L_inf(C, dec_C),
            rel_L_inf(M, dec_M),
            rel_L_inf(mu, dec_mu),
            rel_L_inf(PT, dec_PT)
        );
    }

    // compute maximum absolute value of partial derivatives
    // VTOT
    Type max_dV_dvelocityX() const {
        auto mask = xt::not_equal(VTOT, 0);
        xt::xarray<Type> non_zero_V = xt::filter(VTOT, mask);
        xt::xarray<Type> non_zero_V_vx = xt::filter(velocityX, mask);
        xt::xarray<Type> dV_dvx = non_zero_V_vx / non_zero_V;
        Type res = xt::amax(xt::abs(dV_dvx))();
        return res;
    }
    Type max_dV_dvelocityY() const {
        auto mask = xt::not_equal(VTOT, 0);
        xt::xarray<Type> non_zero_V = xt::filter(VTOT, mask);
        xt::xarray<Type> non_zero_V_vy = xt::filter(velocityY, mask);
        xt::xarray<Type> dV_dvy = non_zero_V_vy / non_zero_V;
        Type res = xt::amax(xt::abs(dV_dvy))();
        return res;
    }
    Type max_dV_dvelocityZ() const {
        auto mask = xt::not_equal(VTOT, 0);
        xt::xarray<Type> non_zero_V = xt::filter(VTOT, mask);
        xt::xarray<Type> non_zero_V_vz = xt::filter(velocityZ, mask);
        xt::xarray<Type> dV_dvz = non_zero_V_vz / non_zero_V;
        Type res = xt::amax(xt::abs(dV_dvz))();
        return res;
    }
    // T
    Type max_dT_dpressure() const {
        xt::xarray<Type> dT_dpres = 1/R / density;
        Type res = xt::amax(xt::abs(dT_dpres))();
        return res;
    }
    Type max_dT_ddensity() const {
        xt::xarray<Type> dT_dden = -1 * pressure/R / xt::square(density);
        Type res = xt::amax(xt::abs(dT_dden))();
        return res;
    }
    // C
    Type max_dC_dpressure() const {
        xt::xarray<Type> dC_dpres = sqrt(gamma)/2 / xt::sqrt(pressure*density);
        Type res = xt::amax(xt::abs(dC_dpres))();
        return res;
    }
    Type max_dC_ddensity() const {
        xt::xarray<Type> dC_dden = pressure/density/2 / xt::sqrt(pressure*density);
        Type res = xt::amax(xt::abs(dC_dden))();
        return res;
    }
    // mu
    Type max_dmu_dpressure() const {
        double cc = mu_r*(Tr+S) / pow(Tr, 1.5);
        xt::xarray<Type> dmu_dpres = cc * (0.5*sqrt(R)*xt::pow(pressure,1.5)*xt::sqrt(density) + 1.5*pow(R,1.5)*S*xt::sqrt(pressure)*xt::pow(density,1.5))
                                        / xt::square(sqrt(R)*pressure*xt::sqrt(density)+pow(R,1.5)*S*xt::pow(density,1.5));
        Type res = xt::amax(xt::abs(dmu_dpres))();
        return res;
    }
    Type max_dmu_ddensity() const {
        double cc = -mu_r*(Tr+S) / pow(Tr, 1.5);
        xt::xarray<Type> dmu_dden = cc * pressure/xt::square(density)/R * (0.5*xt::pow(T,1.5) + 1.5*xt::pow(T,0.5)*S) / xt::square(T+S);
        Type res = xt::amax(xt::abs(dmu_dden))();
        return res;
    }
    // Mach
    Type max_dM_dpressure() const {
        xt::xarray<Type> dM_dpres = -(7*xt::sqrt(xt::square(velocityX) + xt::square(velocityY) + xt::square(velocityZ)))/(10*density*xt::pow((7*pressure)/(5*density),1.5));
        Type res = xt::amax(xt::abs(dM_dpres))();
        return res;
    }
    Type max_dM_ddensity() const {
        xt::xarray<Type> dM_dden = (7*pressure*xt::sqrt(xt::square(velocityX) + xt::square(velocityY) + xt::square(velocityZ)))/(10*xt::square(density)*xt::pow((7*pressure)/(5*density),1.5));
        Type res = xt::amax(xt::abs(dM_dden))();
        return res;
    }
    Type max_dM_dvelocityX() const {
        xt::xarray<Type> dM_dv = xt::sqrt(5*density/(7*pressure));
        Type res = xt::amax(xt::abs(dM_dv))();
        return res;
    }
    Type max_dM_dvelocityY() const {
        xt::xarray<Type> dM_dv = xt::sqrt(5*density/(7*pressure));
        Type res = xt::amax(xt::abs(dM_dv))();
        return res;
    }
    Type max_dM_dvelocityZ() const {
        xt::xarray<Type> dM_dv = xt::sqrt(5*density/(7*pressure));
        Type res = xt::amax(xt::abs(dM_dv))();
        return res;
    }
    // PT
    xt::xarray<Type> compute_VTOT_squared() const {
        return xt::square(velocityX) + xt::square(velocityY) + xt::square(velocityZ);
    }
    Type max_dPT_dpressure() const {
        xt::xarray<Type> VTOT_squared = compute_VTOT_squared();
        xt::xarray<Type> dPT_dpres = xt::pow((density*(VTOT_squared))/(7*pressure) + 1,3.5) - (density * xt::pow((density*(VTOT_squared))/(7*pressure) + 1,2.5)*(VTOT_squared))/(2*pressure);
        Type res = xt::amax(xt::abs(dPT_dpres))();
        return res;
    }   
    Type max_dPT_ddensity() const {
        xt::xarray<Type> VTOT_squared = compute_VTOT_squared();
        xt::xarray<Type> dPT_dden = 0.5 * xt::pow((density*(VTOT_squared))/(7*pressure) + 1,2.5)*(VTOT_squared);
        Type res = xt::amax(xt::abs(dPT_dden))();
        return res;
    } 
    Type max_dPT_dvelocityX() const {
        xt::xarray<Type> dPT_dvx = density * velocityX * xt::pow((density*compute_VTOT_squared())/(7*pressure) + 1,2.5);
        Type res = xt::amax(xt::abs(dPT_dvx))();
        return res;
    } 
    Type max_dPT_dvelocityY() const {
        xt::xarray<Type> dPT_dvy = density * velocityY * xt::pow((density*compute_VTOT_squared())/(7*pressure) + 1,2.5);
        Type res = xt::amax(xt::abs(dPT_dvy))();
        return res;
    }
    Type max_dPT_dvelocityZ() const {
        xt::xarray<Type> dPT_dvz = density * velocityZ * xt::pow((density*compute_VTOT_squared())/(7*pressure) + 1,2.5);
        Type res = xt::amax(xt::abs(dPT_dvz))();
        return res;
    }

    // LP generic error bound derivation method is in "ge_LpSolver.hpp"

    // LP top-2-bottom
    xt::xarray<double> Lp_VTOT_2_vx_vy_vz(double e_f) const {
        std::string lp_type = "VTOT";
        int n_row = 1, n_col = 3;

        std::vector<xt::xarray<Type>> coeffs(n_row);
        coeffs[0] = {max_dV_dvelocityX() * xt_ptp(velocityX), max_dV_dvelocityY() * xt_ptp(velocityY), max_dV_dvelocityZ() * xt_ptp(velocityZ)};

        xt::xarray<Type> bounds = e_f * xt_ptp(VTOT);
        xt::xarray<Type> range_args = {xt_ptp(velocityX), xt_ptp(velocityY), xt_ptp(velocityZ)};

        double optimal = 0;
        auto result = generic_LpDerivation(lp_type, n_row, n_col, e_f, coeffs, bounds, range_args, optimal);
        return result;
    }
    xt::xarray<double> Lp_T_2_pressure_density(double e_f) const {
        std::string lp_type = "T";
        int n_row = 1, n_col = 2;

        std::vector<xt::xarray<Type>> coeffs(n_row);
        coeffs[0] = {max_dT_dpressure() * xt_ptp(pressure), max_dT_ddensity() * xt_ptp(density)};

        xt::xarray<Type> bounds = e_f * xt_ptp(T);
        xt::xarray<Type> range_args = {xt_ptp(pressure), xt_ptp(density)};

        double optimal = 0;
        auto result = generic_LpDerivation(lp_type, n_row, n_col, e_f, coeffs, bounds, range_args, optimal);
        return result;
    }
    xt::xarray<double> Lp_C_2_pressure_density(double e_f) const {
        std::string lp_type = "C";
        int n_row = 1, n_col = 2;

        std::vector<xt::xarray<Type>> coeffs(n_row);
        coeffs[0] = {max_dC_dpressure() * xt_ptp(pressure), max_dC_ddensity() * xt_ptp(density)};

        xt::xarray<Type> bounds = e_f * xt_ptp(C);
        xt::xarray<Type> range_args = {xt_ptp(pressure), xt_ptp(density)};

        double optimal = 0;
        auto result = generic_LpDerivation(lp_type, n_row, n_col, e_f, coeffs, bounds, range_args, optimal);
        return result;
    }
    xt::xarray<double> Lp_mu_2_pressure_density(double e_f) const {
        std::string lp_type = "mu";
        int n_row = 1, n_col = 2;

        std::vector<xt::xarray<Type>> coeffs(n_row);
        coeffs[0] = {max_dmu_dpressure() * xt_ptp(pressure), max_dmu_ddensity() * xt_ptp(density)};

        xt::xarray<Type> bounds = e_f * xt_ptp(mu);
        xt::xarray<Type> range_args = {xt_ptp(pressure), xt_ptp(density)};

        double optimal = 0;
        auto result = generic_LpDerivation(lp_type, n_row, n_col, e_f, coeffs, bounds, range_args, optimal);
        return result;
    }
    xt::xarray<double> Lp_M_2_all(double e_f) const {
        std::string lp_type = "M";
        int n_row = 1, n_col = 5;

        std::vector<xt::xarray<Type>> coeffs(n_row);
        coeffs[0] = {max_dM_dpressure() * xt_ptp(pressure), max_dM_ddensity() * xt_ptp(density), max_dM_dvelocityX() * xt_ptp(velocityX), max_dM_dvelocityY() * xt_ptp(velocityY), max_dM_dvelocityZ() * xt_ptp(velocityZ)};

        xt::xarray<Type> bounds = e_f * xt_ptp(M);
        xt::xarray<Type> range_args = {xt_ptp(pressure), xt_ptp(density), xt_ptp(velocityX), xt_ptp(velocityY), xt_ptp(velocityZ)};

        double optimal = 0;
        auto result = generic_LpDerivation(lp_type, n_row, n_col, e_f, coeffs, bounds, range_args, optimal);
        return result;
    }
    xt::xarray<double> Lp_PT_2_all(double e_f) const {
        std::string lp_type = "PT";
        int n_row = 1, n_col = 5;

        std::vector<xt::xarray<Type>> coeffs(n_row);
        coeffs[0] = {max_dPT_dpressure() * xt_ptp(pressure), max_dPT_ddensity() * xt_ptp(density), max_dPT_dvelocityX() * xt_ptp(velocityX), max_dPT_dvelocityY() * xt_ptp(velocityY), max_dPT_dvelocityZ() * xt_ptp(velocityZ)};

        xt::xarray<Type> bounds = e_f * xt_ptp(PT);
        xt::xarray<Type> range_args = {xt_ptp(pressure), xt_ptp(density), xt_ptp(velocityX), xt_ptp(velocityY), xt_ptp(velocityZ)};

        double optimal = 0;
        auto result = generic_LpDerivation(lp_type, n_row, n_col, e_f, coeffs, bounds, range_args, optimal);
        return result;
    }
    GE_Var_ebList Lp_all_2_all(double e_f) const {
        std::string lp_type = "all";
        int n_row = 6, n_col = 5;

        std::vector<xt::xarray<Type>> coeffs(n_row);
        coeffs[0] = {0, 0, max_dV_dvelocityX() * xt_ptp(velocityX), max_dV_dvelocityY() * xt_ptp(velocityY), max_dV_dvelocityZ() * xt_ptp(velocityZ)};
        coeffs[1] = {max_dT_dpressure() * xt_ptp(pressure), max_dT_ddensity() * xt_ptp(density), 0, 0, 0};
        coeffs[2] = {max_dC_dpressure() * xt_ptp(pressure), max_dC_ddensity() * xt_ptp(density), 0, 0, 0};
        coeffs[3] = {max_dmu_dpressure() * xt_ptp(pressure), max_dmu_ddensity() * xt_ptp(density), 0, 0, 0};
        coeffs[4] = {max_dM_dpressure() * xt_ptp(pressure), max_dM_ddensity() * xt_ptp(density), max_dM_dvelocityX() * xt_ptp(velocityX), max_dM_dvelocityY() * xt_ptp(velocityY), max_dM_dvelocityZ() * xt_ptp(velocityZ)};
        coeffs[5] = {max_dPT_dpressure() * xt_ptp(pressure), max_dPT_ddensity() * xt_ptp(density), max_dPT_dvelocityX() * xt_ptp(velocityX), max_dPT_dvelocityY() * xt_ptp(velocityY), max_dPT_dvelocityZ() * xt_ptp(velocityZ)};

        xt::xarray<Type> bounds = {xt_ptp(VTOT), xt_ptp(T), xt_ptp(C), xt_ptp(mu), xt_ptp(M), xt_ptp(PT)};
        bounds *= e_f;
        xt::xarray<Type> range_args = {xt_ptp(pressure), xt_ptp(density), xt_ptp(velocityX), xt_ptp(velocityY), xt_ptp(velocityZ)};

        double optimal = 0;
        auto result = generic_LpDerivation(lp_type, n_row, n_col, e_f, coeffs, bounds, range_args, optimal);
        GE_Var_ebList eb_list(result[0], result[1], result[2], result[3], result[4]);
        return eb_list;
    }

    // naive generic error bound derivation method
    std::vector<xt::xarray<double>> generic_naiveDerivation(int n_args, const xt::xarray<Type>& max_partial, const xt::xarray<Type>& range_args, Type range_f, double e_f) const {
        double tau = range_f * e_f;
        xt::xarray<double> eb_uniform_rel_err = xt::ones<Type>({n_args}) * (tau / xt::sum(max_partial*range_args)());
        xt::xarray<double> eb_uniform_abs_err = xt::ones<Type>({n_args}) * (tau / xt::sum(max_partial)()) / range_args;
        xt::xarray<double> eb_uniform_component = xt::ones<Type>({n_args}) * (tau/n_args) / max_partial / range_args;
        return {eb_uniform_rel_err, eb_uniform_abs_err, eb_uniform_component};
    }

    // naive top-2-bottom
    std::vector<xt::xarray<double>> naive_VTOT_2_vx_vy_vz(double e_f) const {
        xt::xarray<Type> max_partial = {max_dV_dvelocityX(), max_dV_dvelocityY(), max_dV_dvelocityZ()};
        xt::xarray<Type> range_args = {xt_ptp(velocityX), xt_ptp(velocityY), xt_ptp(velocityZ)};
        Type range_f = xt_ptp(VTOT);
        auto res_V_xyz =generic_naiveDerivation(range_args.size(), max_partial, range_args, range_f, e_f);
        return res_V_xyz;
    }
    std::vector<xt::xarray<double>> naive_T_2_pressure_density(double e_f) const {
        xt::xarray<Type> max_partial = {max_dT_dpressure(), max_dT_ddensity()};
        xt::xarray<Type> range_args = {xt_ptp(pressure), xt_ptp(density)};
        Type range_f = xt_ptp(T);
        auto res_T_pd = generic_naiveDerivation(range_args.size(), max_partial, range_args, range_f, e_f);
        return res_T_pd;
    }
    std::vector<xt::xarray<double>> naive_C_2_pressure_density(double e_f) const {
        xt::xarray<Type> max_partial = {max_dC_dpressure(), max_dC_ddensity()};
        xt::xarray<Type> range_args = {xt_ptp(pressure), xt_ptp(density)};
        Type range_f = xt_ptp(C);
        auto res_C_pd = generic_naiveDerivation(range_args.size(), max_partial, range_args, range_f, e_f);
        return res_C_pd;
    }
    std::vector<xt::xarray<double>> naive_mu_2_pressure_density(double e_f) const {
        xt::xarray<Type> max_partial = {max_dmu_dpressure(), max_dmu_ddensity()};
        xt::xarray<Type> range_args = {xt_ptp(pressure), xt_ptp(density)};
        Type range_f = xt_ptp(mu);
        auto res_mu_pd = generic_naiveDerivation(range_args.size(), max_partial, range_args, range_f, e_f);
        return res_mu_pd;
    }
    std::vector<xt::xarray<double>> naive_M_2_all(double e_f) const {
        xt::xarray<Type> max_partial = {max_dM_dpressure(), max_dM_ddensity(), max_dM_dvelocityX(), max_dM_dvelocityY(), max_dM_dvelocityZ()};
        xt::xarray<Type> range_args = {xt_ptp(pressure), xt_ptp(density), xt_ptp(velocityX), xt_ptp(velocityY), xt_ptp(velocityZ)};
        Type range_f = xt_ptp(M);
        auto res_M_all = generic_naiveDerivation(range_args.size(), max_partial, range_args, range_f, e_f);
        return res_M_all;
    }
    std::vector<xt::xarray<double>> naive_PT_2_all(double e_f) const {
        xt::xarray<Type> max_partial = {max_dPT_dpressure(), max_dPT_ddensity(), max_dPT_dvelocityX(), max_dPT_dvelocityY(), max_dPT_dvelocityZ()};
        xt::xarray<Type> range_args = {xt_ptp(pressure), xt_ptp(density), xt_ptp(velocityX), xt_ptp(velocityY), xt_ptp(velocityZ)};
        Type range_f = xt_ptp(PT);
        auto res_PT_all = generic_naiveDerivation(range_args.size(), max_partial, range_args, range_f, e_f);
        return res_PT_all;
    }
    GE_Var_ebList naive_all_2_all(double e_f, int opt) const {
        xt::xarray<double> eb_xyz_4_V = naive_VTOT_2_vx_vy_vz(e_f)[opt];
        xt::xarray<double> eb_pd_4_T = naive_T_2_pressure_density(e_f)[opt];
        xt::xarray<double> eb_pd_4_C = naive_C_2_pressure_density(e_f)[opt];
        xt::xarray<double> eb_pd_4_mu = naive_mu_2_pressure_density(e_f)[opt];
        xt::xarray<double> eb_all_4_M = naive_M_2_all(e_f)[opt];
        xt::xarray<double> eb_all_4_PT = naive_PT_2_all(e_f)[opt];

        double eb_pres = xt::amin(xt::xarray<double>({eb_pd_4_T(0), eb_pd_4_C(0), eb_pd_4_mu(0), eb_all_4_M(0), eb_all_4_PT(0)}))();
        double eb_den = xt::amin(xt::xarray<double>({eb_pd_4_T(1), eb_pd_4_C(1), eb_pd_4_mu(1), eb_all_4_M(1), eb_all_4_PT(1)}))();
        double eb_vx = xt::amin(xt::xarray<double>({eb_xyz_4_V(0), eb_all_4_M(2), eb_all_4_PT(2)}))();
        double eb_vy = xt::amin(xt::xarray<double>({eb_xyz_4_V(1), eb_all_4_M(3), eb_all_4_PT(3)}))();
        double eb_vz = xt::amin(xt::xarray<double>({eb_xyz_4_V(2), eb_all_4_M(4), eb_all_4_PT(4)}))();

        GE_Var_ebList eb_list(eb_pres, eb_den, eb_vx, eb_vy, eb_vz);
        return eb_list;
    }

private:
    xt::xarray<Type> compute_VTOT(const xt::xarray<Type>& velocityX, const xt::xarray<Type>& velocityY, const xt::xarray<Type>& velocityZ) const {
        return xt::sqrt(xt::square(velocityX) + xt::square(velocityY) + xt::square(velocityZ));
    }
    xt::xarray<Type> compute_T(const xt::xarray<Type>& pressure, const xt::xarray<Type>& density) const {
        return pressure / density / R;
    }
    xt::xarray<Type> compute_C(const xt::xarray<Type>& pressure, const xt::xarray<Type>& density) const {
        return xt::sqrt(gamma * pressure/density);
    }
    xt::xarray<Type> compute_M(const xt::xarray<Type>& VTOT, const xt::xarray<Type>& C) const {
        return VTOT / C;
    }
    xt::xarray<Type> compute_PT(const xt::xarray<Type>& pressure, const xt::xarray<Type>& M) const {
        return pressure * xt::pow((1 + (gamma-1)/2 * xt::square(M)), gamma/(gamma-1));
    }
    xt::xarray<Type> compute_mu(const xt::xarray<Type>& pressure, const xt::xarray<Type>& density) const {
        xt::xarray<Type> T = compute_T(pressure, density);
        return mu_r * xt::pow(T/Tr, 1.5) * (Tr+S) / (T+S);
    }
};

}
#endif