#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include "utils.hpp"
#include "QoIsWizard/QoIWizard4GE_nonLp.hpp"


int main(int argc, char** argv){
    using Type = double;

    int zone_id = 0;
    std::string id_str = std::to_string(zone_id);
    size_t num_elements = 0;
    // auto pressure_vec = MGARD::readfile<Type>(("/Users/xuanwu/github/datasets/GE/data/" + "zone_" + id_str + "_Pressure.dat").c_str(), num_elements);
    // auto density_vec = MGARD::readfile<Type>(("/Users/xuanwu/github/datasets/GE/data/" + "zone_" + id_str + "_Density.dat").c_str(), num_elements);
    // auto velocityX_vec = MGARD::readfile<Type>(("/Users/xuanwu/github/datasets/GE/data/" + "zone_" + id_str + "_VelocityX.dat").c_str(), num_elements);
    // auto velocityY_vec = MGARD::readfile<Type>(("/Users/xuanwu/github/datasets/GE/data/" + "zone_" + id_str + "_VelocityY.dat").c_str(), num_elements);
    // auto velocityZ_vec = MGARD::readfile<Type>(("/Users/xuanwu/github/datasets/GE/data/" + "zone_" + id_str + "_VelocityZ.dat").c_str(), num_elements);
    std::string prefix = "/Users/xliang/Dataset/GE/data/zone_";
    auto pressure_vec = MGARD::readfile<Type>((prefix + id_str + "_Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<Type>((prefix + id_str + "_Density.dat").c_str(), num_elements);
    auto velocityX_vec = MGARD::readfile<Type>((prefix + id_str + "_VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<Type>((prefix + id_str + "_VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<Type>((prefix + id_str + "_VelocityZ.dat").c_str(), num_elements);

    MDR::GEQoI_gen<Type> qoi_gen = MDR::GEQoI_gen<Type>(pressure_vec, density_vec, velocityX_vec, velocityY_vec, velocityZ_vec, num_elements);

    Type max_dT_dpres = qoi_gen.max_dT_dpressure();
    std::cout << "max_dT_dpres = " << max_dT_dpres << std::endl;
    return 0;
}