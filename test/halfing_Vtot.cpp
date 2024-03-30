#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include "utils.hpp"
#include "qoi_utils.hpp"
#include "Reconstructor/Reconstructor.hpp"
#include "Synthesizer4GE.hpp"

using namespace MDR;

std::vector<double> P_ori;
std::vector<double> D_ori;
std::vector<double> Vx_ori;
std::vector<double> Vy_ori;
std::vector<double> Vz_ori;
double * P_dec = NULL;
double * D_dec = NULL;
double * Vx_dec = NULL;
double * Vy_dec = NULL;
double * Vz_dec = NULL;
double * V_TOT_ori = NULL;
double * Temp_ori = NULL;
double * C_ori = NULL;
double * Mach_ori = NULL;
double * PT_ori = NULL;
double * mu_ori = NULL;
std::vector<double> error_V_TOT;
std::vector<double> error_Temp;
std::vector<double> error_C;
std::vector<double> error_Mach;
std::vector<double> error_PT;
std::vector<double> error_mu;
std::vector<double> error_est_V_TOT;
std::vector<double> error_est_Temp;
std::vector<double> error_est_C;
std::vector<double> error_est_Mach;
std::vector<double> error_est_PT;
std::vector<double> error_est_mu;


template<class T>
bool halfing_error_V_TOT(const T * Vx, const T * Vy, const T * Vz, size_t n, const std::vector<unsigned char>& mask, const double tau, std::vector<double>& ebs){
	double eb_Vx = ebs[0];
	double eb_Vy = ebs[1];
	double eb_Vz = ebs[2];
	double max_value = 0;
	int max_index = 0;
	for(int i=0; i<n; i++){
		// error of total velocity square
		double e_V_TOT_2 = 0;
		if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		// error of total velocity
		double e_V_TOT = 0;
		if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		double V_TOT = sqrt(V_TOT_2);
		// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);

		error_est_V_TOT[i] = e_V_TOT;
		error_V_TOT[i] = V_TOT - V_TOT_ori[i];

		if(max_value < error_est_V_TOT[i]){
			max_value = error_est_V_TOT[i];
			max_index = i;
		}

	}
	std::cout << "Vx = " << Vx[max_index] << " Vy = " << Vy[max_index] << " Vz = " << Vz[max_index] << std::endl;
	std::cout << "eb_Vx = " << eb_Vx << " eb_Vy = " << eb_Vy << " eb_Vz = " << eb_Vz << std::endl;
	std::cout << "coeff_x = " << fabs(Vx[max_index])*eb_Vx << " coeff_y = " << fabs(Vy[max_index])*eb_Vy <<  " coeff_z = " << fabs(Vz[max_index])*eb_Vz << std::endl;
	std::cout << names[0] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		std::cout << "comparing: " << max_value << " " << tau << std::endl;
		std::vector<double> new_ebs(3, 0);
		// error of Vtotal
        {
			new_ebs[0] = ebs[0] / 2;
			new_ebs[1] = ebs[1] / 2;
			new_ebs[2] = ebs[2] / 2;            
			// double coeff_0 = fabs(Vy[max_index])*eb_Vy * fabs(Vz[max_index])*eb_Vz;
			// double coeff_1 = fabs(Vx[max_index])*eb_Vx * fabs(Vz[max_index])*eb_Vz;
			// double coeff_2 = fabs(Vx[max_index])*eb_Vx * fabs(Vy[max_index])*eb_Vy;
			// double coeff_0 = fabs(Vx[max_index])*eb_Vx;
			// double coeff_1 = fabs(Vy[max_index])*eb_Vy;
			// double coeff_2 = fabs(Vz[max_index])*eb_Vz;
			// double average = (coeff_0 + coeff_1 + coeff_2)/3;
			// new_ebs[0] = ebs[0] / coeff_0 * average / 2;
			// new_ebs[1] = ebs[1] / coeff_1 * average / 2;
			// new_ebs[2] = ebs[2] / coeff_2 * average / 2;            
        }
		for(int i=0; i<3; i++){
			ebs[i] = new_ebs[i];
		}
		return false;
	}
	return true;
}



int main(int argc, char ** argv){

    using T = double;
    int zone_id = 0;
    std::string id_str = std::to_string(zone_id);
    size_t num_elements = 0;
    std::string prefix = data_file_prefix + "zone_";
    P_ori = MGARD::readfile<T>((prefix + id_str + "_Pressure.dat").c_str(), num_elements);
    D_ori = MGARD::readfile<T>((prefix + id_str + "_Density.dat").c_str(), num_elements);
    Vx_ori = MGARD::readfile<T>((prefix + id_str + "_VelocityX.dat").c_str(), num_elements);
    Vy_ori = MGARD::readfile<T>((prefix + id_str + "_VelocityY.dat").c_str(), num_elements);
    Vz_ori = MGARD::readfile<T>((prefix + id_str + "_VelocityZ.dat").c_str(), num_elements);
    double target_rel_eb = 0.01;
    std::vector<double> ebs;
    ebs.push_back(compute_value_range(Vx_ori));
    ebs.push_back(compute_value_range(Vy_ori));
    ebs.push_back(compute_value_range(Vz_ori));
	int n_variable = ebs.size();

    for(int i=0; i<ebs.size(); i++){
    	ebs[i] *= target_rel_eb;
    }

    std::vector<T> V_TOT(num_elements);
    std::vector<T> Temp(num_elements);
    std::vector<T> C(num_elements);
    std::vector<T> Mach(num_elements);
    std::vector<T> PT(num_elements);
    std::vector<T> mu(num_elements);
    compute_QoIs(Vx_ori.data(), Vy_ori.data(), Vz_ori.data(), P_ori.data(), D_ori.data(), num_elements, V_TOT.data(), Temp.data(), C.data(), Mach.data(), PT.data(), mu.data());

    std::vector<double> tau;
    tau.push_back(compute_value_range(V_TOT)*target_rel_eb);
    tau.push_back(compute_value_range(Temp)*target_rel_eb);
    tau.push_back(compute_value_range(C)*target_rel_eb);
    tau.push_back(compute_value_range(Mach)*target_rel_eb);
    tau.push_back(compute_value_range(PT)*target_rel_eb);
    tau.push_back(compute_value_range(mu)*target_rel_eb);
	V_TOT_ori = V_TOT.data();
	Temp_ori = Temp.data();
	C_ori = C.data();
	Mach_ori = Mach.data();
	PT_ori = PT.data();
	mu_ori = mu.data();

    std::string mask_file = rdata_file_prefix + "mask.bin";
    size_t num_valid_data = 0;
    auto mask = MGARD::readfile<unsigned char>(mask_file.c_str(), num_valid_data);
    std::vector<MDR::ComposedReconstructor<T, MGARDHierarchicalDecomposer<T>, DirectInterleaver<T>, PerBitBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
    for(int i=0; i<n_variable; i++){
        // std::string rdir_prefix = rdata_file_prefix + "zone_" + id_str + "_" + varlist[i];
        std::string rdir_prefix = rdata_file_prefix + varlist[i];
        std::string metadata_file = rdir_prefix + "_refactored_data/metadata.bin";
        std::vector<std::string> files;
        int num_levels = 9;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored_data/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto decomposer = MGARDHierarchicalDecomposer<T>();
        auto interleaver = DirectInterleaver<T>();
        auto encoder = PerBitBPEncoder<T, uint32_t>();
        auto compressor = AdaptiveLevelCompressor(64);
        auto estimator = MaxErrorEstimatorHB<T>();
        auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
        auto retriever = ConcatLevelFileRetriever(metadata_file, files);
        reconstructors.push_back(generateReconstructor<T>(decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever));
        reconstructors.back().load_metadata();
    }    
    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<double>(num_elements));

    int iter = 0;
    int max_iter = 5;
	bool tolerance_met = false;
    while((!tolerance_met) && (iter < max_iter)){
    	iter ++;
	    for(int i=0; i<n_variable; i++){
	        auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
	        if(i < 3){
	            // reconstruct with mask
	            int index = 0;
	            for(int j=0; j<num_elements; j++){
	                if(mask[j]){
	                    reconstructed_vars[i][j] = reconstructed_data[index ++];
	                }
	                else reconstructed_vars[i][j] = 0;
	            }
	        }
	        else{
	            memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
	        }
	    }

	    Vx_dec = reconstructed_vars[0].data();
	    Vy_dec = reconstructed_vars[1].data();
	    Vz_dec = reconstructed_vars[2].data();
	    MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
	    MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
	    MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
	    error_V_TOT = std::vector<double>(num_elements);
	    error_est_V_TOT = std::vector<double>(num_elements);
		std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
	    MDR::print_vec(ebs);
	    tolerance_met = halfing_error_V_TOT(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau[0], ebs);
		std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
	    MDR::print_vec(ebs);
	    std::cout << names[0] << " requested error = " << tau[0] << std::endl;
	    print_max_abs(names[0] + " error", error_V_TOT);
	    print_max_abs(names[0] + " error_est", error_est_V_TOT);   	
    }
	std::cout << "iter = " << iter << std::endl;
   
    return 0;
}