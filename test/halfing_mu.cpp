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
bool halfing_error_mu_uniform(const T * P, const T * D, size_t n, const std::vector<unsigned char>& mask, const double tau, std::vector<double>& ebs){
	double eb_P = ebs[0];
	double eb_D = ebs[1];
	double R = 287.1;
	double gamma = 1.4;
	double mi = 3.5;
	double mu_r = 1.716e-5;
	double T_r = 273.15;
	double S = 110.4;
	double c_1 = 1.0 / R;
	double c_2 = sqrt(gamma * R);
	double c_3 = T_r + S;
	double max_value = 0;
	int max_index = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		double e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
		double Temp = P[i] / (D[i] * R);
		double e_TrS_TS = c_3 * compute_bound_radical(Temp, S, e_T);
		double TrS_TS = c_3 / (Temp + S);
		double e_T_Tr_3 = 3*pow(Temp/T_r, 2)*(e_T/T_r) + 3*Temp/T_r*(e_T/T_r)*(e_T/T_r) + (e_T/T_r)*(e_T/T_r)*(e_T/T_r);
		double T_Tr_3 = pow(Temp/T_r, 3);
		double e_T_Tr_3_sqrt = compute_bound_square_root_x(T_Tr_3, e_T_Tr_3);
		double T_Tr_3_sqrt = sqrt(T_Tr_3);
		double e_mu = mu_r * compute_bound_multiplication(T_Tr_3_sqrt, TrS_TS, e_T_Tr_3_sqrt, e_TrS_TS);
		double mu = mu_r * T_Tr_3_sqrt * TrS_TS;

		error_est_mu[i] = e_mu;
		error_mu[i] = mu - mu_ori[i];

		error_est_mu[i] = e_mu;
		error_mu[i] = mu - mu_ori[i];

		if(max_value < error_est_mu[i]){
			max_value = error_est_mu[i];
			max_index = i;
		}
	}
	std::cout << names[3] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		double estimate_error = max_value;
		double eb_P = ebs[0];
		double eb_D = ebs[1];
		while(estimate_error > tau){
    		std::cout << "uniform decrease\n";
			eb_P = eb_P / 1.5;
			eb_D = eb_D / 1.5;
			double e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
			double Temp = P[i] / (D[i] * R);
			double e_TrS_TS = c_3 * compute_bound_radical(Temp, S, e_T);
			double TrS_TS = c_3 / (Temp + S);
			double e_T_Tr_3 = 3*pow(Temp/T_r, 2)*(e_T/T_r) + 3*Temp/T_r*(e_T/T_r)*(e_T/T_r) + (e_T/T_r)*(e_T/T_r)*(e_T/T_r);
			double T_Tr_3 = pow(Temp/T_r, 3);
			double e_T_Tr_3_sqrt = compute_bound_square_root_x(T_Tr_3, e_T_Tr_3);
			double T_Tr_3_sqrt = sqrt(T_Tr_3);
			estimate_error = mu_r * compute_bound_multiplication(T_Tr_3_sqrt, TrS_TS, e_T_Tr_3_sqrt, e_TrS_TS);			
		}
		ebs[0] = eb_P;
		ebs[1] = eb_D;
		return false;
	}
	return true;
}


int main(int argc, char ** argv){

    using T = double;
    double target_rel_eb = atof(argv[1]);

    size_t num_elements = 0;
    P_ori = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    D_ori = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements);
    Vx_ori = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    Vy_ori = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    Vz_ori = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<double> ebs;
    ebs.push_back(compute_value_range(P_ori));
    ebs.push_back(compute_value_range(D_ori));
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
	std::vector<size_t> total_retrieved_size(n_variable, 0);

    int iter = 0;
    int max_iter = 5;
	bool tolerance_met = false;
    while((!tolerance_met) && (iter < max_iter)){
    	iter ++;
	    for(int i=0; i<n_variable; i++){
	        auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
			total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
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
	    P_dec = reconstructed_vars[0].data();
	    D_dec = reconstructed_vars[1].data();
	    MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
	    MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
	    error_mu = std::vector<double>(num_elements);
	    error_est_mu = std::vector<double>(num_elements);
		std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
	    MDR::print_vec(ebs);
	    tolerance_met = halfing_error_mu_uniform(P_dec, D_dec, num_elements, mask, tau[5], ebs);
		std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
	    MDR::print_vec(ebs);
	    std::cout << names[5] << " requested error = " << tau[5] << std::endl;
	    print_max_abs(names[5] + " error", error_mu);
	    print_max_abs(names[5] + " error_est", error_est_mu);   	
    }
	std::cout << "iter = " << iter << std::endl;
   
   	size_t total_size = std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0);
	double cr = n_variable * num_elements * sizeof(T) * 1.0 / total_size;
	std::cout << "each retrieved size:" << std::endl;
	MDR::print_vec(total_retrieved_size);
	std::cout << "aggregated cr = " << cr << std::endl;

    return 0;
}