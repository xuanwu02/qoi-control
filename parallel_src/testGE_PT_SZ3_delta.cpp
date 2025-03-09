#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include <mpi.h>
#include "utils.hpp"
#include "SZ3/api/sz.hpp"

const std::vector<std::string> var_name_out{"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};
const int n_vars = 5;

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

template <class T>
void compute_QoIs(const T * Vx, const T * Vy, const T * Vz, const T * P, const T * D, size_t n,
					T * V_TOT_, T * Temp_, T * C_, T * Mach_, T * PT_, T * mu_){
	double R = 287.1;
	double gamma = 1.4;
	double mi = 3.5;
	double mu_r = 1.716e-5;
	double T_r = 273.15;
	double S = 110.4;
	for(int i=0; i<n; i++){
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		double V_TOT = sqrt(V_TOT_2);
		double Temp = P[i] / (D[i] * R);
		double C = sqrt(gamma * R * Temp);
		double Mach = V_TOT / C;
		double Mach_tmp = 1 + (gamma-1)/2 * Mach * Mach;
		double Mach_tmp_mi = sqrt(pow(Mach_tmp, 7));
		double PT = P[i] * Mach_tmp_mi;
		double TrS_TS = (T_r + S) / (Temp + S);
		double T_Tr_3 = pow(Temp/T_r, 3);
		double T_Tr_3_sqrt = sqrt(T_Tr_3);
		double mu = mu_r * T_Tr_3_sqrt * TrS_TS;
		V_TOT_[i] = V_TOT;
		Temp_[i] = Temp;
		C_[i] = C;
		Mach_[i] = Mach;
		PT_[i] = PT;
		mu_[i] = mu;
	}
}

template <class T>
T compute_value_range(const std::vector<T>& vec){
	T min = vec[0];
	T max = vec[0];
	for(int i=0; i<vec.size(); i++){
		if(vec[i] < min) min = vec[i];
		if(vec[i] > max) max = vec[i];
	}
	return max - min;
}

// f(x) = x^2
template <class T>
inline double compute_bound_x_square(T x, T eb){
	return 2*fabs(x)*eb + eb*eb;
}

// f(x) = sqrt(x)
template <class T>
inline double compute_bound_square_root_x(T x, T eb){
	if(x == 0) {
		return sqrt(eb);
	}
	if(x > eb){
		return eb / (sqrt(x - eb) + sqrt(x));
	}
	else{
		return eb / sqrt(x);
	}
}

template <class T>
inline double compute_bound_multiplication(T x, T y, T eb_x, T eb_y){
	return fabs(x)*eb_y + fabs(y)*eb_x + eb_x*eb_y;
}

template <class T>
inline double compute_bound_division(T x, T y, T eb_x, T eb_y){
	if(eb_y < fabs(y)){
		double e = fabs(x)*eb_y + fabs(y)*eb_x;
		if(y > 0) return e / (y*(y - eb_y));
		else return e/ (y*(y + eb_y));
	}
	else{
		std::cout << "Warning: cannot control error in x/y\n";
		return 0;
	}
}

template<class T>
bool halfing_error_PT_uniform(const T * Vx, const T * Vy, const T * Vz, const T * P, const T * D, size_t n, const std::vector<unsigned char>& mask, const double tau, std::vector<double>& ebs){
	double eb_Vx = ebs[0];
	double eb_Vy = ebs[1];
	double eb_Vz = ebs[2];
	double eb_P = ebs[3];
	double eb_D = ebs[4];
	double R = 287.1;
	double gamma = 1.4;
	double mi = 3.5;
	double mu_r = 1.716e-5;
	double T_r = 273.15;
	double S = 110.4;
	double c_1 = 1.0 / R;
	double c_2 = sqrt(gamma * R);
	int C7i[8] = {1, 7, 21, 35, 35, 21, 7, 1};
	double max_value = 0;
	int max_index = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		double e_V_TOT_2 = 0;
		if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		double e_V_TOT = 0;
		if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		double V_TOT = sqrt(V_TOT_2);
		double e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
		double Temp = P[i] / (D[i] * R);
		double e_C = c_2*compute_bound_square_root_x(Temp, e_T);
		double C = c_2 * sqrt(Temp);
		double e_Mach = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
		double Mach = V_TOT / C;
		double e_Mach_tmp = (gamma-1) / 2 * compute_bound_x_square(Mach, e_Mach);
		double Mach_tmp = 1 + (gamma-1)/2 * Mach * Mach;
		double e_Mach_tmp_mi = 0;
		for(int i=1; i<=7; i++){
			e_Mach_tmp_mi += C7i[i] * pow(Mach_tmp, 7-i) * pow(e_Mach_tmp, i);
		}
		double Mach_tmp_mi = sqrt(pow(Mach_tmp, 7));
		double e_PT = compute_bound_multiplication(P[i], Mach_tmp_mi, eb_P, e_Mach_tmp_mi);
		double PT = P[i] * Mach_tmp_mi;

		error_est_PT[i] = e_PT;
		error_PT[i] = PT - PT_ori[i];
		if(max_value < error_est_PT[i]){
			max_value = error_est_PT[i];
			max_index = i;
		}
	}
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		double estimate_error = max_value;
		double eb_Vx = ebs[0];
		double eb_Vy = ebs[1];
		double eb_Vz = ebs[2];
		double eb_P = ebs[3];
		double eb_D = ebs[4];
		while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
			eb_Vx = eb_Vx / 1.5;
			eb_Vy = eb_Vy / 1.5;
			eb_Vz = eb_Vz / 1.5; 
			eb_P = eb_P / 1.5;
			eb_D = eb_D / 1.5;
			double e_V_TOT_2 = 0;
			if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
			double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
			double e_V_TOT = 0;
			if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			double V_TOT = sqrt(V_TOT_2);
			double e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
			double Temp = P[i] / (D[i] * R);
			double e_C = c_2*compute_bound_square_root_x(Temp, e_T);
			double C = c_2 * sqrt(Temp);
			double e_Mach = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
			double Mach = V_TOT / C;
			double e_Mach_tmp = (gamma-1) / 2 * compute_bound_x_square(Mach, e_Mach);
			double Mach_tmp = 1 + (gamma-1)/2 * Mach * Mach;
			double e_Mach_tmp_mi = 0;
			for(int i=1; i<=7; i++){
				e_Mach_tmp_mi += C7i[i] * pow(Mach_tmp, 7-i) * pow(e_Mach_tmp, i);
			}
			double Mach_tmp_mi = sqrt(pow(Mach_tmp, 7));
			estimate_error = compute_bound_multiplication(P[i], Mach_tmp_mi, eb_P, e_Mach_tmp_mi);
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		ebs[3] = eb_P;
		ebs[4] = eb_D;
		return false;
	}
	return true;
}

template <class T>
T print_max_abs(int rank, const std::string& name, const std::vector<T>& vec){
	T max = fabs(vec[0]);
	for(int i=1; i<vec.size(); i++){
		if(max < fabs(vec[i])) max = fabs(vec[i]);
	}
	// printf("Processor %d var %s: max absolute value =  %.4f\n", rank, name.c_str(), max);
	return max;
}

template<class T>
char * SZ3_compress(size_t num_elements, T * data, double abs_eb, size_t& compressed_size){
    SZ3::Config conf(num_elements);
    conf.cmprAlgo = SZ3::ALGO_INTERP_LORENZO;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = abs_eb;
    size_t cmpSize = 0;
    char *cmpData = SZ_compress<T>(conf, data, cmpSize);
    compressed_size = cmpSize;
    return cmpData;
}

template<class T>
void SZ3_decompress(char * cmpData, size_t compressed_size, T * dec_data){
    SZ3::Config conf1;
    SZ_decompress<T>(conf1, cmpData, compressed_size, dec_data);
}

inline int find_index(double target_rel_eb, double& rel_eb){
    int i = 0;
    while(target_rel_eb < rel_eb){
        i ++;
        rel_eb /= 10;
    }
    return i;
}

int main(int argc, char ** argv){

	MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    using T = double;
    const int target_level = 8;
    // read_file
    size_t num_elements = 0;
	std::string filename = "/pscratch/xli281_uksr/xliang/GE/sol_4114800_aver_b" + std::to_string(rank) + ".bp/"; 
	Vx_ori = MGARD::readfile<T>((filename + var_name_out[0] + ".dat").c_str(), num_elements);
    Vy_ori = MGARD::readfile<T>((filename + var_name_out[1] + ".dat").c_str(), num_elements);
    Vz_ori = MGARD::readfile<T>((filename + var_name_out[2] + ".dat").c_str(), num_elements);
    P_ori = MGARD::readfile<T>((filename + var_name_out[3] + ".dat").c_str(), num_elements);
    D_ori = MGARD::readfile<T>((filename + var_name_out[4] + ".dat").c_str(), num_elements);

    double target_rel_eb = atof(argv[1]);
    std::vector<double> ebs;
    ebs.push_back(compute_value_range(Vx_ori));
    ebs.push_back(compute_value_range(Vy_ori));
    ebs.push_back(compute_value_range(Vz_ori));
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

    double tau;
	double local_max = -9999, local_min = 9999;
	double global_max = 0, global_min = 0;
	for(int i=0; i<num_elements; i++){
		if(PT[i] > local_max) local_max = PT[i];
		if(PT[i] < local_min) local_min = PT[i];
	}
	MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	double max_value_range = global_max - global_min;
	tau = max_value_range*target_rel_eb;
	PT_ori = PT.data();

    std::vector<double> value_range(n_vars);
    //for(int i=0; i<n_vars; i++){
    //    value_range[i] = compute_value_range(vars_vec[i]);
    //}
    value_range[0] = compute_value_range(Vx_ori);
    value_range[1] = compute_value_range(Vy_ori);
    value_range[2] = compute_value_range(Vz_ori);
    value_range[3] = compute_value_range(P_ori);
    value_range[4] = compute_value_range(D_ori);
    std::string mask_file = "/pscratch/xli281_uksr/xliang/GE/block_" + std::to_string(rank) + "_refactored/mask.bin";
    size_t num_valid_data = 0;
    auto mask = MGARD::readfile<unsigned char>(mask_file.c_str(), num_valid_data);

    int iter = 0;
    int max_iter = 10;
	bool tolerance_met = false;
    size_t total_size = 0;
	double local_elapsed_time = 0, max_time = 0;
    std::vector<size_t> total_retrieved_sizes(n_vars, 0);
    std::vector<int> current_ind(n_vars, -1);
    std::vector<int> prev_ind(n_vars, 0);
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<double>(num_elements));
	T * reconstructed_data = (T *) malloc(num_elements * sizeof(T));
	local_elapsed_time = -MPI_Wtime();
    while((!tolerance_met) && (iter < max_iter)){
    	iter ++;
	    for(int i=0; i<n_variable; i++){
	        std::string rdir_prefix = "/pscratch/xli281_uksr/xliang/GE/block_" + std::to_string(rank) + "_refactored/" + var_name_out[i] + "/";
            double rel_eb = ebs[i]/value_range[i];
            double file_eb = 0.1;
            auto file_ind = find_index(rel_eb, file_eb);
	    	if(file_ind > 9) file_ind = 9;
			if(file_ind > current_ind[i]){
                for(int j=current_ind[i]+1; j<=file_ind; j++){
	                std::string filename = rdir_prefix + "SZ3_delta_eb_" + std::to_string(j) + ".bin";
                    size_t n = 0;
                    auto cmpData = MGARD::readfile<char>(filename.c_str(), n);
                    total_retrieved_sizes[i] += n;
		    		total_size += n;
                    SZ3_decompress(cmpData.data(), n, reconstructed_data);
                    if(i < 3){
                        // reconstruct with mask
                        int index = 0;
                        for(int j=0; j<num_elements; j++){
                            if(mask[j]){
                                reconstructed_vars[i][j] += reconstructed_data[index ++];
                            }
                            else reconstructed_vars[i][j] = 0;
                        }
                    }
                    else{
                        for(int j=0; j<num_elements; j++){
                            reconstructed_vars[i][j] += reconstructed_data[j];
                        }
                    }
                }
                current_ind[i] = file_ind;
            }
	    }
	    Vx_dec = reconstructed_vars[0].data();
	    Vy_dec = reconstructed_vars[1].data();
	    Vz_dec = reconstructed_vars[2].data();
	    P_dec = reconstructed_vars[3].data();
	    D_dec = reconstructed_vars[4].data();
	    // MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
	    // MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
	    // MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
	    // MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
	    // MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
	    error_PT = std::vector<double>(num_elements);
	    error_est_PT = std::vector<double>(num_elements);
	    tolerance_met = halfing_error_PT_uniform(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs);
    }
	int maxit;
	MPI_Reduce(&iter, &maxit, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) std::cout << "max iter = " << maxit << std::endl;
	// std::cout << "rank " << rank << " max iter = " << iter << std::endl;
    free(reconstructed_data);
	local_elapsed_time += MPI_Wtime();
	MPI_Reduce(&local_elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    // printf("Processor %d total size = %d\n", rank, total_size);
    // MPI_Barrier(MPI_COMM_WORLD);
    if(!rank) printf("Target PT error = %.4f\n", tau);
	double max_error = 0;
	max_error = print_max_abs(rank, "PT error", error_PT);
	double max_pt_error = 0;
	MPI_Reduce(&max_error, &max_pt_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("Max aggregated PT error = %.4f\n", max_pt_error);
	max_error = print_max_abs(rank, "PT error", error_est_PT);
	double max_pt_error_est = 0;
	MPI_Reduce(&max_error, &max_pt_error_est, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("Max aggregated PT est error = %.4f\n", max_pt_error_est);
	unsigned long long int total_num = 0;
	MPI_Reduce(&num_elements, &total_num, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	unsigned long long int total_retrieved = 0;
	MPI_Reduce(&total_size, &total_retrieved, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	if(!rank) printf("Aggregated bitrate = %.4f, retrieved_size = %ld, total_num_elements = %ld\n", 8*total_retrieved * 1.0 / (total_num * n_variable), total_retrieved, total_num);
	if(!rank) printf("elapsed_time = %.6f\n", max_time);
	std::vector<char> data_buffer(num_elements*sizeof(double));
	std::vector<unsigned long long int> offsets;
	std::vector<size_t> counts;
    for(int i=0; i<n_variable; i++){
		std::string rdir_prefix = "/pscratch/xli281_uksr/xliang/GE/block_" + std::to_string(rank) + "_refactored/" + var_name_out[i] + "/";
    	// copy local data
    	size_t local_size = 0;
    	char * data_buffer_pos = data_buffer.data();
		std::vector<size_t> sizes(current_ind[i]+1, 0);
    	for(int j=prev_ind[i]; j<=current_ind[i]; j++){
			size_t n = 0;
			std::string filename = rdir_prefix + "SZ3_delta_eb_" + std::to_string(j) + ".bin";
			auto cmpData = MGARD::readfile<char>(filename.c_str(), n); 
			memcpy(data_buffer_pos, cmpData.data(), n);	
			data_buffer_pos += n;
			local_size += n;
			sizes[j] = n;
    	}
		std::string sizes_file = "/project/xli281_uksr/xwu/test_GE/sizes/psz3delta_" + std::to_string(rank) + "_" + std::to_string(i) + "_" + std::string(argv[1]);
		MGARD::writefile(sizes_file.c_str(), sizes.data(), sizes.size());

    	unsigned long long int count = local_size;
		unsigned long long int offset = 0;
		unsigned long long int buffer = 0;
		for(int j=0; j<size; j++){
			if(j == rank){
				if(j != 0) {
					MPI_Recv(&offset, 1, MPI_UNSIGNED_LONG_LONG, j-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				buffer = offset + count;
				if(j != size - 1) MPI_Send(&buffer, 1, MPI_UNSIGNED_LONG_LONG, j+1, 0, MPI_COMM_WORLD);
			}
		}
		offsets.push_back(offset);
		counts.push_back(local_size);
		MPI_File file;
		std::string filename = "/pscratch/xli281_uksr/xliang/GE/PSZ3delta/PT_PSZ3delta_" + std::to_string(target_rel_eb) + "_" + var_name_out[i] + "_aggregated.dat";
		MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
		MPI_File_write_at(file, offset, data_buffer.data(), local_size, MPI_SIGNED_CHAR, MPI_STATUS_IGNORE);
		MPI_File_close(&file);
    }
	// std::string offsets_file = "/project/xli281_uksr/xwu/test_GE/offsets/psz3delta_" + std::to_string(rank) + "_" + argv[1];
	// MGARD::writefile(offsets_file.c_str(), offsets.data(), offsets.size());
	// std::string count_file = "/project/xli281_uksr/xwu/test_GE/count/psz3delta_" + std::to_string(rank) + "_" + argv[1];
	// MGARD::writefile(count_file.c_str(), counts.data(), counts.size());

    MPI_Finalize();
    return 0;
}
