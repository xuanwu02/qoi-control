#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <chrono>
#include <dirent.h>
#include <zstd.h>
#include "adios2.h"
#include "mpi.h"
#include "Reconstructor/Reconstructor.hpp"
#include "Refactor/Refactor.hpp"
#include "SZ3/api/sz.hpp"

const std::vector<std::string> var_name{"U_aver", "V_aver", "W_aver", "Pressure", "Rho"};
const std::vector<std::string> var_name_out{"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};
const int n_vars = 5;

using namespace MDR; 

template <class T>
void print_info(int p, const std::string& name, const std::vector<T>& vec){
    T max = vec[0];
    T min = vec[0];
    for(int i=1; i<vec.size(); i++){
        if(max < vec[i]) max = vec[i];
        if(min > vec[i]) min = vec[i];
    }
    printf("Processor %d var %s (size = %d): min = %.4f, max = %.4f\n", p, name.c_str(), vec.size(), min, max);
}

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class Writer>
MDR::ComposedRefactor<T, Decomposer, Interleaver, Encoder, Compressor, ErrorCollector, Writer> generateRefactor(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer){
    auto refactor = MDR::ComposedRefactor<T, Decomposer, Interleaver, Encoder, Compressor, ErrorCollector, Writer>(decomposer, interleaver, encoder, compressor, collector, writer);
    return refactor;
}

template<class Type>
void refactor_singleZone(int rank, int num_elements, std::vector<Type>& velocityX_vec, std::vector<Type>& velocityY_vec, std::vector<Type>& velocityZ_vec, std::vector<Type>& pressure_vec, std::vector<Type>& density_vec){
    std::string filename = "/pscratch/xli281_uksr/xliang/GE/sol_4114800_aver_b" + std::to_string(rank) + ".bp/";
    MGARD::writefile((filename + var_name_out[0] + ".dat").c_str(), velocityX_vec.data(), velocityX_vec.size());
    MGARD::writefile((filename + var_name_out[1] + ".dat").c_str(), velocityY_vec.data(), velocityX_vec.size());
    MGARD::writefile((filename + var_name_out[2] + ".dat").c_str(), velocityZ_vec.data(), velocityX_vec.size());
    MGARD::writefile((filename + var_name_out[3] + ".dat").c_str(), pressure_vec.data(), velocityX_vec.size());
    MGARD::writefile((filename + var_name_out[4] + ".dat").c_str(), density_vec.data(), velocityX_vec.size());
    std::vector<uint32_t> dims;
    dims.push_back(num_elements);
    // compute masks
    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){
            mask[i] = 1;
            num_valid_data ++;
        }
    }
    std::string mask_file = "/pscratch/xli281_uksr/xliang/GE/block_" + std::to_string(rank) + "_refactored/mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    std::vector<std::vector<Type>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<uint32_t> dims_masked;
    dims_masked.push_back(num_valid_data);
    std::vector<Type> buffer(num_valid_data);

    bool use_negabinary = false;
    using Decomposer = MGARDHierarchicalDecomposer<Type>;
    using Interleaver = DirectInterleaver<Type>;
    using Encoder = PerBitBPEncoder<Type, uint32_t>;
    using Compressor = AdaptiveLevelCompressor;
    using ErrorCollector = SquaredErrorCollector<Type>;
    using Writer = ConcatLevelFileWriter;
    const int target_level = 8;
    const int num_bitplanes = 60;

    // using Decomposer = MGARDOrthoganalDecomposer<Type>;
    // using Encoder = NegaBinaryBPEncoder<Type, uint32_t>;
    // use_negabinary = true;

    for(int i=0; i<n_vars; i++){
        std::string rdir_prefix = "/pscratch/xli281_uksr/xliang/GE/block_" + std::to_string(rank) + "_refactored/" + var_name_out[i] + "/";
        std::string metadata_file = rdir_prefix + "metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
    	// std::cout << metadata_file << std::endl;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "level_" + std::to_string(i) + ".bin";
	        // std::cout << filename << "\n";
            files.push_back(filename);
        }
        auto decomposer = Decomposer();
        auto interleaver = Interleaver();
        auto encoder = Encoder();
        auto compressor = Compressor(64);
        auto collector = ErrorCollector();
        auto writer = Writer(metadata_file, files);
        auto refactor = generateRefactor<Type>(decomposer, interleaver, encoder, compressor, collector, writer);
        refactor.negabinary = use_negabinary;
        if(i < 3){
            // use masked refactoring for vx vy vz
            int index = 0;
            for(int j=0; j<num_elements; j++){
                if(mask[j]){
                    buffer[index ++] = vars_vec[i][j];
                }
            }
            refactor.refactor(buffer.data(), dims_masked, target_level, num_bitplanes);
        }
        else{
            refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes);
        }
    }
}

int main(int argc, char **argv) {
    MPI_Init(NULL, NULL);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    //if(size != 64){
    //    std::cout << "Need to use 64 processors for 64 blocks\n";
    //    MPI_Abort(MPI_COMM_WORLD, -1);
    //}

    const int target_level = 8;

    std::string filename = "/pscratch/xli281_uksr/xliang/GE/sol_4114800_aver_b" + std::to_string(rank) + ".bp";
    adios2::ADIOS ad; 
    adios2::IO reader_io = ad.DeclareIO("Input");
    adios2::Engine reader = reader_io.Open(filename, adios2::Mode::Read);

    while (true) {
        // Begin step
        adios2::StepStatus read_status = reader.BeginStep(adios2::StepMode::Read, 10.0f);
        if (read_status == adios2::StepStatus::NotReady) {
            // std::cout << "Stream not ready yet. Waiting...\n";
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            continue;
        }
        else if (read_status != adios2::StepStatus::OK) {
            break;
        }

        adios2::Variable<double> var_ad2;
        std::vector<std::vector<double>> vars_vec(n_vars, std::vector<double>());

        for (int i=0; i<n_vars; i++) {
            var_ad2 = reader_io.InquireVariable<double>(var_name[i]);
            std::vector<std::size_t> shape = var_ad2.Shape();
            reader.Get(var_ad2, vars_vec[i], adios2::Mode::Sync);
            reader.PerformGets();
            //print_info(rank, var_name[i], vars_vec[i]);
        }
	    refactor_singleZone(rank, vars_vec[0].size(), vars_vec[0], vars_vec[1], vars_vec[2], vars_vec[3], vars_vec[4]); 
        reader.EndStep();
    }
    reader.Close();

    MPI_Finalize();
    
    return 0;
}
