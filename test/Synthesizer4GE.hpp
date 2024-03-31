#ifndef _MDR_GE_SYNTHESIZER_HPP
#define _MDR_GE_SYNTHESIZER_HPP

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>
#include <cmath>
#include <numeric>
#include "Reconstructor/Reconstructor.hpp"
#include "Refactor/Refactor.hpp"
#include "SZ3/api/sz.hpp"

const std::string data_file_prefix = "/Users/xliang/Dataset/GE/data/";
const std::string rdata_file_prefix = "/Users/xliang/Dataset/GE/refactor/";
const std::vector<std::string> varlist = {"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};
const int n_vars = 5;

namespace MDR {

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class Writer>
MDR::ComposedRefactor<T, Decomposer, Interleaver, Encoder, Compressor, ErrorCollector, Writer> generateRefactor(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer){
    auto refactor = MDR::ComposedRefactor<T, Decomposer, Interleaver, Encoder, Compressor, ErrorCollector, Writer>(decomposer, interleaver, encoder, compressor, collector, writer);
    return refactor;
}

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever> generateReconstructor(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    return reconstructor;
}

const int target_level = 8;
const int num_bitplanes = 60;

template<class Type>
void refactor_singleZone(int zone_id){
    std::string id_str = std::to_string(zone_id);

    size_t num_elements = 0;
    auto pressure_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_Density.dat").c_str(), num_elements);
    auto velocityX_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityZ.dat").c_str(), num_elements);
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
    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    std::vector<std::vector<Type>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<uint32_t> dims_masked;
    dims_masked.push_back(num_valid_data);
    std::vector<Type> buffer(num_valid_data);

    bool use_negabinary = false;
    using Decomposer = MGARDHierarchicalDecomposer<Type>;
    using Interleaver = DirectInterleaver<Type>;
    // using Encoder = PerBitBPEncoder<Type, uint32_t>;
    using Compressor = AdaptiveLevelCompressor;
    using ErrorCollector = SquaredErrorCollector<Type>;
    using Writer = ConcatLevelFileWriter;

    // using Decomposer = MGARDOrthoganalDecomposer<Type>;
    using Encoder = NegaBinaryBPEncoder<Type, uint32_t>;
    use_negabinary = true;

    for(int i=0; i<n_vars; i++){
        // std::string rdir_prefix = rdata_file_prefix + "zone_" + id_str + "_" + varlist[i];
        std::string rdir_prefix = rdata_file_prefix + varlist[i];
        std::string metadata_file = rdir_prefix + "_refactored_data/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored_data/level_" + std::to_string(i) + ".bin";
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
            std::cout << "index = " << index << std::endl;
            refactor.refactor(buffer.data(), dims_masked, target_level, num_bitplanes);            
        } 
        else{
            refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes);            
        }
    }
}

template <class T>
static T compute_vr(const std::vector<T>& vec){
    T min = vec[0];
    T max = vec[0];
    for(int i=0; i<vec.size(); i++){
        if(vec[i] < min) min = vec[i];
        if(vec[i] > max) max = vec[i];
    }
    return max - min;
}

template<class Type>
void test_singleZone(int zone_id){

    std::string id_str = std::to_string(zone_id);

    size_t num_elements = 0;
    auto pressure_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_Density.dat").c_str(), num_elements);
    auto velocityX_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityZ.dat").c_str(), num_elements);

    std::vector<std::vector<Type>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    using Decomposer = MGARDHierarchicalDecomposer<Type>;
    using Interleaver = DirectInterleaver<Type>;
    // using Encoder = PerBitBPEncoder<Type, uint32_t>;
    using Compressor = AdaptiveLevelCompressor;
    using ErrorEstimator = MaxErrorEstimatorHB<Type>;
    using SizeInterpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<Type>>;
    using Retriever = ConcatLevelFileRetriever;

    // using Decomposer = MGARDOrthoganalDecomposer<Type>;
    // using ErrorEstimator = MaxErrorEstimatorOB<Type>;
    // using SizeInterpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorOB<Type>>;
    using Encoder = NegaBinaryBPEncoder<Type, uint32_t>;

    std::vector<MDR::ComposedReconstructor<Type, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>> reconstructors;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    size_t num_valid_data = 0;
    auto mask = MGARD::readfile<unsigned char>(mask_file.c_str(), num_valid_data);

    for(int i=0; i<n_vars; i++){
        // std::string rdir_prefix = rdata_file_prefix + "zone_" + id_str + "_" + varlist[i];
        std::string rdir_prefix = rdata_file_prefix + varlist[i];
        std::string metadata_file = rdir_prefix + "_refactored_data/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored_data/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto decomposer = Decomposer();
        auto interleaver = Interleaver();
        auto encoder = Encoder();
        auto compressor = Compressor(64);
        auto estimator = ErrorEstimator();
        auto interpreter = SizeInterpreter(estimator);
        auto retriever = Retriever(metadata_file, files);
        reconstructors.push_back(generateReconstructor<Type>(decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever));
        reconstructors.back().load_metadata();
    }
    
    std::vector<std::vector<double>> reconstructed_vars(n_vars, std::vector<double>(num_elements));
    std::vector<double> value_range(n_vars);
    for(int i=0; i<n_vars; i++){
        value_range[i] = compute_vr(vars_vec[i]);
        std::cout << "value_range = " << value_range[i] << std::endl;
    }    
    double rel_eb = 0.1;
    while(rel_eb > 1e-5){
        std::cout << "rel_eb = " << rel_eb << std::endl;
        for(int i=0; i<n_vars; i++){
            auto reconstructed_data = reconstructors[i].progressive_reconstruct(rel_eb*value_range[i], -1);
            size_t size = reconstructors[i].get_retrieved_size();
            if(i < 3){
                // reconstruct with mask
                int index = 0;
                for(int j=0; j<num_elements; j++){
                    if(mask[j]){
                        reconstructed_vars[i][j] = reconstructed_data[index ++];
                    }
                    else reconstructed_vars[i][j] = 0;
                }
                std::cout << varlist[i] << " bitrate = " <<  size * 1.0 / num_valid_data << std::endl;
            }
            else{
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(Type));
                std::cout << varlist[i] << " bitrate = " <<  size * 1.0 / num_elements << std::endl;
            }
            MGARD::print_statistics(vars_vec[i].data(), reconstructed_vars[i].data(), num_elements);
            // MGARD::writefile("reconstructed.dat", reconstructed_data, num_elements);
        }        
        rel_eb /= 2;
    }
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

template<class Type>
void refactor_singleZone_SZ3(int zone_id){
    std::string id_str = std::to_string(zone_id);

    size_t num_elements = 0;
    auto pressure_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_Density.dat").c_str(), num_elements);
    auto velocityX_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityZ.dat").c_str(), num_elements);
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
    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    std::vector<std::vector<Type>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<uint32_t> dims_masked;
    dims_masked.push_back(num_valid_data);
    std::vector<Type> buffer(num_valid_data);

    std::vector<double> value_range(n_vars);
    for(int i=0; i<n_vars; i++){
        value_range[i] = compute_vr(vars_vec[i]);
        std::cout << "value_range = " << value_range[i] << std::endl;
    }    
    std::vector<double> rel_ebs;
    const int num_snapshot = 8;
    double eb = 1.0;
    for(int i=0; i<num_snapshot; i++){
        eb /= 10;
        rel_ebs.push_back(eb);
    }
    for(int i=0; i<n_vars; i++){
        std::string rdir_prefix = rdata_file_prefix + varlist[i];
        if(i < 3){
            // use masked refactoring for vx vy vz
            int index = 0;
            for(int j=0; j<num_elements; j++){
                if(mask[j]){
                    buffer[index ++] = vars_vec[i][j];
                }
            }
            for(int j=0; j<num_snapshot; j++){
                std::string filename = rdir_prefix + "_refactored_data/SZ3_eb_" + std::to_string(j) + ".bin";
                size_t compressed_size = 0;
                auto compressed_data = SZ3_compress(num_valid_data, buffer.data(), rel_ebs[j]*value_range[i], compressed_size);
                MGARD::writefile(filename.c_str(), compressed_data, compressed_size);
                free(compressed_data);
            }
            std::cout << "index = " << index << std::endl;
        } 
        else{
            for(int j=0; j<num_snapshot; j++){
                std::string filename = rdir_prefix + "_refactored_data/SZ3_eb_" + std::to_string(j) + ".bin";
                size_t compressed_size = 0;
                auto compressed_data = SZ3_compress(num_elements, vars_vec[i].data(), rel_ebs[j]*value_range[i], compressed_size);
                MGARD::writefile(filename.c_str(), compressed_data, compressed_size);
                free(compressed_data);
            }
        }
    }
}

inline int find_index(double target_rel_eb, double& rel_eb){
    int i = 0;
    while(target_rel_eb < rel_eb){
        i ++;
        rel_eb /= 10;
    }
    return i;
}

template<class Type>
void test_singleZone_SZ3(int zone_id){

    std::string id_str = std::to_string(zone_id);

    size_t num_elements = 0;
    auto pressure_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_Density.dat").c_str(), num_elements);
    auto velocityX_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityZ.dat").c_str(), num_elements);

    std::vector<std::vector<Type>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};

    std::string mask_file = rdata_file_prefix + "mask.bin";
    size_t num_valid_data = 0;
    auto mask = MGARD::readfile<unsigned char>(mask_file.c_str(), num_valid_data);
    
    std::vector<std::vector<double>> reconstructed_vars(n_vars, std::vector<double>(num_elements));
    std::vector<double> value_range(n_vars);
    for(int i=0; i<n_vars; i++){
        value_range[i] = compute_vr(vars_vec[i]);
        std::cout << "value_range = " << value_range[i] << std::endl;
    }    
    double rel_eb = 0.1;
    std::vector<size_t> total_retrieved_sizes(n_vars, 0);
    std::vector<int> current_ind(n_vars, -1);

    Type * reconstructed_data = (Type *) malloc(num_elements * sizeof(Type));
    while(rel_eb > 1e-5){
        std::cout << "rel_eb = " << rel_eb << std::endl;
        for(int i=0; i<n_vars; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i];
            double file_eb = 0.1;
            auto file_ind = find_index(rel_eb, file_eb);
            std::cout << "Requested tolerance = " << rel_eb * value_range[i] << ", expected tolerance = " << file_eb * value_range[i] << "\n"; 
            std::string filename = rdir_prefix + "_refactored_data/SZ3_eb_" + std::to_string(file_ind) + ".bin";
            size_t n = 0;
            auto cmpData = MGARD::readfile<char>(filename.c_str(), n);
            SZ3_decompress(cmpData.data(), n, reconstructed_data);
            if(file_ind > current_ind[i]){
                current_ind[i] = file_ind;
                total_retrieved_sizes[i] += n;
            }
            if(i < 3){
                // reconstruct with mask
                int index = 0;
                for(int j=0; j<num_elements; j++){
                    if(mask[j]){
                        reconstructed_vars[i][j] = reconstructed_data[index ++];
                    }
                    else reconstructed_vars[i][j] = 0;
                }
                std::cout << varlist[i] << " bitrate = " <<  total_retrieved_sizes[i] * 1.0 / num_valid_data << std::endl;
            }
            else{
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(Type));
                std::cout << varlist[i] << " bitrate = " <<  total_retrieved_sizes[i] * 1.0 / num_elements << std::endl;
            }
            MGARD::print_statistics(vars_vec[i].data(), reconstructed_vars[i].data(), num_elements);
            // MGARD::writefile("reconstructed.dat", reconstructed_data, num_elements);
        }        
        rel_eb /= 2;
    }
    free(reconstructed_data);
}

template<class Type>
void refactor_singleZone_SZ3_delta(int zone_id){
    std::string id_str = std::to_string(zone_id);

    size_t num_elements = 0;
    auto pressure_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_Density.dat").c_str(), num_elements);
    auto velocityX_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityZ.dat").c_str(), num_elements);
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
    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    std::vector<std::vector<Type>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<uint32_t> dims_masked;
    dims_masked.push_back(num_valid_data);
    std::vector<Type> buffer(num_valid_data);

    std::vector<double> value_range(n_vars);
    for(int i=0; i<n_vars; i++){
        value_range[i] = compute_vr(vars_vec[i]);
        std::cout << "value_range = " << value_range[i] << std::endl;
    }    
    std::vector<double> rel_ebs;
    const int num_snapshot = 8;
    double eb = 1.0;
    for(int i=0; i<num_snapshot; i++){
        eb /= 10;
        rel_ebs.push_back(eb);
    }
    for(int i=0; i<n_vars; i++){
        std::string rdir_prefix = rdata_file_prefix + varlist[i];
        if(i < 3){
            // use masked refactoring for vx vy vz
            int index = 0;
            for(int j=0; j<num_elements; j++){
                if(mask[j]){
                    buffer[index ++] = vars_vec[i][j];
                }
            }
            std::vector<Type> data_buffer(buffer);
            std::vector<Type> dec_data_buffer(buffer);
            for(int j=0; j<num_snapshot; j++){
                std::string filename = rdir_prefix + "_refactored_data/SZ3_delta_eb_" + std::to_string(j) + ".bin";
                size_t compressed_size = 0;
                auto compressed_data = SZ3_compress(num_valid_data, data_buffer.data(), rel_ebs[j]*value_range[i], compressed_size);
                MGARD::writefile(filename.c_str(), compressed_data, compressed_size);
                SZ3_decompress(compressed_data, compressed_size, dec_data_buffer.data());
                for(int i=0; i<num_valid_data; i++){
                    data_buffer[i] = data_buffer[i] - dec_data_buffer[i];
                }
                free(compressed_data);
            }
            std::cout << "index = " << index << std::endl;
        } 
        else{
            std::vector<Type> data_buffer(vars_vec[i]);
            std::vector<Type> dec_data_buffer(vars_vec[i]);
            for(int j=0; j<num_snapshot; j++){
                std::string filename = rdir_prefix + "_refactored_data/SZ3_delta_eb_" + std::to_string(j) + ".bin";
                size_t compressed_size = 0;
                auto compressed_data = SZ3_compress(num_elements, data_buffer.data(), rel_ebs[j]*value_range[i], compressed_size);
                MGARD::writefile(filename.c_str(), compressed_data, compressed_size);
                SZ3_decompress(compressed_data, compressed_size, dec_data_buffer.data());
                for(int i=0; i<num_elements; i++){
                    data_buffer[i] = data_buffer[i] - dec_data_buffer[i];
                }
                free(compressed_data);
            }
        }
    }
}

template<class Type>
void test_singleZone_SZ3_delta(int zone_id){

    std::string id_str = std::to_string(zone_id);

    size_t num_elements = 0;
    auto pressure_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_Density.dat").c_str(), num_elements);
    auto velocityX_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityZ.dat").c_str(), num_elements);

    std::vector<std::vector<Type>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};

    std::string mask_file = rdata_file_prefix + "mask.bin";
    size_t num_valid_data = 0;
    auto mask = MGARD::readfile<unsigned char>(mask_file.c_str(), num_valid_data);
    
    std::vector<std::vector<double>> reconstructed_vars(n_vars, std::vector<double>(num_elements));
    std::vector<double> value_range(n_vars);
    for(int i=0; i<n_vars; i++){
        value_range[i] = compute_vr(vars_vec[i]);
        std::cout << "value_range = " << value_range[i] << std::endl;
    }    
    double rel_eb = 0.1;
    std::vector<size_t> total_retrieved_sizes(n_vars, 0);
    std::vector<int> current_ind(n_vars, -1);

    Type * reconstructed_data = (Type *) malloc(num_elements * sizeof(Type));
    while(rel_eb > 1e-5){
        std::cout << "rel_eb = " << rel_eb << std::endl;
        for(int i=0; i<n_vars; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i];
            double file_eb = 0.1;
            auto file_ind = find_index(rel_eb, file_eb);
            std::cout << "Requested tolerance = " << rel_eb * value_range[i] << ", expected tolerance = " << file_eb * value_range[i] << "\n"; 
            if(file_ind > current_ind[i]){
                for(int j=current_ind[i]+1; j<=file_ind; j++){
                    std::string filename = rdir_prefix + "_refactored_data/SZ3_delta_eb_" + std::to_string(j) + ".bin";
                    size_t n = 0;
                    auto cmpData = MGARD::readfile<char>(filename.c_str(), n);
                    total_retrieved_sizes[i] += n;
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
            if(i < 3){
                std::cout << varlist[i] << " bitrate = " <<  total_retrieved_sizes[i] * 1.0 / num_valid_data << std::endl;
            }
            else{
                std::cout << varlist[i] << " bitrate = " <<  total_retrieved_sizes[i] * 1.0 / num_elements << std::endl;
            }
            MGARD::print_statistics(vars_vec[i].data(), reconstructed_vars[i].data(), num_elements);
            // MGARD::writefile("reconstructed.dat", reconstructed_data, num_elements);
        }        
        rel_eb /= 2;
    }
    free(reconstructed_data);
}

}
#endif