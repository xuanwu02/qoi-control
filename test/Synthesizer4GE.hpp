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
        auto decomposer = MDR::MGARDHierarchicalDecomposer<Type>();
        auto interleaver = MDR::DirectInterleaver<Type>();
        auto encoder = MDR::PerBitBPEncoder<Type, uint32_t>();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto collector = MDR::SquaredErrorCollector<Type>();
        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        auto refactor = generateRefactor<Type>(decomposer, interleaver, encoder, compressor, collector, writer);
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
    std::vector<MDR::ComposedReconstructor<Type, MGARDHierarchicalDecomposer<Type>, DirectInterleaver<Type>, PerBitBPEncoder<Type, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<Type>>, MaxErrorEstimatorHB<Type>, ConcatLevelFileRetriever>> reconstructors;
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
        auto decomposer = MGARDHierarchicalDecomposer<Type>();
        auto interleaver = DirectInterleaver<Type>();
        auto encoder = PerBitBPEncoder<Type, uint32_t>();
        auto compressor = AdaptiveLevelCompressor(64);
        auto estimator = MaxErrorEstimatorHB<Type>();
        auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<Type>>(estimator);
        auto retriever = ConcatLevelFileRetriever(metadata_file, files);
        reconstructors.push_back(generateReconstructor<Type>(decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever));
        reconstructors.back().load_metadata();
    }
    
    std::vector<std::vector<double>> reconstructed_vars(n_vars, std::vector<double>(num_elements));
    std::vector<double> var_eb_vec(n_vars, 1);
    for(int i=0; i<n_vars; i++){
        auto reconstructed_data = reconstructors[i].progressive_reconstruct(var_eb_vec[i], -1);
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
            memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(Type));
        }
        MGARD::print_statistics(vars_vec[i].data(), reconstructed_vars[i].data(), num_elements);
        // MGARD::writefile("reconstructed.dat", reconstructed_data, num_elements);
    }

}

}
#endif