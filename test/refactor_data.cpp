#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include "utils.hpp"
#include "Synthesizer4GE.hpp"

using namespace MDR;

int main(int argc, char** argv){

    using T = double;
    int arg_pos = 1;
    std::string data = argv[arg_pos++];

    if(data == "GE"){
        refactor_GE<T>();
        refactor_GE_SZ3<T>();
        refactor_GE_SZ3_delta<T>();
    }
    else{
        uint32_t n1 = atoi(argv[arg_pos++]);
        uint32_t n2 = atoi(argv[arg_pos++]);
        uint32_t n3 = atoi(argv[arg_pos++]);
        if(data == "S3D"){
            refactor_S3D<T>(n1, n2, n3);
            refactor_S3D_SZ3<T>(n1, n2, n3);
            refactor_S3D_SZ3_delta<T>(n1, n2, n3);
        }
        else if(data == "nyx" || data == "hurricane"){
            refactor_velocities_3d<T>(data, n1, n2, n3);
        }
    }
    
    return 0;

}