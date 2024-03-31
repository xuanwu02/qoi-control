#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include "utils.hpp"
#include "Synthesizer4GE.hpp"

using namespace std;

int main(int argc, char** argv){

    using T = double;
    int arg_pos = 1;
    int id = atoi(argv[arg_pos++]);
    bool refactor = false;
    if(argc > 2){
        refactor = true;
    }
    if(refactor){
        std::cout << "start refactoring\n";
        // MDR::refactor_singleZone<T>(id);
        // MDR::refactor_singleZone_SZ3<T>(id);
        MDR::refactor_singleZone_SZ3_delta<T>(id);
    }
    else{
        std::cout << "start recontructing\n";
        // MDR::test_singleZone<T>(id);
        // MDR::test_singleZone_SZ3<T>(id);
        MDR::test_singleZone_SZ3_delta<T>(id);
    }

    return 0;

}