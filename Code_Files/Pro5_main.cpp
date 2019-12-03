#include <iostream>
#include <armadillo>
#include <omp.h>
#include "time.h"
#include <random>
#include <chrono>

#include "Pro5_Functions.h"

using namespace std;
using namespace arma;

int main(){

    string Tests;
    cout << "Do you want to run tests? y/n \n";
    cin >> Tests;
    if (Tests == "y"){
        test_sampling();
        test_stability();
        test_stability_savings();
        test_stability_taxes();
    }


    int Ex = 0;
    int Cycles = 0;
    cout << "Number of experiment runs? \n";
    cin >> Ex;
    cout << "Number of MC Cycles? \n";
    cin >> Cycles;

    Financial_analysis(Ex,Cycles,"Test_sampling","Test_smpl_median",0,1);

    //Financial_analysis(Ex,Cycles,"Tax_test","Tax_Median",1.);
    //savings(Ex,Cycles);

    return 0;
} // end main
