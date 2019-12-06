#include <iostream>
#include <fstream>
#include <armadillo>
#include <omp.h>
#include <stdio.h>
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



    task_d(Ex,Cycles);


    return 0;
} // end main
