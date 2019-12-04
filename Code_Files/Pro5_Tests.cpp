#include <iostream>
#include <armadillo>
#include <omp.h>
#include "time.h"
#include <random>
#include <chrono>

#include "Pro5_Functions.h"

using namespace std;
using namespace arma;

void test_sampling(){
    cout << "Running sampling index test" << endl;
    cout << "----------------------" << endl;
    vec M = vec(10,fill::randu);
    mat c = mat(M.n_elem,M.n_elem,fill::zeros);

    int fail = 0;
    for (int i = 0;i<100;i++){
        vector<int> index = Sampling_Rule(M,c,0,2);
        if (index[0] < 0 or index[0] >= M.n_elem or index[1] < 0 or index[1] >= M.n_elem or index[0] == index[1]){
            cout << "INDEX ERROR: \n ";
            cout << index[0] << "   " << index[1] << endl;
            fail = 1;
        }
    }
    if (fail == 0){
        cout << "Test passed, no index errors" << endl;
    }
    cout << "matrix for #transactions: \n ----------------------" << endl;
    cout << c << endl << sum(c) << endl;
    cout << "----------------------" << endl;
} //end test_sampling

void test_stability(){
    cout << "Running stability test with no savings" << endl;
    cout << "----------------------" << endl;

    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    int N = 10;
    int m0 = 100;
    vec M1;
    M1 = m0*vec(N,fill::ones);
    vec M2 = M1;
    mat c = mat(M1.n_elem,M1.n_elem,fill::zeros);

    int counter = 0;
    while (counter < 10){

    vector<int> index = Sampling_Rule(M2,c);
    transaction(index[0],index[1],0,M2);

    counter += 1;
    }
    cout << M1 << endl << endl << M2 << endl;

    if (sum(M1) == sum(M2)){
        cout << "Total money is conserved, and equal to " << sum(M1) << endl;
    }
    else{
        cout << "TOTAL MONEY IS NOT CONSERVED" << endl << "BEFORE: " << sum(M1) << endl << "AFTER: " << sum(M2) << endl;
    }
    cout << "----------------------" << endl;
}// end stability test

void test_stability_savings(){
    cout << "Running stability test with savings" << endl;
    cout << "----------------------" << endl;

    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    int N = 10;
    int m0 = 100;
    vec M1;
    M1 = m0*vec(N,fill::ones);
    vec M2 = M1;


    for (int i = 0; i < 10;i++){
        int i1 = static_cast<int>( (generate_canonical< double, 128 > (rng))*static_cast<double>(N/2) );
        int i2 = static_cast<int>( (generate_canonical< double, 128 > (rng))*static_cast<double>(N/2) + N/2 );
        transaction(i1,i2,0.5,M2);
    }
    cout << M1 << endl << endl << M2 << endl;

    if (sum(M1) == sum(M2)){
        cout << "Total money is conserved, and equal to " << sum(M1) << endl;
    }
    else{
        cout << "TOTAL MONEY IS NOT CONSERVED" << endl << "BEFORE: " << sum(M1) << endl << "AFTER: " << sum(M2) << endl;
    }
    cout << "----------------------" << endl;
} //end stability_savings test

void test_stability_taxes(){
    cout << "Running stability test with taxes" << endl;
    cout << "----------------------" << endl;

    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    int N = 10;
    int m0 = 100;
    vec M1;
    M1 = m0*vec(N,fill::ones);
    vec M2 = M1;


    for (int i = 0; i < 10;i++){
        int i1 = static_cast<int>( (generate_canonical< double, 128 > (rng))*static_cast<double>(N/2) );
        int i2 = static_cast<int>( (generate_canonical< double, 128 > (rng))*static_cast<double>(N/2) + N/2 );
        transaction_taxes(i1,i2,0.25,M2);
    }
    cout << M1 << endl << endl << M2 << endl;

    if (sum(M1) == sum(M2)){
        cout << "Total money is conserved, and equal to " << sum(M1) << endl;
    }
    else{
        cout << "TOTAL MONEY IS NOT CONSERVED" << endl << "BEFORE: " << sum(M1) << endl << "AFTER: " << sum(M2) << endl;
    }
    cout << "----------------------" << endl;
} //end stability_taxes test

