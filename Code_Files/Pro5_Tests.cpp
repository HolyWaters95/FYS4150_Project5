#include <iostream>
#include <armadillo>
#include <omp.h>
#include "time.h"
#include <random>
#include <chrono>

#include "Pro5_Functions.h"

using namespace std;
using namespace arma;

void test_stability(){
    cout << "Running stability test with no savings" << endl;
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
        transaction(i1,i2,M2);
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
        transaction_savings(i1,i2,0.5,M2);
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

